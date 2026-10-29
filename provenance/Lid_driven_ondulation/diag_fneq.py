# -*- coding: utf-8 -*-
"""
diag_fneq.py
============
INDEPENDENT estimation of the viscous dissipation from the non-equilibrium
moments of the LBM distribution, for the MRT collision operator of lbm_wavy.

Goal: provide D*_fneq to compare with D*_nabla_u (tensor S by finite
differences, computed in diagnostics.py). Both must use the same
normalisation and the same fluid mask.

Conventions (identical to lbm_wavy.py):
  - MRT collision: m = M f ; m* = m - S_diag (m - m_eq) ; f* = M_INV m*
  - shear moments 7,8 relaxed at omega = 1/tau
  - nu = c_s^2 (tau - 1/2), c_s^2 = 1/3 ; mu = rho * nu

Computation chain:
  f^neq = f - f^eq,  f^eq = M_INV m_eq(rho, u_lat)
  Pi^neq_ab = sum_q f^neq_q e_qa e_qb
  S_lat = -Pi^neq / (2 rho c_s^2 tau)   (MRT Chapman-Enskog coefficient)
  S* = S_lat * N / U_ref
  D*_fneq = sum_fluid 2 S*:S* (1/N)^2     (= same formula as D*_nabla_u)

The transfer coefficient is -2 rho c_s^2 TAU (not -2 rho nu = -2 rho c_s^2
(tau-1/2)). Verified empirically by node-by-node fitting (src/_probe_fneq)
on the flat cavity Re=100, N=64: slope -0.4355 vs -2 rho c_s^2 tau = -0.4333
(0.5%), whereas -2 rho nu = -0.100. The 'full' and 'dev' (trace removed)
versions are returned, as well as an interior mode (ring near the walls excluded).
"""

from __future__ import annotations

import numpy as np

from lbm_wavy import M, M_INV
from geometry import EX, EY

_EX2 = (EX.astype(np.float64) ** 2)
_EY2 = (EY.astype(np.float64) ** 2)
_EXY = (EX.astype(np.float64) * EY.astype(np.float64))


def equilibrium_populations(rho, ux, uy):
    """f^eq = M_INV m_eq for a (rho, ux, uy) field in lattice units."""
    r = rho
    u2 = ux * ux + uy * uy
    meq = np.empty((9,) + r.shape, dtype=np.float64)
    meq[0] = r
    meq[1] = -2.0 * r + 3.0 * r * u2
    meq[2] = r - 3.0 * r * u2
    meq[3] = r * ux
    meq[4] = -r * ux
    meq[5] = r * uy
    meq[6] = -r * uy
    meq[7] = r * (ux * ux - uy * uy)
    meq[8] = r * ux * uy
    return np.tensordot(M_INV, meq, axes=([1], [0]))     # (9, N, N)


def fneq_stress(f, rho, ux, uy):
    """Return non-equilibrium (Pi_xx, Pi_yy, Pi_xy), lattice units."""
    feq = equilibrium_populations(rho, ux, uy)
    fneq = f - feq
    Pi_xx = np.einsum("qij,q->ij", fneq, _EX2)
    Pi_yy = np.einsum("qij,q->ij", fneq, _EY2)
    Pi_xy = np.einsum("qij,q->ij", fneq, _EXY)
    return Pi_xx, Pi_yy, Pi_xy


def _interior_mask(fluid, width):
    """Fluid mask trimmed by a ring of width `width` (in nodes)."""
    if width <= 0:
        return fluid.astype(bool)
    m = fluid.astype(bool).copy()
    N = m.shape[0]
    for w in range(1, width + 1):
        m[w:N - w, w:N - w] &= fluid[w:N - w, w:N - w].astype(bool)
        m[w:N - w, :w] &= False
        m[w:N - w, N - w:] &= False
        m[:w, :] &= False
        m[N - w:, :] &= False
    return m


def fneq_dissipation(solver, interior_width=0):
    """D*_fneq and associated diagnostics for a converged WavyMRT solver.

    Returns a dict with D_fneq_full, D_fneq_dev (trace removed),
    D_nabla_u (reference, same mask), ratio, and S* fields.
    """
    geom = solver.geom
    N = geom["N"]
    d = geom["dx"]
    fluid = geom["fluid"].astype(bool)
    mask = _interior_mask(geom["fluid"], interior_width)

    rho = solver.rho                       # lattice (approx 1)
    ux_lat = solver.ux * solver.U_ref      # solver.ux is non-dimensionalised by U_ref
    uy_lat = solver.uy * solver.U_ref
    f = solver.f_final

    Pi_xx, Pi_yy, Pi_xy = fneq_stress(f, rho, ux_lat, uy_lat)

    # MRT Chapman-Enskog coefficient: Pi^neq = -(2 rho c_s^2 tau) S_lat.
    # (tau, not tau-1/2; see header + src/_probe_fneq.)
    coeff = 2.0 * rho * (1.0 / 3.0) * solver.tau
    S_xx = -Pi_xx / coeff * N / solver.U_ref
    S_yy = -Pi_yy / coeff * N / solver.U_ref
    S_xy = -Pi_xy / coeff * N / solver.U_ref

    # full
    D_full = float(np.sum(2.0 * (S_xx ** 2 + S_yy ** 2 + 2.0 * S_xy ** 2)[mask]) * d * d)
    # deviatoric (trace removed)
    tr = S_xx + S_yy
    sxx = S_xx - 0.5 * tr
    syy = S_yy - 0.5 * tr
    D_dev = float(np.sum(2.0 * (sxx ** 2 + syy ** 2 + 2.0 * S_xy ** 2)[mask]) * d * d)
    # trace (control: should be small in incompressible flow)
    D_trace = float(np.sum(2.0 * (0.5 * tr ** 2 + 0.5 * tr ** 2)[mask]) * d * d)

    # recall D*_nabla_u on the SAME mask
    from diagnostics import _gradients
    U = np.where(fluid, solver.ux, 0.0)
    V = np.where(fluid, solver.uy, 0.0)
    dudx, dudy, dvdx, dvdy, _, _ = _gradients(U, V, geom, solver.U_lid_norm)
    S11, S22, S12 = dudx, dvdy, 0.5 * (dudy + dvdx)
    D_nab = float(np.sum(2.0 * (S11 ** 2 + S22 ** 2 + 2.0 * S12 ** 2)[mask]) * d * d)

    return {
        "D_fneq_full": D_full,
        "D_fneq_dev": D_dev,
        "D_trace": D_trace,
        "D_nabla_u": D_nab,
        "ratio_full": D_full / D_nab if D_nab > 0 else np.nan,
        "ratio_dev": D_dev / D_nab if D_nab > 0 else np.nan,
        "interior_width": interior_width,
        "n_mask": int(mask.sum()),
        "ce_coeff": float(np.median(coeff[fluid])) if fluid.any() else np.nan,
    }


if __name__ == "__main__":
    import sys
    print("diag_fneq: f^neq diagnostic module (use test_fneq.py)")
