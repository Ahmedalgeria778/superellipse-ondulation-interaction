# -*- coding: utf-8 -*-
"""
lbm_wavy.py
===========
MRT-LBM D2Q9 for the lid-driven cavity with a wavy lower wall.

Strict reworking of the validated MRT solver (Paper 2/3, lbm_mrt_reference_square.py) :
  - same MRT collision operator (moment space, diagonal S),
    S = (0, 1.4, 1.4, 0, 1.2, 0, 1.2, 1/tau, 1/tau)
  - same PULL streaming + half-way bounce-back
  - same momentum injection at the lid (q = 7,8)

Modifications applied (plan, sections 1, 3, 5) :
  - solid/fluid mask for the wavy lower wall y_b(x)
  - curved boundary condition via interpolated bounce-back (Bouzidi et al. 2001,
    linear interpolation, see geometry.build_geometry for the fraction q)
  - classic (uniform) or corner-regularized lid (plan 3.2)
  - collision/streaming only on fluid nodes

LBM parameter scaling (plan, section 5.2) : a target tau is set per
Reynolds number and U_0,LB = (tau-1/2)/3 * Re / N is deduced from it, so that
Re = U_0,LB N / nu and nu = c_s^2 (tau - 1/2) with c_s^2 = 1/3.
"""

from __future__ import annotations

import os
import time
import numpy as np
from numba import njit, prange

from geometry import (
    EX, EY, OPP, build_geometry, lid_profile_regularized, derived_params,
)

W = np.array([4 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 36, 1 / 36, 1 / 36, 1 / 36],
             dtype=np.float64)

# M transformation matrix (9 moments) — identical to the reference solver
M = np.array([
    [1, 1, 1, 1, 1, 1, 1, 1, 1],
    [-4, -1, -1, -1, -1, 2, 2, 2, 2],
    [4, -2, -2, -2, -2, 1, 1, 1, 1],
    [0, 1, 0, -1, 0, 1, -1, -1, 1],
    [0, -2, 0, 2, 0, 1, -1, -1, 1],
    [0, 0, 1, 0, -1, 1, 1, -1, -1],
    [0, 0, -2, 0, 2, 1, 1, -1, -1],
    [0, 1, -1, 1, -1, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 1, -1, 1, -1]
], dtype=np.float64)
M_INV = np.linalg.inv(M)

# target tau per Reynolds number (plan, section 5.2) ; U_0,LB follows from it for any N
TARGET_TAU = {100: 0.65, 400: 0.60, 1000: 0.54}
DEFAULT_TAU = 0.60


@njit(parallel=True, cache=True)
def mrt_step_wavy(f, M, M_INV, S_diag, EX, EY, W, OPP, U_lid2d, t_idx,
                  fluid, bb_q, bb_i, bb_j, bb_qf, bb_low, bb_bi, bb_bj, N, n_bb):
    Q = 9
    rho = np.zeros((N, N))
    ux = np.zeros((N, N))
    uy = np.zeros((N, N))
    f_post = np.empty_like(f)

    # 1) Macros + MRT collision on fluid nodes only
    for i in prange(N):
        for j in range(N):
            if fluid[i, j] == 0:
                continue
            r = 0.0
            mx = 0.0
            my = 0.0
            for q in range(Q):
                fq = f[q, i, j]
                r += fq
                mx += EX[q] * fq
                my += EY[q] * fq
            rho[i, j] = max(r, 1e-10)
            ux_ = mx / rho[i, j]
            uy_ = my / rho[i, j]
            ux[i, j] = ux_
            uy[i, j] = uy_
            u2 = ux_ * ux_ + uy_ * uy_

            m_eq = np.zeros(Q)
            m_eq[0] = r
            m_eq[1] = -2.0 * r + 3.0 * r * u2
            m_eq[2] = r - 3.0 * r * u2
            m_eq[3] = r * ux_
            m_eq[4] = -r * ux_
            m_eq[5] = r * uy_
            m_eq[6] = -r * uy_
            m_eq[7] = r * (ux_ * ux_ - uy_ * uy_)
            m_eq[8] = r * ux_ * uy_

            m = np.zeros(Q)
            for q in range(Q):
                acc = 0.0
                for k in range(Q):
                    acc += M[q, k] * f[k, i, j]
                m[q] = acc
            m_star = np.zeros(Q)
            for q in range(Q):
                m_star[q] = m[q] - S_diag[q] * (m[q] - m_eq[q])
            for q in range(Q):
                acc = 0.0
                for k in range(Q):
                    acc += M_INV[q, k] * m_star[k]
                f_post[q, i, j] = acc

    # 2) PULL streaming + half-way bounce-back (side walls, lid)
    f_new = np.empty_like(f)
    for i in prange(N):
        for j in range(N):
            if fluid[i, j] == 0:
                for q in range(Q):
                    f_new[q, i, j] = f[q, i, j]
                continue
            for q in range(Q):
                ex = int(EX[q])
                ey = int(EY[q])
                prev_i = i - ey
                prev_j = j - ex
                if 0 <= prev_i < N and 0 <= prev_j < N and fluid[prev_i, prev_j] == 1:
                    f_new[q, i, j] = f_post[q, prev_i, prev_j]
                else:
                    f_new[q, i, j] = f_post[OPP[q], i, j]
                    if prev_i == N:  # moving lid: momentum injection
                        uw = U_lid2d[t_idx, j]
                        rho_w = rho[i, j]
                        if q == 7:
                            f_new[7, i, j] -= 6.0 * W[5] * rho_w * uw
                        elif q == 8:
                            f_new[8, i, j] += 6.0 * W[6] * rho_w * uw

    # 3) Wavy lower wall: interpolated bounce-back (Bouzidi, linear)
    for k in prange(n_bb):
        q = bb_q[k]
        i = bb_i[k]
        j = bb_j[k]
        qf = bb_qf[k]
        if bb_low[k] == 1:
            bi = bb_bi[k]
            bj = bb_bj[k]
            f_new[q, i, j] = (2.0 * qf * f_post[OPP[q], i, j]
                              + (1.0 - 2.0 * qf) * f_post[OPP[q], bi, bj])
        else:
            f_new[q, i, j] = ((1.0 / (2.0 * qf)) * f_post[OPP[q], i, j]
                              + ((2.0 * qf - 1.0) / (2.0 * qf)) * f_post[q, i, j])

    return f_new, rho, ux, uy


class WavyMRT:
    """D2Q9 MRT-LBM solver, square cavity, wavy lower wall."""

    def __init__(self, Re, N, h_over_H=0.0, Nw=1, lid='regularized', eps_lid=0.02,
                 flat_height=None, U_ref=None, target_tau=None, max_iter=200000,
                 tol=1e-9, verbose=True, log=None):
        self.Re = int(Re)
        self.N = int(N)
        self.h_over_H = float(h_over_H)
        self.Nw = int(Nw)
        self.lid = lid
        self.eps_lid = float(eps_lid)
        self.flat_height = flat_height
        self.max_iter = int(max_iter)
        self.tol = float(tol)
        self.verbose = verbose
        self.log = log

        # --- LBM scaling (plan 5.2) ---
        if target_tau is None:
            target_tau = TARGET_TAU.get(self.Re, DEFAULT_TAU)
        self.target_tau = float(target_tau)
        self.tau = float(target_tau)
        if U_ref is not None:
            self.U_ref = float(U_ref)
        else:
            self.U_ref = (self.tau - 0.5) / 3.0 * self.Re / self.N
        self.nu = self.U_ref * self.N / self.Re
        # consistency: nu = 1/3 (tau - 1/2)
        self.tau = 3.0 * self.nu + 0.5
        self.omega = 1.0 / self.tau
        self.Ma = np.sqrt(3.0) * self.U_ref

        self.S_diag = np.array(
            [0.0, 1.4, 1.4, 0.0, 1.2, 0.0, 1.2, self.omega, self.omega],
            dtype=np.float64)

        # --- geometry ---
        self.geom = build_geometry(self.N, self.h_over_H, self.Nw,
                                   flat_height=self.flat_height)
        self.fluid = self.geom["fluid"]
        self.solid = self.geom["solid"]

        # --- lid profile ---
        xc = self.geom["xc"]
        if lid == 'uniform':
            self.U_lid_norm = np.ones(self.N)
        elif lid == 'regularized':
            self.U_lid_norm = lid_profile_regularized(xc, self.eps_lid)
        elif lid == 'smoothstep_mean1':
            # protocol B: same mean velocity as the uniform lid.
            # the profile is renormalized by its discrete mean so that
            # mean(U_lid/U_ref) = 1 exactly on the simulation mesh.
            g_lid = lid_profile_regularized(xc, self.eps_lid)
            self.U_lid_norm = g_lid / g_lid.mean()
            self.lid_mean = float(self.U_lid_norm.mean())
        else:
            raise ValueError(f"unknown lid: {lid}")
        if not hasattr(self, "lid_mean"):
            self.lid_mean = float(np.mean(self.U_lid_norm))
        self.U_lid2d = (self.U_ref * self.U_lid_norm)[None, :]

        # --- initial state: fluid at rest (plan 3.3) ---
        self.f = np.ascontiguousarray(
            W[:, None, None] * np.ones((self.N, self.N))[None], dtype=np.float64)
        self.f_init = None
        self.thermalized = False

    # ------------------------------------------------------------------
    def _advance(self, f, n_steps, t_offset=0):
        """Advances n_steps steps; returns final (f, rho, ux, uy)."""
        g = self.geom
        period = self.U_lid2d.shape[0]  # = 1 (steady lid)
        for k in range(n_steps):
            f, rho, ux, uy = mrt_step_wavy(
                f, M, M_INV, self.S_diag, EX, EY, W, OPP, self.U_lid2d, 0,
                self.fluid, g["bb_q"], g["bb_i"], g["bb_j"], g["bb_qf"],
                g["bb_low"], g["bb_bi"], g["bb_bj"], self.N, g["n_bb"])
        return f, rho, ux, uy

    def run(self):
        N = self.N
        f = (self.f.copy() if self.f_init is None else self.f_init.copy())
        t0 = time.time()
        ux_prev = np.zeros((N, N))
        self.conv_hist = []

        # lid spin-up over t* = 5 (plan 3.3), smooth ramp
        ramp_iters = int(round(5.0 / self.U_ref / 1.0))
        ramp_iters = max(ramp_iters, 1)
        if self.verbose and self.log:
            self.log.write(f"   lid ramp: {ramp_iters} steps (t*=5)\n")
        self.U_lid2d = (self.U_lid_norm)[None, :]
        for k in range(ramp_iters):
            amp = (k + 1) / ramp_iters
            self.U_lid2d = (self.U_ref * amp * self.U_lid_norm)[None, :]
            f, rho, ux, uy = mrt_step_wavy(
                f, M, M_INV, self.S_diag, EX, EY, W, OPP, self.U_lid2d, 0,
                self.fluid, self.geom["bb_q"], self.geom["bb_i"], self.geom["bb_j"],
                self.geom["bb_qf"], self.geom["bb_low"], self.geom["bb_bi"],
                self.geom["bb_bj"], N, self.geom["n_bb"])
        self.U_lid2d = (self.U_ref * self.U_lid_norm)[None, :]

        # steady-state loop
        ok = 0
        mask = (self.fluid == 1)
        self.converged = False
        it = 0
        for it in range(1, self.max_iter + 1):
            f, rho, ux, uy = mrt_step_wavy(
                f, M, M_INV, self.S_diag, EX, EY, W, OPP, self.U_lid2d, 0,
                self.fluid, self.geom["bb_q"], self.geom["bb_i"], self.geom["bb_j"],
                self.geom["bb_qf"], self.geom["bb_low"], self.geom["bb_bi"],
                self.geom["bb_bj"], N, self.geom["n_bb"])
            if it % 2000 == 0:
                err = float(np.max(np.abs((ux - ux_prev)[mask])))
                self.conv_hist.append((it, err, time.time() - t0))
                if self.verbose and it % 20000 == 0 and self.log:
                    self.log.write(f"   it={it:7d} dU={err:.3e} t={time.time()-t0:.1f}s\n")
                    self.log.flush()
                if err < self.tol and it > 5000:
                    ok += 1
                    if ok >= 2:
                        self.converged = True
                        break
                else:
                    ok = 0
                ux_prev = ux.copy()

        self.iters = it
        self.last_residual = err if 'err' in locals() else float('inf')
        self.status = 'converged' if self.converged else 'max_iter'
        self.ux = ux / self.U_ref
        self.uy = uy / self.U_ref
        self.rho = rho
        self.ux = np.where(mask, self.ux, 0.0)
        self.uy = np.where(mask, self.uy, 0.0)
        self.rho = np.where(mask, self.rho, 1.0)
        self.f_final = f.copy()
        # fluid mass (relative fluctuation)
        self.mass = float(self.rho[mask].sum())
        self.mass_drift = abs(self.mass - mask.sum()) / mask.sum()
        if self.verbose and self.log:
            self.log.write(
                f"   Re={self.Re} N={N} h/H={self.h_over_H} Nw={self.Nw} "
                f"lid={self.lid} | U={self.U_ref:.5f} tau={self.tau:.4f} "
                f"Ma={self.Ma:.3f} | it={it} resid={self.last_residual:.2e} "
                f"mass_drift={self.mass_drift:.2e} | {time.time()-t0:.1f}s\n")
        return self


if __name__ == "__main__":
    import sys
    Re = int(sys.argv[1]) if len(sys.argv) > 1 else 100
    N = int(sys.argv[2]) if len(sys.argv) > 2 else 64
    h = float(sys.argv[3]) if len(sys.argv) > 3 else 0.0
    nw = int(sys.argv[4]) if len(sys.argv) > 4 else 1
    s = WavyMRT(Re, N, h, nw, lid='regularized', max_iter=80000, tol=1e-8,
                verbose=True, log=sys.stdout)
    s.run()
    print("resid =", s.last_residual, "mass_drift =", s.mass_drift)
