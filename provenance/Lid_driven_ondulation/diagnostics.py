# -*- coding: utf-8 -*-
"""
diagnostics.py
==============
Extraction of the plan quantities (sections 9 and 13) for each case :

  - topology: streamfunction psi, vorticity omega, detection of
    vortex centers (local extrema of psi), number of recirculations ;
  - kinetic energy  K  and K* = K / (L H U0^2) ;
  - mechanical power of the lid  P_lid = mu int U_lid (du/dy)|lid dx ;
  - viscous dissipation  D_nu = int 2 mu S:S dA ;
  - energy balance  P_lid vs D_nu (must satisfy P ~ D at steady state).

Derivatives are computed by finite differences on the fluid field, with
an explicit correction of the normal derivative at the walls (distance to the wall).
"""

from __future__ import annotations

import os as _os

import numpy as np
from scipy import ndimage

from geometry import wall_height, wall_slope

# --- isolated D*_fneq branch (out-of-equilibrium moments), DISABLED by default ---
# Explicit activation: LBM_DIAG_FNEQ=1 (interior width: LBM_DIAG_FNEQ_INTERIOR)
# The computation is done at run time (requires solver.f_final) ; it is absent
# from reanalyze, which only has the saved fields. By default, no
# behavior change for the current / upcoming campaign.
FNEQ_ENABLED = _os.environ.get("LBM_DIAG_FNEQ", "0") == "1"
FNEQ_INTERIOR = int(_os.environ.get("LBM_DIAG_FNEQ_INTERIOR", "0"))


def _gradients(U, V, geom, U_lid_norm):
    """Velocity gradient components, with wall correction."""
    N = geom["N"]
    d = geom["dx"]
    fluid = geom["fluid"].astype(bool)
    solid = ~fluid

    dudx = np.gradient(U, d, axis=1)
    dudy = np.gradient(U, d, axis=0)
    dvdx = np.gradient(V, d, axis=1)
    dvdy = np.gradient(V, d, axis=0)

    # --- lid (i = N-1), wall 0.5*d above the node ---
    for j in range(N):
        if fluid[N - 1, j]:
            dudy[N - 1, j] = (U_lid_norm[j] - U[N - 1, j]) / (0.5 * d)
            dvdy[N - 1, j] = (0.0 - V[N - 1, j]) / (0.5 * d)

    # --- side walls (j = 0 and j = N-1), wall at 0.5*d ---
    for i in range(N):
        if fluid[i, 0]:
            dudx[i, 0] = (0.0 - U[i, 0]) / (0.5 * d)
            dvdx[i, 0] = (0.0 - V[i, 0]) / (0.5 * d)
        if fluid[i, N - 1]:
            dudx[i, N - 1] = (0.0 - U[i, N - 1]) / (0.5 * d)
            dvdx[i, N - 1] = (0.0 - V[i, N - 1]) / (0.5 * d)

    # --- wavy lower wall: lowest fluid node per column ---
    yb = geom["yb"]
    for j in range(N):
        col = np.nonzero(fluid[:, j])[0]
        if col.size == 0:
            continue
        i0 = col[0]
        ybj = yb[0, j]
        dist = (i0 + 0.5) * d - ybj           # vertical distance to the wall
        if dist > 1e-12:
            dudy[i0, j] = (0.0 - U[i0, j]) / dist
            dvdy[i0, j] = (0.0 - V[i0, j]) / dist
    return dudx, dudy, dvdx, dvdy, fluid, solid


def streamfunction(U, geom):
    """psi such that u = d psi/dy, v = -d psi/dx, psi = 0 on the lower wall."""
    N = geom["N"]
    d = geom["dx"]
    fluid = geom["fluid"].astype(bool)
    psi = np.zeros((N, N), dtype=np.float64)
    for j in range(N):
        col = np.nonzero(fluid[:, j])[0]
        if col.size == 0:
            continue
        i0 = col[0]
        acc = 0.0
        psi[i0, j] = 0.0
        for i in range(i0 + 1, N):
            if not fluid[i, j]:
                break
            acc += 0.5 * (U[i, j] + U[i - 1, j]) * d
            psi[i, j] = acc
    return psi


def _subgrid_extremum(p, i, j):
    """Sub-cell refinement (separable parabola) of an extremum of p at (i,j).

    Returns (di, dj), offsets in index units, bounded by +-0.5. Allows
    comparing centers to a reference on a different grid (e.g. Ghia j/N)
    without the systematic half-cell bias of the (j+0.5)/N convention.
    """
    ni, nj = p.shape
    di = dj = 0.0
    if 1 <= i <= ni - 2:
        a, b, c = p[i - 1, j], p[i, j], p[i + 1, j]
        den = a - 2.0 * b + c
        if abs(den) > 1e-300:
            di = max(-0.5, min(0.5, 0.5 * (a - c) / den))
    if 1 <= j <= nj - 2:
        a, b, c = p[i, j - 1], p[i, j], p[i, j + 1]
        den = a - 2.0 * b + c
        if abs(den) > 1e-300:
            dj = max(-0.5, min(0.5, 0.5 * (a - c) / den))
    return di, dj


def detect_vortices(psi, geom, rel_thresh=1e-4, min_sep=5, margin=0.02):
    """Vortex centers = significant local extrema of psi.

    Filters applied:
      - amplitude |psi| >= rel_thresh * |psi|_main (relative amplitude criterion) ;
      - interior extrema: corner/wall mask, outside a `margin` band of
        the lid (y=1), the side walls (x=0,1) and the
        lower wall (y=yb(x)) -- removes corner artifacts ;
      - deduplication: minimum distance min_sep (in nodes) between centers ;
      - parabolic sub-cell refinement of each extremum.
    Returns a list of dicts {x, y, psi, sign, amp} sorted by decreasing |psi|.
    """
    N = geom["N"]
    d = geom["dx"]
    fluid = geom["fluid"].astype(bool)
    yb = geom["yb"]
    p = np.where(fluid, psi, 0.0)

    mx = ndimage.maximum_filter(p, size=min_sep)
    mn = ndimage.minimum_filter(p, size=min_sep)

    cand = (p == mx) | (p == mn)
    idx = np.argwhere(cand & fluid)
    vals = p[idx[:, 0], idx[:, 1]]

    pts = []
    for (i, j), v in zip(idx, vals):
        yb_j = float(yb[0, j])
        di, dj = _subgrid_extremum(p, int(i), int(j))
        x = (j + 0.5 + dj) * d
        y = (i + 0.5 + di) * d
        # corner/wall mask: minimum distance to any wall
        dist = min(1.0 - y, x, 1.0 - x, y - yb_j)
        if dist < margin:
            continue
        pts.append({"i": int(i), "j": int(j), "x": float(x), "y": float(y),
                    "psi": float(v), "sign": 1 if v > 0 else -1,
                    "amp": float(abs(v)), "wall_dist": float(dist)})
    if not pts:
        return []
    pts.sort(key=lambda r: -r["amp"])
    main_amp = pts[0]["amp"]
    thr = rel_thresh * main_amp

    out = []
    for r in pts:
        if r["amp"] < thr:
            continue
        if all((r["x"] - o["x"]) ** 2 + (r["y"] - o["y"]) ** 2
               >= (min_sep * d) ** 2 for o in out):
            out.append(r)
    return out


def compute_diagnostics(solver):
    """Diagnostics of a converged LBM solver (wrapper)."""
    out = compute_diagnostics_fields(
        solver.ux, solver.uy,
        _grid_geom(solver),
        solver.U_lid_norm,
        case_id=getattr(solver, "case_id", None),
        Re=solver.Re, N=solver.N,
        h_over_H=solver.h_over_H, Nw=solver.Nw, lid=solver.lid,
        eps_lid=solver.eps_lid, U_ref=solver.U_ref, tau=solver.tau,
        Ma=solver.Ma, mass_drift=solver.mass_drift,
        rho=solver.rho,
        last_residual=solver.last_residual)
    _maybe_fneq(solver, out)
    return out


def _maybe_fneq(solver, out):
    """Isolated branch (opt-in): adds D*_fneq and associated columns to `out`.

    Does nothing if FNEQ_ENABLED is false (behavior unchanged). On
    error, the exception is LOGGED explicitly (D_star_fneq_error) and
    D_star_fneq is set to NaN -- we do not silently hide a
    numerical problem, but neither do we fail the whole run.
    """
    if not FNEQ_ENABLED:
        return
    out["D_star_fneq"] = float("nan")
    out["D_star_fneq_dev"] = float("nan")
    out["D_star_fneq_ratio"] = float("nan")
    out["D_star_fneq_error"] = None
    try:
        from diag_fneq import fneq_dissipation
        r = fneq_dissipation(solver, interior_width=FNEQ_INTERIOR)
        out["D_star_fneq"] = float(r["D_fneq_full"])
        out["D_star_fneq_dev"] = float(r["D_fneq_dev"])
        out["D_star_fneq_ratio"] = float(r["ratio_full"])
        if not np.isfinite(out["D_star_fneq"]):
            raise FloatingPointError(
                f"D*_fneq non-finite<｜begin▁of▁sentence｜><{out['D_star_fneq']}>")
    except Exception as exc:                             # noqa: BLE001
        out["D_star_fneq_error"] = f"{type(exc).__name__}: {exc}"


def _grid_geom(solver):
    """Geometry as seen by the diagnostics (regular grid, also useful for FEM)."""
    N = solver.N
    dx = solver.geom["dx"]
    return {
        "N": N, "dx": dx, "fluid": solver.geom["fluid"],
        "solid": solver.geom["solid"],
        "xc": solver.geom["xc"], "yc": solver.geom["yc"],
        "X": solver.geom["X"], "Y": solver.geom["Y"],
        "yb": solver.geom["yb"],
    }


def compute_diagnostics_fields(ux, uy, geom, U_lid_norm,
                               case_id=None, Re=None, N=None,
                               h_over_H=None, Nw=None, lid=None, eps_lid=0.02,
                               U_ref=1.0, tau=None, Ma=None, rho=None,
                               ma_threshold=0.1,
                               mass_drift=None, last_residual=None):
    """Diagnostics from fields ux, uy (dimensionless with respect to U0) and geom."""
    N = geom["N"]
    d = geom["dx"]
    fluid = geom["fluid"].astype(bool)
    U = np.where(fluid, ux, 0.0)
    V = np.where(fluid, uy, 0.0)

    dudx, dudy, dvdx, dvdy, fluid, solid = _gradients(U, V, geom, U_lid_norm)
    omega = dvdx - dudy
    psi = streamfunction(U, geom)

    # Dimensionless fields u*, S*: u* = u/U0 (already the case in solver.ux).
    # Plan normalizations (sections 9.2, 9.3) :
    #   K*  = 0.5 int |u*|^2 dA*
    #   D*  = 2 int S*:S* dA*            (= D_nu / (mu U0^2))
    #   P*  = int U*_lid (du*/dy*)|lid dx*
    S11, S22, S12 = dudx, dvdy, 0.5 * (dudy + dvdx)
    Dstar = float(np.sum(2.0 * (S11 ** 2 + S22 ** 2 + 2.0 * S12 ** 2)[fluid]) * d * d)
    K = float(0.5 * np.sum((U ** 2 + V ** 2)[fluid]) * d * d)
    Z = float(np.sum((omega ** 2)[fluid]) * d * d)

    P = 0.0
    for j in range(N):
        if fluid[N - 1, j]:
            P += U_lid_norm[j] * dudy[N - 1, j]
    P = float(P * d)

    area = float(fluid.sum()) * d * d
    vortices = detect_vortices(psi, geom)
    main = vortices[0] if vortices else None

    # --- Mach / compressibility quality check (plan 3.3) --------------
    # U_ref is the lid velocity in lattice units; ux,uy are
    # dimensionless w.r.t. U_ref, so we reconstruct the lattice velocity.
    speed_lat = np.sqrt(U ** 2 + V ** 2) * U_ref
    Ma_max_local = float(np.sqrt(3.0) * speed_lat.max()) if speed_lat.size else 0.0
    max_drho = float(np.max(np.abs(rho - 1.0))) if rho is not None else None
    Ma_ok = bool(Ma_max_local <= ma_threshold)
    if not Ma_ok:
        import warnings
        warnings.warn(
            f"local Mach {Ma_max_local:.4f} > threshold {ma_threshold:.3f} for "
            f"{case_id} (N={N}): Mach regime exceeded.", RuntimeWarning)

    out = {
        "case": case_id,
        "Re": Re, "N": N,
        "h_over_H": h_over_H, "Nw": Nw, "lid": lid,
        "eps_lid": eps_lid,
        "U_ref": U_ref, "tau": tau, "Ma": Ma,
        "Ma_max_local": Ma_max_local, "max_drho": max_drho,
        "Ma_threshold": ma_threshold, "Ma_ok": Ma_ok,
        "params": geom.get("params"),
        "fluid_area": area,
        "mass_drift": mass_drift,
        "last_residual": last_residual,
        "K": K, "K_star": K,
        "P_star": P,
        "D_star": Dstar,
        "Z": Z,
        "balance": (P - Dstar) / Dstar if Dstar > 0 else np.nan,
        "n_vortices": len(vortices),
        "main_vortex": main,
        "vortices": vortices,
    }
    return out
