# -*- coding: utf-8 -*-
"""
06_neural_operator/physics_observables.py

Phase 6 (official) - physics observables over the 30-case main grid plus the
4 numerical-control diagnostics, mapped over (n, Re).

Frozen conventions (see 05_dataset/PROVENANCE.md):
    * dimensionless velocity     u_t  = u_lattice / U_lid
    * reference length unit a = 1 (L_c = 2 a = 2)   nu_t = 2 / Re
    * streamfunction            grad^2 psi = -omega,  psi = 0 on the boundary
    * vorticity                 omega = dv/dx - du/dy   (clockwise -> omega < 0)
    * vortex filtering          |psi_v| > eta * |psi_primary|,
                                eta in {0.5%, 1%, 2%, 5%}
    * lid stress                sigma = -p I + 2 nu_t S,  p = (rho-1)/(3 U_lid^2)
    * lid_power (CSV)           power INTO the fluid = +int u_lid.(sigma n) ds
                                (the energy balance P_lid = D is expected;
                                the user's formula -int u.sigma n ds gives the
                                power ON the lid, i.e. -lid_power)
    * force convention          F_x = -int_lid (sigma n)_x ds  (force on the lid)

Outputs (English labels, dpi >= 300):
    results/physics_observables.csv
    results/topology_map.png
    results/main_vortex_maps.png
    results/dissipation_lid.png
    results/energy_balance.png
    results/symmetry.png
    results/streamfunction_montage.png
"""

from __future__ import annotations

import csv
import os
import sys
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve
from scipy.ndimage import minimum_filter, maximum_filter, convolve

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "01_geometry"))
sys.path.insert(0, str(ROOT / "02_LBM"))

from superellipse import Superellipse
from figure_style import save_fig, style_colorbar

DATA_DIR = ROOT / "05_dataset" / "data"
DIAG_DIR = ROOT / "05_dataset" / "diagnostics"
OUT_DIR = ROOT / "06_neural_operator" / "results"

ETA_DEFAULT = 0.01
ETA_SWEEP = (0.005, 0.01, 0.02, 0.05)
RE_GRID = (100.0, 400.0, 1000.0)
N_GRID_MAIN = (1.25, 1.5, 1.75, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0)


# ----------------------------------------------------------------------
#  Grid handling
# ----------------------------------------------------------------------
def align_grid(X: np.ndarray, Y: np.ndarray, *fields):
    """Reorder so that axis0 = x (increasing), axis1 = y (increasing)."""
    fields = list(fields)
    if abs(X[1, 0] - X[0, 0]) < abs(Y[1, 0] - Y[0, 0]):
        X = X.T
        Y = Y.T
        fields = [np.asarray(f).T for f in fields]
    if X[1, 0] < X[0, 0]:
        X = X[::-1, :]
        fields = [f[::-1, :] for f in fields]
    if Y[0, 1] < Y[0, 0]:
        Y = Y[:, ::-1]
        fields = [f[:, ::-1] for f in fields]
    return X, Y, fields


def deep_mask(mask: np.ndarray) -> np.ndarray:
    """Fluid nodes whose 4 direct neighbours are fluid too (interior size)."""
    m = mask > 0
    return m[:-2, 1:-1] & m[2:, 1:-1] & m[1:-1, :-2] & m[1:-1, 2:]


def loose_interior(fluid: np.ndarray, min_neigh: int = 6) -> np.ndarray:
    """Fluid nodes with >= min_neigh fluid neighbours among the 8."""
    f = fluid > 0
    k = np.ones((3, 3), dtype=float)
    k[1, 1] = 0.0
    cnt = convolve(f.astype(float), k, mode="constant", cval=0.0)
    return f & (cnt >= min_neigh)


def embed_full(interior: np.ndarray, N: int) -> np.ndarray:
    """Embed an (N-2)x(N-2) interior mask into a full N x N mask."""
    full = np.zeros((N, N), dtype=bool)
    full[1:-1, 1:-1] = interior > 0
    return full


def fd_fields(ux: np.ndarray, uy: np.ndarray, dx: float, dy: float,
              deep_full: np.ndarray):
    """Omega and strain components on deep-interior nodes (NaN elsewhere)."""
    N = ux.shape[0]
    duxdx = np.gradient(ux, dx, axis=0)
    duxdy = np.gradient(ux, dy, axis=1)
    duydx = np.gradient(uy, dx, axis=0)
    duydy = np.gradient(uy, dy, axis=1)
    omega = duydx - duxdy
    sxx = duxdx
    syy = duydy
    sxy = 0.5 * (duxdy + duydx)
    for f in (omega, sxx, syy, sxy):
        f[~deep_full] = np.nan
    return {"omega": omega, "sxx": sxx, "syy": syy, "sxy": sxy}


def fill_omega(omega: np.ndarray) -> np.ndarray:
    """Replace NaN omega (near-wall) with the mean of finite 8-neighbours."""
    ones = np.ones((3, 3), dtype=float)
    ok = ~np.isnan(omega)
    num = convolve(np.where(ok, omega, 0.0), ones, mode="constant", cval=0.0)
    den = convolve(ok.astype(float), ones, mode="constant", cval=0.0)
    om = omega.copy()
    bad = ~ok
    ratio = np.zeros_like(num)
    np.divide(num, den, out=ratio, where=den > 0)
    om[bad] = np.where(den[bad] > 0, ratio[bad], np.nan)
    return om


def solve_streamfunction(omega: np.ndarray, mask_fluid, dx: float, dy: float) -> np.ndarray:
    """grad^2 psi = -omega on fluid nodes, psi = 0 on solid nodes.

    Discrete Laplacian (unit cell):  4 psi - sum(neighbours) = omega * h^2.
    Links are only present between two fluid nodes; out-of-range = Dirichlet 0.
    """
    N = omega.shape[0]
    n = N * N
    fluid = (mask_fluid > 0).ravel()
    h = 0.5 * (dx + dy)

    r, c, v = [], [], []
    for di, dj in ((0, 1), (1, 0)):
        Isrc, Jsrc = np.meshgrid(np.arange(N), np.arange(N), indexing="ij")
        Inb, Jnb = Isrc + di, Jsrc + dj
        ok = (Inb >= 0) & (Inb < N) & (Jnb >= 0) & (Jnb < N)
        i_src = (Isrc[ok] * N + Jsrc[ok]).ravel()
        i_nb = (Inb[ok] * N + Jnb[ok]).ravel()
        both = fluid[i_src] & fluid[i_nb]
        a, b = i_src[both], i_nb[both]
        r.append(np.concatenate([a, b]))
        c.append(np.concatenate([b, a]))
        v.append(np.full(2 * both.sum(), -1.0))

    r.append(np.arange(n))
    c.append(np.arange(n))
    v.append(np.where(fluid, 4.0, 1.0))

    A = sparse.coo_matrix(
        (np.concatenate(v), (np.concatenate(r), np.concatenate(c))),
        shape=(n, n),
    ).tocsr()

    rhs = np.zeros(n)
    rhs[fluid] = h * h * np.nan_to_num(omega.ravel()[fluid])
    sol = spsolve(A, rhs)
    return sol.reshape(N, N)


def find_extrema(psi: np.ndarray, emask: np.ndarray, X, Y, omega):
    """Local extrema of psi (8-neighbourhood) restricted to `emask` (full size)."""
    mx = maximum_filter(psi, size=3)
    mn = minimum_filter(psi, size=3)
    is_ext = (psi == mx) | (psi == mn)
    is_ext &= emask
    cands = []
    for i, j in zip(*np.where(is_ext)):
        cands.append({
            "i": int(i), "j": int(j),
            "x": float(X[i, j]), "y": float(Y[i, j]),
            "psi": float(psi[i, j]),
            "omega": float(omega[i, j]) if np.isfinite(omega[i, j]) else np.nan,
        })
    return cands


def nearest_nodes(x, y, X, Y, dx, dy, mask):
    """Nearest valid fluid node (i, j) for a set of boundary points."""
    i = np.rint((np.atleast_1d(x) - X[0, 0]) / dx).astype(int)
    j = np.rint((np.atleast_1d(y) - Y[0, 0]) / dy).astype(int)
    N = X.shape[0]
    out_i = np.full(i.size, -1)
    out_j = np.full(i.size, -1)
    for k in range(i.size):
        bestd = np.inf
        for di in (-1, 0, 1):
            for dj in (-1, 0, 1):
                ii, jj = i[k] + di, j[k] + dj
                if 0 <= ii < N and 0 <= jj < N and mask[ii, jj] > 0:
                    d = abs(ii - i[k]) + abs(jj - j[k])
                    if d < bestd:
                        bestd = d
                        out_i[k], out_j[k] = ii, jj
    return out_i, out_j


def sample_near(field: np.ndarray, mask: np.ndarray, x, y, X, Y, dx, dy):
    """Finite field value at (x,y) from the nearest valid fluid grid node."""
    i = np.rint((np.atleast_1d(x) - X[0, 0]) / dx).astype(int)
    j = np.rint((np.atleast_1d(y) - Y[0, 0]) / dy).astype(int)
    N = field.shape[0]
    out = np.full(i.size, np.nan)
    for k in range(i.size):
        best = None
        bestd = np.inf
        for di in (-1, 0, 1):
            for dj in (-1, 0, 1):
                ii, jj = i[k] + di, j[k] + dj
                if 0 <= ii < N and 0 <= jj < N and mask[ii, jj] > 0:
                    if np.isfinite(field[ii, jj]):
                        d = abs(ii - i[k]) + abs(jj - j[k])
                        if d < bestd:
                            bestd = d
                            best = float(field[ii, jj])
        if best is not None:
            out[k] = best
    return out


def mirror_symmetry(ux: np.ndarray, uy: np.ndarray, fluid: np.ndarray):
    """Mirror-residual about x = 0 (axis0): u -> -u, v -> v if flow were
    mirror-symmetric. The tangential lid breaks this symmetry (chiral flow),
    so the residual is expected to be O(1) by physics, not noise."""
    f = fluid > 0
    uxm = ux[::-1, :]
    uym = uy[::-1, :]
    valid = f & f[::-1, :]
    e_u = np.linalg.norm(ux[valid] - uxm[valid]) / np.linalg.norm(ux[valid])
    e_v = np.linalg.norm(uy[valid] + uym[valid]) / np.linalg.norm(uy[valid])
    return e_u, e_v


# ----------------------------------------------------------------------
#  Core observable computation
# ----------------------------------------------------------------------
def compute_observables(ux, uy, rho, fluid, X, Y, U_lid, Re, n, has_pressure=True):
    N = ux.shape[0]
    dx = (X[-1, 0] - X[0, 0]) / (N - 1)
    dy = (Y[0, -1] - Y[0, 0]) / (N - 1)
    dA = dx * dy
    nu_t = 2.0 / Re

    # NOTE: ux/uy are stored already dimensionless (in units of U_lid).
    ux_t = np.asarray(ux, dtype=np.float64)
    uy_t = np.asarray(uy, dtype=np.float64)

    deep_full = embed_full(deep_mask(fluid), N)
    g = fd_fields(ux_t, uy_t, dx, dy, deep_full)

    psi = solve_streamfunction(np.nan_to_num(g["omega"]), fluid, dx, dy)

    valid = ~np.isnan(g["omega"])
    twoSS = (
        2.0 * g["sxx"] * g["sxx"]
        + 2.0 * g["syy"] * g["syy"]
        + 4.0 * g["sxy"] * g["sxy"]
    )
    om2 = g["omega"] * g["omega"]
    diss_strain = 2.0 * nu_t * np.nansum(twoSS[valid]) * dA
    diss_omega = nu_t * np.nansum(om2[valid]) * dA
    diss_consistency = diss_strain / diss_omega if diss_omega > 0 else np.nan

    ff = fluid > 0
    kin = 0.5 * np.nansum((ux_t * ux_t + uy_t * uy_t)[ff]) * dA
    umax = float(np.nanmax(np.hypot(ux_t[ff], uy_t[ff])))

    emask = loose_interior(fluid)
    om_fill = fill_omega(g["omega"].copy())
    cands = find_extrema(psi, emask, X, Y, om_fill)

    if not cands:
        primary = None
        n_v = n_cw = n_ccw = {e: 0 for e in ETA_SWEEP}
    else:
        primary = max(cands, key=lambda c: abs(c["psi"]))
        pmag = abs(primary["psi"])
        n_v, n_cw, n_ccw = {}, {}, {}
        for e in ETA_SWEEP:
            thr = e * pmag
            cw = [c for c in cands if c["psi"] < 0 and abs(c["psi"]) > thr]
            ccw = [c for c in cands if c["psi"] >= 0 and abs(c["psi"]) > thr]
            n_cw[e] = len(cw)
            n_ccw[e] = len(ccw)
            n_v[e] = len(cw) + len(ccw)

    lid_power = np.nan
    f_lid_x = np.nan
    f_lid_y = np.nan
    if has_pressure:
        gd = Superellipse(a=1.0, b=1.0, n=n)
        bpts = gd.boundary_points(4096)
        lm = bpts[:, 1] > 1e-12
        lpts = bpts[lm]
        if lpts.shape[0] > 1:
            ds = np.hypot(np.diff(lpts[:, 0]), np.diff(lpts[:, 1]))
            ds = np.concatenate([ds, [0.0]])
            lnx, lny = gd.normal(lpts[:, 0], lpts[:, 1])
            tx, ty = lny, -lnx
            if np.mean(tx) < 0.0:
                tx, ty = -tx, -ty
            # --- wall-shear reconstruction: du_t/dn at the wall ----------
            gi, gj = nearest_nodes(lpts[:, 0], lpts[:, 1], X, Y, dx, dy, fluid)
            valid = gi >= 0
            if valid.any():
                xn = X[gi[valid], gj[valid]]
                yn = Y[gi[valid], gj[valid]]
                delta = np.hypot(lpts[valid, 0] - xn, lpts[valid, 1] - yn)
                u_f = ux_t[gi[valid], gj[valid]]
                v_f = uy_t[gi[valid], gj[valid]]
                u_t_f = u_f * tx[valid] + v_f * ty[valid]
                ok = (delta > 0) & np.isfinite(u_t_f)
                dudn = (1.0 - u_t_f[ok]) / delta[ok]
                lid_power = 2.0 * nu_t * np.sum(dudn * ds[valid][ok])
            # --- pressure strain for the lid force -------------------------
            p_t = (rho - 1.0) / (3.0 * U_lid * U_lid)
            ps = sample_near(p_t, fluid, lpts[:, 0], lpts[:, 1], X, Y, dx, dy)
            Sxxs = sample_near(g["sxx"], fluid, lpts[:, 0], lpts[:, 1], X, Y, dx, dy)
            Syys = sample_near(g["syy"], fluid, lpts[:, 0], lpts[:, 1], X, Y, dx, dy)
            Sxys = sample_near(g["sxy"], fluid, lpts[:, 0], lpts[:, 1], X, Y, dx, dy)
            okf = ~(np.isnan(ps) | np.isnan(Sxxs) | np.isnan(Syys) | np.isnan(Sxys))
            if okf.any():
                Sn_x = Sxxs[okf] * lnx[okf] + Sxys[okf] * lny[okf]
                Sn_y = Sxys[okf] * lnx[okf] + Syys[okf] * lny[okf]
                f_lid_x = np.sum(ps[okf] * lnx[okf] * ds[okf]) - 2.0 * nu_t * np.sum(Sn_x * ds[okf])
                f_lid_y = np.sum(ps[okf] * lny[okf] * ds[okf]) - 2.0 * nu_t * np.sum(Sn_y * ds[okf])

    e_sym_u, e_sym_v = mirror_symmetry(ux_t, uy_t, fluid)
    eps_e = np.abs(lid_power - diss_strain) / diss_strain if diss_strain > 0 else np.nan

    return {
        "psi": psi, "cands": cands, "primary": primary,
        "diss_strain": diss_strain, "diss_omega": diss_omega,
        "diss_consistency": diss_consistency, "kin": kin, "umax": umax,
        "lid_power": lid_power, "eps_e": eps_e,
        "f_lid_x": f_lid_x, "f_lid_y": f_lid_y,
        "e_sym_u": e_sym_u, "e_sym_v": e_sym_v,
        "n_v": n_v, "n_cw": n_cw, "n_ccw": n_ccw,
        "X": X, "Y": Y, "fluid": fluid,
        "area": np.count_nonzero(ff) * dA,
    }


# ----------------------------------------------------------------------
#  Case loaders
# ----------------------------------------------------------------------
CSV_FIELDS = [
    "case_id", "tag", "n", "Re", "N", "nrefs",
    "area", "umax",
    "psi_primary", "x_primary", "y_primary", "omega_primary",
    "n_vortices", "n_cw", "n_ccw",
    "n_v_eta_0p005", "n_v_eta_0p01", "n_v_eta_0p02", "n_v_eta_0p05",
    "dissipation", "diss_omega", "diss_consistency", "kinetic_energy",
    "lid_power", "energy_balance_error", "force_lid_x", "force_lid_y",
    "e_sym_u", "e_sym_v",
    "psi_primary_fem", "x_primary_fem", "y_primary_fem", "dx_primary",
    "iters",
]


def _base(z, N, area, r, tag, nrefs=np.nan, case_id=""):
    p = r["primary"]
    d = {
        "case_id": str(case_id), "tag": tag,
        "n": float(z["n"]), "Re": float(z["Re"]), "N": N, "nrefs": nrefs,
        "area": area, "umax": r["umax"],
        "psi_primary": p["psi"] if p else np.nan,
        "x_primary": p["x"] if p else np.nan,
        "y_primary": p["y"] if p else np.nan,
        "omega_primary": p["omega"] if p else np.nan,
        "n_vortices": r["n_v"][ETA_DEFAULT], "n_cw": r["n_cw"][ETA_DEFAULT],
        "n_ccw": r["n_ccw"][ETA_DEFAULT],
        "n_v_eta_0p005": r["n_v"][0.005], "n_v_eta_0p01": r["n_v"][0.01],
        "n_v_eta_0p02": r["n_v"][0.02], "n_v_eta_0p05": r["n_v"][0.05],
        "dissipation": r["diss_strain"], "diss_omega": r["diss_omega"],
        "diss_consistency": r["diss_consistency"], "kinetic_energy": r["kin"],
        "lid_power": r["lid_power"], "energy_balance_error": r["eps_e"],
        "force_lid_x": r["f_lid_x"], "force_lid_y": r["f_lid_y"],
        "e_sym_u": r["e_sym_u"], "e_sym_v": r["e_sym_v"],
        "iters": int(z["iters_lbm"]) if "iters_lbm" in z.files else np.nan,
    }
    return d


def main_case_row(path: Path) -> dict:
    z = np.load(path, allow_pickle=True)
    n, Re, N = float(z["n"]), float(z["Re"]), int(z["N"])
    U_lid = float(z["U_lid"])
    ux, uy, rho, fluid = (np.asarray(z[k]) for k in ("ux", "uy", "rho", "fluid"))
    Xr, Yr = np.asarray(z["X"]), np.asarray(z["Y"])

    X, Y, [ux, uy, rho, fluid] = align_grid(Xr, Yr, ux, uy, rho, fluid)
    r = compute_observables(ux, uy, rho, fluid, X, Y, U_lid, Re, n)

    fx, fy = np.asarray(z["fem_ux"]), np.asarray(z["fem_uy"])
    Xf, Yf, [fx, fy] = align_grid(Xr, Yr, fx, fy)
    rf = compute_observables(fx, fy, None, fluid, Xf, Yf, 1.0, Re, n,
                             has_pressure=False)
    pf = rf["primary"]

    d = _base(z, N, r["area"], r, "lbm_N128",
              nrefs=int(z["fem_nrefs"]), case_id=str(z["case_id"]))
    d["psi_primary_fem"] = pf["psi"] if pf else np.nan
    d["x_primary_fem"] = pf["x"] if pf else np.nan
    d["y_primary_fem"] = pf["y"] if pf else np.nan
    d["dx_primary"] = (d["x_primary"] - d["x_primary_fem"]) if (d["psi_primary"] == d["psi_primary"] and pf) else np.nan
    d["_row"] = r
    d["_row_fem"] = rf
    return d


def diag_lbm_row(path: Path) -> dict:
    z = np.load(path, allow_pickle=True)
    n, Re, N = float(z["n"]), float(z["Re"]), int(z["N"])
    U_lid = float(z["nu_lbm"]) * Re / (N - 1)
    ux, uy, rho, fluid = (np.asarray(z[k]) for k in ("ux", "uy", "rho", "fluid"))
    X, Y = np.asarray(z["X"]), np.asarray(z["Y"])
    X, Y, [ux, uy, rho, fluid] = align_grid(X, Y, ux, uy, rho, fluid)
    r = compute_observables(ux, uy, rho, fluid, X, Y, U_lid, Re, n)
    d = _base(z, N, r["area"], r, "lbm_N256_refinement",
              case_id=f"n{n:g}_Re{Re:g}_N{N}")
    d["_row"] = r
    return d


def diag_fem_row(path: Path) -> dict:
    z = np.load(path, allow_pickle=True)
    n, Re, nrefs = float(z["n"]), float(z["Re"]), int(z["nrefs"])
    ux, uy = np.asarray(z["fem_ux_N128"]), np.asarray(z["fem_uy_N128"])
    fluid = np.asarray(z["fluid_N128"]).astype(np.uint8)
    X, Y = np.asarray(z["X_N128"]), np.asarray(z["Y_N128"])
    X, Y, [ux, uy, fluid] = align_grid(X, Y, ux, uy, fluid)
    r = compute_observables(ux, uy, None, fluid, X, Y, 1.0, Re, n,
                            has_pressure=False)
    p = r["primary"]
    d = _base(z, 128, r["area"], r, f"fem_nrefs{nrefs}_refinement",
              nrefs=nrefs, case_id=f"n{n:g}_Re{Re:g}_fem_nrefs{nrefs}")
    d["psi_primary"] = np.nan
    d["x_primary"] = np.nan
    d["y_primary"] = np.nan
    d["omega_primary"] = np.nan
    d["psi_primary_fem"] = p["psi"] if p else np.nan
    d["x_primary_fem"] = p["x"] if p else np.nan
    d["y_primary_fem"] = p["y"] if p else np.nan
    d["iters"] = np.nan
    d["_row"] = r
    return d


# ----------------------------------------------------------------------
#  CSV + console report
# ----------------------------------------------------------------------
def write_csv(rows):
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    path = OUT_DIR / "physics_observables.csv"
    tmp = OUT_DIR / "physics_observables.csv.tmp"
    with open(tmp, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=CSV_FIELDS, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, np.nan) for k in CSV_FIELDS})
    for _ in range(5):
        try:
            os.replace(tmp, path)
            return str(path)
        except OSError:
            pass
    os.replace(tmp, path)
    return str(path)


def report(rows, data_rows):
    def hdr():
        s = "        "
        for n in N_GRID_MAIN:
            s += f"{n:>7.2f}"
        return s

    print("\n  Nv (eta=1%%) matrix:  rows=Re, cols=n")
    print(hdr())
    for re in RE_GRID:
        s = f"Re={re:>5.0f} "
        for n in N_GRID_MAIN:
            s += f"{data_rows[(n, re)]['n_vortices']:>7d}"
        print(s)

    print("\n  CCW count (eta=1%%):  rows=Re, cols=n")
    print(hdr())
    for re in RE_GRID:
        s = f"Re={re:>5.0f} "
        for n in N_GRID_MAIN:
            s += f"{data_rows[(n, re)]['n_ccw']:>7d}"
        print(s)

    print("\n  Primary vortex: x_p / y_p / psi_p  (main vortex)")
    print("        " + "  ".join(f"{n:>6.2f}" for n in N_GRID_MAIN))
    for re in RE_GRID:
        s = f"Re={re:>5.0f} "
        for n in N_GRID_MAIN:
            r = data_rows[(n, re)]
            s += f"({r['x_primary']:.2f},{r['y_primary']:.2f})"
        print(s)

    print("\n  Dissipation D  (2 nu_t int S:S, nu_t=2/Re):")
    print("        " + "  ".join(f"{n:>7.3f}" for n in N_GRID_MAIN))
    for re in RE_GRID:
        s = f"Re={re:>5.0f} "
        for n in N_GRID_MAIN:
            s += f"{float(data_rows[(n, re)]['dissipation']):>7.3f}"
        print(s)

    print("\n  Energy-balance relative error eps_E (%):")
    print("        " + "  ".join(f"{n:>7.2f}" for n in N_GRID_MAIN))
    for re in RE_GRID:
        s = f"Re={re:>5.0f} "
        for n in N_GRID_MAIN:
            v = float(data_rows[(n, re)]["energy_balance_error"])
            s += f"{100*v:>7.1f}" if v == v else "     na "
        print(s)

    print("\n  Symmetry residual E_sym,u (%% of ||u||) [chiral lid -> O(1)]:")
    print("        " + "  ".join(f"{n:>7.2f}" for n in N_GRID_MAIN))
    for re in RE_GRID:
        s = f"Re={re:>5.0f} "
        for n in N_GRID_MAIN:
            s += f"{100*float(data_rows[(n, re)]['e_sym_u']):>7.1f}"
        print(s)


# ----------------------------------------------------------------------
#  Figures (English, dpi >= 300)
# ----------------------------------------------------------------------
def grid_matrix(data_rows, key):
    A = np.zeros((len(RE_GRID), len(N_GRID_MAIN)))
    for i, re in enumerate(RE_GRID):
        for j, n in enumerate(N_GRID_MAIN):
            A[i, j] = float(data_rows[(n, re)][key])
    return A


def imap(ax, A, title, fmt="%.2f", label="", xl="n (-)", yl="Re (-)", cmap="Blues"):
    import matplotlib.pyplot as plt
    im = ax.imshow(A, aspect="auto", origin="lower", cmap=cmap,
                   extent=[N_GRID_MAIN[0] - 0.3, N_GRID_MAIN[-1] + 0.3,
                           RE_GRID[0] - 60, RE_GRID[-1] + 60])
    cb = ax.figure.colorbar(im, ax=ax)
    style_colorbar(cb, label=label)
    for i in range(A.shape[0]):
        for j in range(A.shape[1]):
            ax.text(N_GRID_MAIN[j], RE_GRID[i], fmt % A[i, j],
                    ha="center", va="center", fontsize=8, color="black")
    ax.set_xlabel(xl, fontsize=13)
    ax.set_ylabel(yl, fontsize=13)
    ax.set_title(title, fontsize=13, fontweight="bold")
    ax.tick_params(labelsize=11)
    return im


def make_figures(rows, data_rows, main_only):
    import matplotlib
    import matplotlib.pyplot as plt

    matplotlib.use("Agg")

    # --- topology map ------------------------------------------------
    A_nv = grid_matrix(data_rows, "n_vortices")
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    imap(axes[0], A_nv, "Nv: number of recirculation regions",
         fmt="%d", label="Nv (-)", cmap="YlGnBu")
    # eta sensitivity: per-Re band over n
    etas = np.asarray(ETA_SWEEP)
    for i, re in enumerate(RE_GRID):
        M = np.zeros((len(N_GRID_MAIN), len(etas)))
        for j, n in enumerate(N_GRID_MAIN):
            nv = data_rows[(n, re)]["_row"]["n_v"]
            for k, e in enumerate(etas):
                M[j, k] = nv[e]
        lo, med, hi = M.min(0), np.median(M, 0), M.max(0)
        axes[1].fill_between(etas, lo, hi,
                             alpha=0.12, color=f"C{i}")
        axes[1].plot(etas, med, marker="o" if i == 0 else (
            "s" if i == 1 else "^"), ms=6, label=f"Re={re:.0f} (-)")
        axes[1].set_xlabel("threshold eta (-)")
        axes[1].set_ylabel("Nv (-)")
        axes[1].set_title("Nv sensitivity to the strengths threshold eta")
        axes[1].grid(True, ls=":", lw=0.6, alpha=0.5)
        axes[1].legend(fontsize=10)
    fig.tight_layout()
    save_fig(fig, str(OUT_DIR / "topology_map.png"), dpi=300)
    print("  figure -> results/topology_map.png")

    # --- primary vortex maps -----------------------------------------
    keys = [("x_primary", "cx (x-coordinate of the primary vortex)", "%.2f"),
            ("y_primary", "cy (y-coordinate of the primary vortex)", "%.2f"),
            ("psi_primary", "|psi| of the primary vortex", "%.3f"),
            ("omega_primary", "vorticity at the primary vortex", "%.1f")]
    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    for ax, (k, t, f) in zip(axes.ravel(), keys):
        A = grid_matrix(data_rows, k)
        imap(ax, np.abs(A) if "psi" in k else A, t, fmt=f,
             label="(-)", cmap="YlGnBu")
    fig.suptitle("Primary vortex: position, strength, vorticity", fontsize=15)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    save_fig(fig, str(OUT_DIR / "main_vortex_maps.png"), dpi=300)
    print("  figure -> results/main_vortex_maps.png")

    # --- dissipation / lid -------------------------------------------
    keys = [("dissipation", "Dissipation D = 2 nu_t int S:S", "%.3f"),
            ("kinetic_energy", "Kinetic energy K = 0.5 int |u|^2", "%.3f"),
            ("lid_power", "Lid power input P_lid", "%.3f"),
            ("force_lid_x", "Lid force F_x", "%.3f")]
    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    for ax, (k, t, f) in zip(axes.ravel(), keys):
        A = grid_matrix(data_rows, k)
        imap(ax, A, t, fmt=f, label="(-)", cmap="YlGnBu")
    fig.suptitle("Global balances across (n, Re)", fontsize=15)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    save_fig(fig, str(OUT_DIR / "dissipation_lid.png"), dpi=300)
    print("  figure -> results/dissipation_lid.png")

    # --- energy balance ----------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    ax = axes[0]
    for i, re in enumerate(RE_GRID):
        _x = [data_rows[(n, re)]["dissipation"] for n in N_GRID_MAIN]
        _y = [data_rows[(n, re)]["lid_power"] for n in N_GRID_MAIN]
        ax.loglog(_x, _y, marker="o" if i == 0 else ("s" if i == 1 else "^"),
                  ls="", label=f"Re={re:.0f} (-)")
    vmin = min(min(data_rows[(n, re)]["dissipation"] for n in N_GRID_MAIN) for re in RE_GRID)
    vmax = max(max(data_rows[(n, re)]["lid_power"] for n in N_GRID_MAIN) for re in RE_GRID)
    ax.loglog([vmin, vmax], [vmin, vmax], "k--", lw=1.2, label="P_lid = D")
    ax.set_xlabel("dissipation D (-)")
    ax.set_ylabel("lid power P_lid (-)")
    ax.set_title("Energy balance: P_lid vs dissipation (steady state)")
    ax.grid(True, ls=":", lw=0.6, alpha=0.5)
    ax.legend(fontsize=10)
    ax = axes[1]
    for i, re in enumerate(RE_GRID):
        eps = [100 * float(data_rows[(n, re)]["energy_balance_error"])
               for n in N_GRID_MAIN]
        eps = [e for e in eps if e == e]
        ax.plot(np.full(len(eps), i) + np.random.uniform(-0.12, 0.12, len(eps)),
                eps, marker="o", ls="",
                label=f"Re={re:.0f} (-)")
    ax.set_xticks(range(len(RE_GRID)))
    ax.set_xticklabels([f"Re={r:.0f}" for r in RE_GRID])
    ax.set_ylabel("|P_lid - D| / D   (%)")
    ax.set_title("Energy-balance residual (discretization + boundary sampling)")
    ax.grid(True, ls=":", lw=0.6, alpha=0.5)
    ax.legend(fontsize=10)
    fig.tight_layout()
    save_fig(fig, str(OUT_DIR / "energy_balance.png"), dpi=300)
    print("  figure -> results/energy_balance.png")

    # --- symmetry ----------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.5))
    imap(axes[0], grid_matrix(data_rows, "e_sym_u"),
         "Mirror residual E_sym,u (x -> -x)", fmt="%.2f", label="E_sym,u (-)",
         cmap="OrRd")
    imap(axes[1], grid_matrix(data_rows, "e_sym_v"),
         "Mirror residual E_sym,v (v -> -v)", fmt="%.2f", label="E_sym,v (-)",
         cmap="OrRd")
    fig.suptitle("Mirror symmetry check (tangential lid breaks reflection -> O(1))",
                 fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    save_fig(fig, str(OUT_DIR / "symmetry.png"), dpi=300)
    print("  figure -> results/symmetry.png")

    # --- streamfunction montage --------------------------------------
    fig, axes = plt.subplots(len(RE_GRID), len(N_GRID_MAIN),
                             figsize=(24, 7))
    cm = plt.get_cmap("Greys")
    for i, re in enumerate(RE_GRID):
        for j, n in enumerate(N_GRID_MAIN):
            r = data_rows[(n, re)]
            ax = axes[i, j]
            X, Y = r["_row"]["X"], r["_row"]["Y"]
            psi = r["_row"]["psi"]
            flu = r["_row"]["fluid"] > 0
            vmax = np.nanmax(np.abs(psi[flu]))
            ax.pcolormesh(X, Y, np.where(flu, psi, np.nan), cmap=cm,
                          vmin=-vmax, vmax=vmax, shading="auto")
            lv = np.linspace(-0.95 * vmax, 0.95 * vmax, 11)[1:-1]
            lv = lv[np.abs(lv) > 1e-9]
            ax.contour(X, Y, psi, levels=lv, colors="k", linewidths=0.35)
            if r["_row"]["primary"] is not None:
                p = r["_row"]["primary"]
                ax.plot(p["x"], p["y"], "xk", ms=6)
            for c in r["_row"]["cands"]:
                if abs(c["psi"]) > ETA_DEFAULT * abs(
                        r["_row"]["primary"]["psi"]) if r["_row"]["primary"] else False:
                    ax.plot(c["x"], c["y"], "ok", ms=2.5, mfc="none")
            ax.set_xticks([])
            ax.set_yticks([])
            if i == 0:
                ax.set_title(f"n={n:g}", fontsize=10)
            if j == 0:
                ax.set_ylabel(f"Re={re:.0f}", fontsize=10)
            ax.set_xlim(-1.05, 1.05)
            ax.set_ylim(-1.05, 1.05)
    fig.suptitle("Streamfunction ~psi (grayscale) with vortex markers "
                 "(x=primary, o=secondary, >1% of primary)",
                 fontsize=16)
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    save_fig(fig, str(OUT_DIR / "streamfunction_montage.png"), dpi=300)
    print("  figure -> results/streamfunction_montage.png")


# ----------------------------------------------------------------------
#  Main
# ----------------------------------------------------------------------
def build_rows():
    rows = []
    data_rows = {}
    for f in sorted(DATA_DIR.glob("se_*.npz")):
        r = main_case_row(f)
        rows.append(r)
        data_rows[(r["n"], r["Re"])] = r
    for f in sorted(DIAG_DIR.glob("diag_lbm_*.npz")):
        rows.append(diag_lbm_row(f))
    for f in sorted(DIAG_DIR.glob("diag_ref_*.npz")):
        rows.append(diag_fem_row(f))
    return rows, data_rows


def main():
    print("Building physics-observable rows (30 main + 4 diagnostics)...")
    rows, data_rows = build_rows()
    path = write_csv(rows)
    print(f"  CSV -> {path}")
    report(rows, data_rows)
    print("\nBuilding figures (dpi=300, English labels)...")
    make_figures(rows, data_rows, main_only=None)
    print("\nAll done.")


if __name__ == "__main__":
    main()