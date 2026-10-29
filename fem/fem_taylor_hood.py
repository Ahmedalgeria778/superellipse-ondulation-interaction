# -*- coding: utf-8 -*-
"""
fem_taylor_hood.py
==================
Taylor-Hood P2/P1 finite element solver (velocity/pressure) for the
steady Navier-Stokes equations in a lid-driven cavity with a corrugated
bottom wall (section 8).

  - triangular elements: P2 velocity (6 nodes), P1 pressure (3 nodes);
  - domain conforming to the wall y_b(x), structured mesh with boundary-layer
    refinement towards the wall and the lid;
  - Picard linearisation (initial Stokes step then iterations);
  - BCs: u = U0 g(x) on the lid (regularised profile, section 3.2),
         u = 0 elsewhere (no-slip);
  - projection onto the regular LBM grid for identical post-processing
    (diagnostics.compute_diagnostics_fields).

Reference: Taylor & Hood (1973); stable P2/P1 spacing (LBB).
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from geometry import wall_height, lid_profile_regularized, build_geometry
from diagnostics import compute_diagnostics_fields


# ----------------------------------------------------------------------
#  Quadrature on the reference triangle (Duffy transform)
# ----------------------------------------------------------------------
def _gl_01(n):
    xs, ws = np.polynomial.legendre.leggauss(n)
    return (xs + 1.0) / 2.0, ws / 2.0


def _tri_quadrature(n=6):
    g1, w1 = _gl_01(n)
    g2, w2 = _gl_01(n)
    pts, wts = [], []
    for a, wa in zip(g1, w1):
        for b, wb in zip(g2, w2):
            pts.append((a, (1.0 - a) * b))
            wts.append(wa * wb * (1.0 - a))
    return np.array(pts), np.array(wts)


QUAD_PTS, QUAD_WTS = _tri_quadrature(6)


def _shape_p2(x, y):
    """P2 shape functions and their (reference) gradients at point (x,y)."""
    l1, l2, l3 = 1.0 - x - y, x, y
    N = np.array([l1 * (2.0 * l1 - 1.0), l2 * (2.0 * l2 - 1.0),
                  l3 * (2.0 * l3 - 1.0), 4.0 * l1 * l2,
                  4.0 * l2 * l3, 4.0 * l3 * l1])
    dNx = np.array([-(4.0 * l1 - 1.0), (4.0 * l2 - 1.0), 0.0,
                    4.0 * (l1 - l2), 4.0 * l3, -4.0 * l3])
    dNy = np.array([-(4.0 * l1 - 1.0), 0.0, (4.0 * l3 - 1.0),
                    -4.0 * l2, 4.0 * l2, 4.0 * (l1 - l3)])
    return N, np.stack([dNx, dNy], axis=0)


def _shape_p1(x, y):
    l1, l2, l3 = 1.0 - x - y, x, y
    N = np.array([l1, l2, l3])
    g = np.array([[-1.0, 1.0, 0.0],
                  [-1.0, 0.0, 1.0]])
    return N, g


# ----------------------------------------------------------------------
#  Conforming mesh
# ----------------------------------------------------------------------
def wall_values(x, h_over_H, Nw, flat_height=None):
    x = np.asarray(x, dtype=np.float64)
    if flat_height is not None:
        return np.full_like(x, float(flat_height))
    return wall_height(x, h_over_H, Nw)


def generate_mesh(nx, ny, h_over_H, Nw, flat_height=None):
    """Conforming structured mesh: vertical columns x_k, y nodes
    graded towards the wall and the lid (cosine)."""
    verts = []
    for k in range(nx + 1):
        x = k / nx
        yb = wall_values(x, h_over_H, Nw, flat_height)
        for j in range(ny + 1):
            s = (1.0 - np.cos(np.pi * j / ny)) / 2.0
            verts.append((x, yb + (1.0 - yb) * s))
    verts = np.asarray(verts)
    nv = (nx + 1) * (ny + 1)

    def v(k, j):
        return k * (ny + 1) + j

    mid = {}
    mids = []

    def _mid(a, b):
        key = frozenset((a, b))
        if key in mid:
            return mid[key]
        nid = nv + len(mids)
        mids.append(0.5 * (verts[a] + verts[b]))
        mid[key] = nid
        return nid

    elems = []
    for k in range(nx):
        for j in range(ny):
            a, b = v(k, j), v(k + 1, j)
            c, d = v(k, j + 1), v(k + 1, j + 1)
            elems.append([a, b, c, _mid(a, b), _mid(b, c), _mid(c, a)])
            elems.append([b, d, c, _mid(b, d), _mid(d, c), _mid(c, b)])
    nodes = np.vstack([verts, np.array(mids)])
    return {"nodes": nodes, "nv": nv, "elems": np.asarray(elems, dtype=np.int64),
            "nx": nx, "ny": ny, "h_over_H": h_over_H, "Nw": Nw,
            "flat_height": flat_height, "v": v}


def _box_bc(mesh, U0, eps):
    """Dirichlet values (ux, uy) per node; None for interior nodes."""
    nx, ny = mesh["nx"], mesh["ny"]
    nodes = mesh["nodes"]
    M = nodes.shape[0]
    ux = np.full(M, np.nan)
    uy = np.full(M, np.nan)
    for k in range(nx + 1):
        x = k / nx
        yb = wall_values(x, mesh["h_over_H"], mesh["Nw"], mesh["flat_height"])
        for j in range(ny + 1):
            nid = k * (ny + 1) + j
            y = nodes[nid, 1]
            on_lid = abs(y - 1.0) < 1e-9
            on_wall = (x < 1e-9 or x > 1.0 - 1e-9 or abs(y - yb) < 1e-9)
            if on_lid:
                ux[nid] = U0 * lid_profile_regularized(np.array([x]), eps)[0]
                uy[nid] = 0.0
            elif on_wall:
                ux[nid] = 0.0
                uy[nid] = 0.0
    return ux, uy, M


def _mid_bc(mesh, U0, eps):
    """BCs for midpoint (edge) nodes lying on the boundary."""
    nodes = mesh["nodes"]
    out = {}
    for e in mesh["elems"]:
        for q in range(3):
            a, b, m = int(e[q]), int(e[(q + 1) % 3]), int(e[q + 3])
            pa, pb = nodes[a], nodes[b]
            if abs(pa[1] - 1.0) < 1e-9 and abs(pb[1] - 1.0) < 1e-9:
                xm = 0.5 * (pa[0] + pb[0])
                out[m] = (U0 * lid_profile_regularized(np.array([xm]), eps)[0], 0.0)
            elif (abs(pa[0]) < 1e-9 and abs(pb[0]) < 1e-9) or \
                 (abs(pa[0] - 1.0) < 1e-9 and abs(pb[0] - 1.0) < 1e-9):
                out[m] = (0.0, 0.0)
            else:
                yba = wall_values(pa[0], mesh["h_over_H"], mesh["Nw"],
                                  mesh["flat_height"])
                ybb = wall_values(pb[0], mesh["h_over_H"], mesh["Nw"],
                                  mesh["flat_height"])
                if abs(pa[1] - yba) < 1e-9 and abs(pb[1] - ybb) < 1e-9:
                    out[m] = (0.0, 0.0)
    return out


# ----------------------------------------------------------------------
#  Assembly (linear Picard)
# ----------------------------------------------------------------------
def assemble(mesh, nu, U, M, P):
    """Assemble the linearised system; U: previous-iteration vector."""
    nodes = mesh["nodes"]
    elems = mesh["elems"]
    dof_total = 2 * M + P
    rows, cols, vals = [], [], []
    JT = np.zeros((2, 2))
    ux_p = U[:M] if U is not None else None
    uy_p = U[M:] if U is not None else None

    # quadrature arrays (precomputed)
    Qx_, Qy_ = QUAD_PTS[:, 0], QUAD_PTS[:, 1]
    NQ = np.empty((len(QUAD_PTS), 6))
    GQX = np.empty((len(QUAD_PTS), 6))
    GQY = np.empty((len(QUAD_PTS), 6))
    NpQ = np.empty((len(QUAD_PTS), 3))
    for q, (xq, yq) in enumerate(zip(Qx_, Qy_)):
        NQ[q], gr = _shape_p2(xq, yq)
        GQX[q] = gr[0]
        GQY[q] = gr[1]
        NpQ[q], _ = _shape_p1(xq, yq)

    for e in elems:
        ids = [int(x) for x in e]
        VA, VB, VC = nodes[ids[0]], nodes[ids[1]], nodes[ids[2]]
        JT[:, 0] = VB - VA
        JT[:, 1] = VC - VA
        det = JT[0, 0] * JT[1, 1] - JT[0, 1] * JT[1, 0]
        vol = abs(det) / 2.0
        # J = [[dxB,dxC],[dyB,dyC]] ; J^{-T} = [[dyC,-dyB],[-dxC,dxB]]/det
        JinvT = np.array([[JT[1, 1], -JT[1, 0]],
                          [-JT[0, 1], JT[0, 0]]]) / det
        px = JinvT[0, 0] * GQX + JinvT[0, 1] * GQY
        py = JinvT[1, 0] * GQX + JinvT[1, 1] * GQY
        wq = QUAD_WTS * abs(det)                # (Q,) ; sum(wq) = volume

        K = nu * (np.einsum("qi,qj,q->ij", px, px, wq)
                  + np.einsum("qi,qj,q->ij", py, py, wq))
        if U is not None:
            upxq = (ux_p[ids] @ NQ.T)          # (Q,)
            upyq = (uy_p[ids] @ NQ.T)
            UQ = upxq[:, None] * px + upyq[:, None] * py     # (Q,6) u_prev.grad Nj
            C = np.einsum("qi,qj,q->ij", NQ, UQ, wq)
        else:
            C = np.zeros((6, 6))
        Qx = -np.einsum("qa,qi,q->ia", NpQ, px, wq)
        Qy = -np.einsum("qa,qi,q->ia", NpQ, py, wq)

        Ael = K + C
        for i in range(6):
            gi = ids[i]
            for j in range(6):
                gj = ids[j]
                rows.append(gi); cols.append(gj); vals.append(Ael[i, j])
                rows.append(M + gi); cols.append(M + gj); vals.append(Ael[i, j])
            for a in range(3):
                ga = ids[a]
                rows.append(2 * M + ga); cols.append(gi); vals.append(-Qx[i, a])
                rows.append(2 * M + ga); cols.append(M + gi); vals.append(-Qy[i, a])
                rows.append(gi); cols.append(2 * M + ga); vals.append(Qx[i, a])
                rows.append(M + gi); cols.append(2 * M + ga); vals.append(Qy[i, a])

    A = sp.csr_matrix((vals, (rows, cols)), shape=(dof_total, dof_total))
    return A


def _apply_bc(A, rhs, ux_v, uy_v, mids, M):
    """Penalisation of Dirichlet dofs (vertex + midpoint)."""
    big = 1e12 * (np.abs(A.data).max() if A.data.size else 1.0)
    b_rows, b_cols, b_vals = [], [], []
    for nid in np.where(~np.isnan(ux_v))[0]:
        b_rows.append(nid); b_cols.append(nid); b_vals.append(big)
        rhs[nid] += big * ux_v[nid]
        dd = nid + M
        b_rows.append(dd); b_cols.append(dd); b_vals.append(big)
        rhs[dd] += big * uy_v[nid]
    for m, (ux0, uy0) in mids.items():
        b_rows.append(m); b_cols.append(m); b_vals.append(big)
        rhs[m] += big * ux0
        dd = m + M
        b_rows.append(dd); b_cols.append(dd); b_vals.append(big)
        rhs[dd] += big * uy0
    A = A + sp.coo_matrix((b_vals, (b_rows, b_cols)),
                          shape=A.shape).tocsr()
    A[2 * M, 2 * M] += big          # zero reference pressure
    return A, rhs


def solve_picard(mesh, Re, U0=1.0, eps=0.02, max_it=60, rtol=1e-10,
                 relax=0.85, verbose=False, log=None):
    """Solve steady NS P2/P1 by Picard."""
    nu = U0 / Re
    M = mesh["nodes"].shape[0]
    P = mesh["nv"]
    ux_v, uy_v, _ = _box_bc(mesh, U0, eps)
    mids = _mid_bc(mesh, U0, eps)

    # initialization: Stokes solution (curl = 0, no convection)
    U = np.zeros(2 * M + P)
    U[:M] = np.where(np.isnan(ux_v), 0.0, ux_v)
    U[M:2 * M] = np.where(np.isnan(uy_v), 0.0, uy_v)
    t0 = time.time()
    rel_hist = []
    Unew = U.copy()
    for it in range(1, max_it + 1):
        A = assemble(mesh, nu, U, M, P)
        rhs = np.zeros(2 * M + P)
        A, rhs = _apply_bc(A, rhs, ux_v, uy_v, mids, M)
        try:
            Unew = spsolve(A.tocsr(), rhs)
        except RuntimeError as exc:
            raise RuntimeError(f"spsolve failure it={it}: {exc}")
        diff = float(np.max(np.abs(Unew[: 2 * M] - U[: 2 * M])))
        nrm = float(np.max(np.abs(Unew[: 2 * M])))
        rel = diff / nrm if nrm > 0 else 1.0
        rel_hist.append(rel)
        U = relax * Unew + (1.0 - relax) * U
        if verbose and (it <= 3 or it % 10 == 0) and log:
            log.write(f"   fem it={it:3d} rel={rel:.2e}\n")
        if rel < rtol:
            break
    U = Unew.copy()
    return {"U": U, "M": M, "P": P, "iters": it, "rel": rel,
            "rel_hist": rel_hist, "elapsed": time.time() - t0,
            "nu": nu, "mesh": mesh}


# ----------------------------------------------------------------------
#  Interpolation onto the regular grid + diagnostics
# ----------------------------------------------------------------------
def _edge_mid_cache(mesh):
    """Dict (a,b) -> midpoint node id for all mesh edges."""
    cm = {}
    for e in mesh["elems"]:
        for q in range(3):
            a, b = int(e[q]), int(e[(q + 1) % 3])
            cm[frozenset((a, b))] = int(e[q + 3])
    return cm


def interp_to_grid(mesh, U, N, h_over_H, Nw, flat_height=None, mask=None):
    """Project the FEM solution onto the regular grid (cell centres) N x N."""
    nodes = mesh["nodes"]
    nx, ny = mesh["nx"], mesh["ny"]
    M = mesh["nodes"].shape[0]
    ux_p, uy_p = U[:M], U[M:]
    e2m = _edge_mid_cache(mesh)

    xc = (np.arange(N) + 0.5) / N
    yc = (np.arange(N) + 0.5) / N
    v = mesh["v"]
    kc = np.minimum((xc * nx).astype(int), nx - 1)

    # triangular cells: for each (k,j), the triangles (a,b,c)
    tris = {}
    for k in range(nx):
        for j in range(ny):
            tris[(k, j)] = [
                (v(k, j), v(k + 1, j), v(k, j + 1)),
                (v(k + 1, j), v(k + 1, j + 1), v(k, j + 1)),
            ]

    # mean y-ordinates of the layers for localization (per column k)
    ym = np.empty((ny + 1, nx))
    for k in range(nx):
        for j in range(ny + 1):
            ym[j, k] = 0.5 * (nodes[v(k, j), 1] + nodes[v(k + 1, j), 1])

    Ux = np.zeros((N, N))
    Uy = np.zeros((N, N))
    for i in range(N):
        y = yc[i]
        for j in range(N):
            if mask is not None and not mask[i, j]:
                continue
            k = kc[j]
            jl = np.searchsorted(ym[:, k], y, side="right") - 1
            if jl < 0:
                jl = 0
            if jl >= ny:
                jl = ny - 1
            found = False
            for (a, b, c) in tris[(k, jl)]:
                lam2, lam3 = _bary(nodes[a], nodes[b], nodes[c], xc[j], y)
                lam1 = 1.0 - lam2 - lam3
                if min(lam1, lam2, lam3) >= -1e-10:
                    Ux[i, j] = _p2_eval(nodes, e2m, ux_p, a, b, c, lam1, lam2, lam3)
                    Uy[i, j] = _p2_eval(nodes, e2m, uy_p, a, b, c, lam1, lam2, lam3)
                    found = True
                    break
            if not found:
                best = None
                for (a, b, c) in tris[(k, jl)]:
                    lam2, lam3 = _bary(nodes[a], nodes[b], nodes[c], xc[j], y)
                    lam1 = 1.0 - lam2 - lam3
                    score = -min(lam1, lam2, lam3)
                    if best is None or score < best[0]:
                        best = (score, a, b, c, lam1, lam2, lam3)
                _, a, b, c, lam1, lam2, lam3 = best
                Ux[i, j] = _p2_eval(nodes, e2m, ux_p, a, b, c, lam1, lam2, lam3)
                Uy[i, j] = _p2_eval(nodes, e2m, uy_p, a, b, c, lam1, lam2, lam3)
    return Ux, Uy


def _bary(A, B, C, x, y):
    dxB, dyB = B[0] - A[0], B[1] - A[1]
    dxC, dyC = C[0] - A[0], C[1] - A[1]
    det = dxB * dyC - dyB * dxC
    if abs(det) < 1e-30:
        return 1 / 3.0, 1 / 3.0
    ex, ey = x - A[0], y - A[1]
    lam2 = (ex * dyC - ey * dxC) / det
    lam3 = (dxB * ey - dyB * ex) / det
    return lam2, lam3


def _p2_eval(nodes, e2m, vals, a, b, c, l1, l2, l3):
    Nv, _ = _shape_p2(l2, l3)
    mids = []
    for (p, q) in ((a, b), (b, c), (c, a)):
        m = e2m.get(frozenset((p, q)), -1)
        mids.append(vals[m] if m >= 0 else 0.5 * (vals[p] + vals[q]))
    return Nv @ (vals[a], vals[b], vals[c], mids[0], mids[1], mids[2])


def solve_to_diagnostics(case, nx, ny, N=128, U0=1.0, eps=0.02,
                         max_it=60, rtol=1e-10, relax=0.85,
                         verbose=False, log=None, outdir="results"):
    """Full driver: mesh -> Picard -> projection -> diagnostics."""
    h, nw, Re = case["h_over_H"], case["Nw"], case["Re"]
    fh = case.get("flat_height")
    mesh = generate_mesh(nx, ny, h, nw, flat_height=fh)
    t0 = time.time()
    sol = solve_picard(mesh, Re, U0, eps, max_it, rtol, relax,
                       verbose=verbose, log=log)
    geom = build_geometry(N, h, nw, flat_height=fh)
    mask = geom["fluid"].astype(bool)
    Ux, Uy = interp_to_grid(mesh, sol["U"][: 2 * sol["M"]], N, h, nw,
                            flat_height=fh, mask=mask)
    U_lid_norm = lid_profile_regularized(geom["xc"], eps) \
        if case["lid"] == "regularized" else np.ones(N)
    diag = compute_diagnostics_fields(
        Ux, Uy, _grid_geom_from_geometry(geom), U_lid_norm,
        case_id=case["id"], Re=Re, N=N, h_over_H=h, Nw=nw,
        lid=case["lid"], eps_lid=eps, U_ref=U0, tau=None, Ma=None)
    diag["iters"] = sol["iters"]
    diag["rel_picard"] = sol["rel"]
    diag["elapsed_s"] = sol["elapsed"]
    diag["nx"] = nx
    diag["ny"] = ny
    return diag, sol, (Ux, Uy, geom)


def _grid_geom_from_geometry(g):
    return {
        "N": g["N"], "dx": g["dx"], "fluid": g["fluid"], "solid": g["solid"],
        "xc": g["xc"], "yc": g["yc"], "X": g["X"], "Y": g["Y"], "yb": g["yb"],
    }


def main(argv=None):
    import cases as C
    ap = argparse.ArgumentParser(description="Taylor-Hood FEM P2/P1 solver")
    ap.add_argument("spec", help="FEM case identifier, e.g. FEM-G06-R400")
    ap.add_argument("--nx", type=int, default=48)
    ap.add_argument("--ny", type=int, default=48)
    ap.add_argument("--N", type=int, default=128)
    ap.add_argument("--max_it", type=int, default=60)
    ap.add_argument("--rtol", type=float, default=1e-10)
    ap.add_argument("--outdir", default="results")
    args = ap.parse_args(argv)

    target = args.spec.upper()
    if "." in target:
        target, sub = target.split(".", 1)
    cases = C.build(["fem"])
    case = None
    for c in cases:
        if c["id"].upper() == target:
            case = c
            break
    if case is None:
        hid = [k for k, v in C.GID.items() if v == target.replace("FEM-", "")]
        if not hid:
            raise SystemExit(f"unknown FEM spec: {args.spec}")
        case = {"id": target, "group": "fem", "h_over_H": hid[0][0],
                "Nw": hid[0][1], "Re": int(args.spec.rsplit("-R", 1)[1]),
                "lid": "regularized", "eps_lid": 0.02, "flat_height": None}

    os.makedirs(args.outdir, exist_ok=True)
    log = open(os.path.join(args.outdir, "fem.log"), "a", encoding="utf-8")
    diag, sol, _ = solve_to_diagnostics(
        case, args.nx, args.ny, N=args.N, max_it=args.max_it, rtol=args.rtol,
        verbose=True, log=log, outdir=args.outdir)
    log.write(f"DONE {case['id']} K*={diag['K_star']:.6f} "
              f"P*={diag['P_star']:.4f} D*={diag['D_star']:.4f} "
              f"vort={diag['n_vortices']} it={sol['iters']} "
              f"({sol['elapsed']:.1f}s)\n")
    log.close()
    od = os.path.join(args.outdir, "fem")
    os.makedirs(od, exist_ok=True)
    meta = dict(diag)
    meta.pop("vortices", None)
    meta["vortices"] = diag.get("vortices")
    with open(os.path.join(od, f"{case['id']}.json"), "w", encoding="utf-8") as fh:
        json.dump(meta, fh, indent=1, default=str)
    print(f"FEM solver finished: {case['id']} it={sol['iters']} "
          f"K*={diag['K_star']:.6f}")
    return diag


if __name__ == "__main__":
    main()