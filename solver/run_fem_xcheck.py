# -*- coding: utf-8 -*-
"""
solver/run_fem_xcheck.py
========================
Independent-solver cross-check of the superelliptic campaign.

The Taylor-Hood P2/P1 FEM solver in fem/fem_taylor_hood.py natively uses the
pure-cosine wall of the earlier plan.  This driver monkeypatches its wall
function and its geometry builder so that it solves the SAME superelliptic
geometry as the MRT-LBM campaign (Eq. y_b of the main text), then projects
onto the identical regular grid and evaluates the identical diagnostics
(D_star, P_star, K_star).  The result is compared, case by case, with the
fixed-tau LBM scalar JSON.

Usage (from the project root):
    python solver/run_fem_xcheck.py --Re 400 --n 2 --alpha 0.0 --nx 64 --ny 64 --N 128
"""

from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
_FEM = os.path.join(_ROOT, "fem")
_PROV = os.path.join(_ROOT, "provenance", "Lid_driven_ondulation")
for _p in (_HERE, _FEM, _PROV):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import combined_geometry as CG          # noqa: E402
import fem_taylor_hood as F             # noqa: E402

A = B = 0.5


def _superellipse_yb(x, n, alpha, Nw):
    """Bottom wall in the FEM [0,1] x-convention (same as LBM campaign)."""
    x = np.asarray(x, dtype=np.float64)
    return CG.y_b_centered(x - A, n, alpha, Nw, A, B)


def install_superellipse(n, alpha, Nw):
    """Monkeypatch the FEM wall function and geometry builder in place."""
    def _wall_values(x, h_over_H=None, Nw_=None, flat_height=None):
        return _superellipse_yb(x, n, alpha, Nw)

    def _build_geometry(N, h_over_H=None, Nw_=None, flat_height=None, **kw):
        return CG.build_combined_geometry(N, n, alpha, Nw, a=A, b=B)

    F.wall_values = _wall_values
    F.build_geometry = _build_geometry


def _close(a, b):
    if b == 0:
        return abs(a)
    return abs(a - b) / abs(b)


def fem_mesh_area(mesh):
    """Exact fluid area of the FEM mesh: sum of element triangles.

    Independent of any projection onto the regular grid, so it can be
    checked against the analytic superellipse area (O(h^2) convergence).
    """
    nodes = mesh["nodes"]
    area = 0.0
    for e in mesh["elems"]:
        ids = [int(x) for x in e]
        VA, VB, VC = nodes[ids[0]], nodes[ids[1]], nodes[ids[2]]
        det = (VB[0] - VA[0]) * (VC[1] - VA[1]) \
            - (VB[1] - VA[1]) * (VC[0] - VA[0])
        area += abs(det) / 2.0
    return area


def fem_native_dissipation(mesh, U, M):
    """D* = 2 int S:S dA computed directly from the P2 velocity field.

    No projection onto the LBM grid and no LBM wall-gradient correction are
    used: the strain rate is evaluated from the finite-element velocity
    derivatives inside each element, with the element Jacobian.  The
    integrand S:S is quadratic in (x,y) within a P2 triangle, so the 6-point
    quadrature already used for assembly integrates it to high accuracy.
    """
    nodes = mesh["nodes"]
    elems = mesh["elems"]
    ux = U[:M]
    uy = U[M:]
    D = 0.0
    for e in elems:
        ids = [int(x) for x in e]
        VA, VB, VC = nodes[ids[0]], nodes[ids[1]], nodes[ids[2]]
        JT = np.column_stack((VB - VA, VC - VA))
        det = JT[0, 0] * JT[1, 1] - JT[0, 1] * JT[1, 0]
        JinvT = np.array([[JT[1, 1], -JT[1, 0]],
                          [-JT[0, 1], JT[0, 0]]]) / det
        uxe = ux[ids]
        uye = uy[ids]
        for (xq, yq), wq in zip(F.QUAD_PTS, F.QUAD_WTS):
            _, gr = F._shape_p2(xq, yq)
            px = JinvT[0, 0] * gr[0] + JinvT[0, 1] * gr[1]   # dN/dx
            py = JinvT[1, 0] * gr[0] + JinvT[1, 1] * gr[1]   # dN/dy
            dudx = float(px @ uxe)
            dudy = float(py @ uxe)
            dvdx = float(px @ uye)
            dvdy = float(py @ uye)
            S11, S22, S12 = dudx, dvdy, 0.5 * (dudy + dvdx)
            D += 2.0 * (S11 ** 2 + S22 ** 2 + 2.0 * S12 ** 2) \
                * wq * abs(det)
    return D


def _lid_edges(mesh):
    """Lid-boundary P2 edges as (ia, ib, imid, ids_of_parent_triangle).

    The lid is the physical line y = 1.  A lid edge is a triangle side whose
    two vertices both satisfy y = 1; imid is the P2 midpoint node on that
    side.  ids are the three linear-vertex node indices of the parent
    triangle, ordered to match the reference triangle (0,0),(1,0),(0,1).
    """
    nodes = mesh["nodes"]
    edges = []
    for e in mesh["elems"]:
        ids = [int(x) for x in e]
        p = nodes[ids[:3]]
        on = np.isclose(p[:, 1], 1.0, atol=1e-12)
        if int(on.sum()) != 2:
            continue
        two = np.where(on)[0]
        ia, ib = ids[two[0]], ids[two[1]]
        imid = None
        for m in ids[3:]:
            mp = nodes[m]
            if np.isclose(mp[1], 1.0, atol=1e-12):
                imid = m
                break
        edges.append((ia, ib, imid, ids))
    return edges


def fem_native_power(mesh, U, M, eps=0.02):
    """P* = int_lid U_lid(x) (du/dy)|_lid dx from the P2 field (no projection).

    Evaluated on the lid boundary edges.  At each Gauss point on the edge the
    parent P2 element gradient is obtained from the affine reference map
    (orientation-independent), and the lid velocity is the regularized
    production profile.
    """
    nodes = mesh["nodes"]
    ux = U[:M]
    edges = _lid_edges(mesh)
    gp = np.array([-np.sqrt(3.0 / 5.0), 0.0, np.sqrt(3.0 / 5.0)])
    gw = np.array([5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0])
    P = 0.0
    for (ia, ib, imid, ids) in edges:
        VA = nodes[ia]
        VB = nodes[ib]
        P1, P2, P3 = nodes[ids[0]], nodes[ids[1]], nodes[ids[2]]
        JT = np.column_stack((P2 - P1, P3 - P1))
        det = JT[0, 0] * JT[1, 1] - JT[0, 1] * JT[1, 0]
        Jinv = np.array([[JT[1, 1], -JT[0, 1]],
                         [-JT[1, 0], JT[0, 0]]]) / det
        for s, w in zip(gp, gw):
            # physical point on the edge (linear parametrisation is exact
            # because the edge is a straight lid segment)
            X = 0.5 * (VA + VB) + 0.5 * (VB - VA) * s
            r = Jinv @ (X - P1)          # reference coords in the triangle
            _, gr = F._shape_p2(r[0], r[1])
            # physical gradient = Jinv^T @ reference gradient
            dudy = float((gr[0] * Jinv[0, 1] + gr[1] * Jinv[1, 1]) @ ux[ids])
            ux_lid = float(F.lid_profile_regularized(np.array([X[0]]), eps)[0])
            P += ux_lid * dudy * w * 0.5 * np.hypot(*(VB - VA))
    return P


def fem_native_kinetic(mesh, U, M):
    """K* = 0.5 int |u|^2 dA from the P2 velocity field (element quadrature)."""
    nodes = mesh["nodes"]
    ux = U[:M]
    uy = U[M:]
    K = 0.0
    for e in mesh["elems"]:
        ids = [int(x) for x in e]
        P1, P2, P3 = nodes[ids[0]], nodes[ids[1]], nodes[ids[2]]
        JT = np.column_stack((P2 - P1, P3 - P1))
        det = JT[0, 0] * JT[1, 1] - JT[0, 1] * JT[1, 0]
        uxe = ux[ids]
        uye = uy[ids]
        for (xq, yq), wq in zip(F.QUAD_PTS, F.QUAD_WTS):
            Nv, _ = F._shape_p2(xq, yq)
            uu = float(Nv @ uxe)
            vv = float(Nv @ uye)
            K += 0.5 * (uu ** 2 + vv ** 2) * wq * abs(det)
    return K


def main(argv=None):
    ap = argparse.ArgumentParser(description="FEM superellipse cross-check")
    ap.add_argument("--Re", type=float, default=400.0)
    ap.add_argument("--n", type=float, default=2.0)
    ap.add_argument("--alpha", type=float, default=0.0)
    ap.add_argument("--Nw", type=int, default=1)
    ap.add_argument("--nx", type=int, default=64)
    ap.add_argument("--ny", type=int, default=64)
    ap.add_argument("--N", type=int, default=128)
    ap.add_argument("--max_it", type=int, default=80)
    ap.add_argument("--rtol", type=float, default=1e-11)
    ap.add_argument("--outdir", default=os.path.join(_ROOT, "results", "fem"))
    args = ap.parse_args(argv)

    install_superellipse(args.n, args.alpha, args.Nw)

    case = {
        "id": "FEMX-n%g-a%.2f-R%g-nx%d" % (args.n, args.alpha, args.Re, args.nx),
        "group": "xcheck",
        "h_over_H": None,
        "Nw": args.Nw,
        "Re": int(args.Re),
        "lid": "regularized",
        "eps_lid": 0.02,
        "flat_height": None,
    }

    os.makedirs(args.outdir, exist_ok=True)
    log = open(os.path.join(args.outdir, "fem_xcheck.log"), "a", encoding="utf-8")
    diag, sol, (Ux, Uy, geom) = F.solve_to_diagnostics(
        case, args.nx, args.ny, N=args.N, max_it=args.max_it,
        rtol=args.rtol, verbose=True, log=log, outdir=args.outdir)
    log.write("DONE %s K*=%.6f P*=%.4f D*=%.4f it=%d (%.1fs)\n"
              % (case["id"], diag["K_star"], diag["P_star"], diag["D_star"],
                 sol["iters"], sol["elapsed"]))
    log.close()

    np.savez(os.path.join(args.outdir, case["id"] + "_field.npz"),
             ux=Ux, uy=Uy, fluid=geom["fluid"], xc=geom["xc"], yc=geom["yc"],
             yb=geom["yb"], N=args.N, Re=args.Re, n=args.n, alpha=args.alpha,
             nx=args.nx, ny=args.ny)

    # --- native FEM field (mesh + solution vector), for reproducibility ---
    mesh = sol["mesh"]
    np.savez(os.path.join(args.outdir, case["id"] + "_natfield.npz"),
             nodes=mesh["nodes"], elems=mesh["elems"], U=sol["U"], M=sol["M"],
             nx=args.nx, ny=args.ny, Re=args.Re, n=args.n, alpha=args.alpha,
             eps_lid=case["eps_lid"])

    D_native = fem_native_dissipation(mesh, sol["U"], sol["M"])
    diag["D_star_native"] = D_native
    P_native = fem_native_power(mesh, sol["U"], sol["M"],
                                eps=case["eps_lid"])
    diag["P_star_native"] = P_native
    K_native = fem_native_kinetic(mesh, sol["U"], sol["M"])
    diag["K_star_native"] = K_native
    diag["fluid_area_mesh"] = fem_mesh_area(mesh)
    diag["balance_native"] = (P_native - D_native) / D_native \
        if D_native != 0 else None

    meta = {k: v for k, v in diag.items() if k != "vortices"}
    meta["vortices"] = diag.get("vortices")
    meta["nx"] = args.nx
    meta["ny"] = args.ny
    meta["iters"] = sol["iters"]
    meta["elapsed_s"] = sol["elapsed"]
    with open(os.path.join(args.outdir, case["id"] + ".json"), "w",
              encoding="utf-8") as fh:
        json.dump(meta, fh, indent=1, default=str)

    # --- compare with the LBM campaign scalar, if present ---
    tag = "se_n%g_w%.2f_Nw%d_Re%g_N%d_tol1e-08.json" % (
        args.n, args.alpha, args.Nw, args.Re, args.N)
    lbm_path = os.path.join(_ROOT, "results", "scalars", tag)
    print("FEM  %s: K*=%.6f  P*=%.5f  D*=%.6f (proj)  D*=%.6f (native)"
          "  P*=%.5f (native)  K*=%.6f (native)  (it=%d, %.1fs)"
          % (case["id"], diag["K_star"], diag["P_star"], diag["D_star"],
             D_native, P_native, K_native, sol["iters"], sol["elapsed"]))
    print("  native balance P/D-1 = %+.4f  |  mesh area = %.6f  |  nx=%d"
          % (diag["balance_native"], diag["fluid_area_mesh"], args.nx))
    if os.path.exists(lbm_path):
        lbm = json.load(open(lbm_path, encoding="utf-8"))
        for key in ("K_star", "P_star", "D_star"):
            fv, lv = diag[key], lbm.get(key)
            if lv is not None:
                print("  %-10s FEM=%.6f  LBM=%.6f  rel=%.3f%%"
                      % (key, fv, lv, 100.0 * _close(fv, lv)))
        lv = lbm.get("D_star")
        if lv is not None:
            print("  %-10s FEM=%.6f  LBM=%.6f  rel=%.3f%%"
                  % ("D_star_nat", D_native, lv,
                     100.0 * _close(D_native, lv)))
    else:
        print("  (no LBM scalar found at %s)" % tag)
    return diag


if __name__ == "__main__":
    main()
