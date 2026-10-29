# -*- coding: utf-8 -*-
"""Exact square-cavity Ghia benchmark (Point 3 of the review).

Runs the flat-bottom, uniform-lid lid-driven cavity (classical Ghia case)
and dumps the two centerline profiles for comparison with Ghia et al. (1982):

    u(y) along the vertical centreline x = 0.5
    v(x) along the horizontal centreline y = 0.5

Geometry: assembled directly here as a truly flat bottom (y_b = 0) so that
no superelliptic relief is present.  The lid is 'uniform' (no corner
regularisation) to match the Ghia boundary condition exactly.

Usage:
    python run_square_ghia.py <Re> <N> [lid]
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
for _p in (_HERE, _ROOT, os.path.join(_ROOT, "provenance", "Lid_driven_ondulation")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from solver.lbm_combined import CombinedMRT  # noqa: E402
from solver import combined_geometry as cg   # noqa: E402

OUT = os.path.join(_ROOT, "results", "ghia")
LOG = os.path.join(_ROOT, "results", "logs")


def build_flat_geometry(N, a=0.5, b=0.0, H=1.0):
    """Square box, flat bottom at y=0 (solid row removed by q distance)."""
    L = 2.0 * a
    dx = L / N
    dy = H / N
    xc = (np.arange(N) + 0.5) * dx
    yc = (np.arange(N) + 0.5) * dy
    X, Y = np.meshgrid(xc, yc)
    yb = np.zeros_like(X)                 # truly flat bottom
    solid = (Y <= yb).astype(np.uint8)    # all fluid (Y>0 everywhere)
    fluid = (1 - solid).astype(np.uint8)
    bb_q, bb_i, bb_j, bb_qf, bb_low, bb_bi, bb_bj = [], [], [], [], [], [], []
    for i in range(N):
        for j in range(N):
            if not fluid[i, j]:
                continue
            xf = (j + 0.5) * dx
            yf = (i + 0.5) * dy
            for q in range(1, 9):
                pi_ = i - cg.EY[q]
                pj_ = j - cg.EX[q]
                inside = (0 <= pi_ < N) and (0 <= pj_ < N)
                if inside and fluid[pi_, pj_]:
                    continue
                if (pj_ < 0) or (pj_ >= N):
                    continue
                if pi_ >= N:
                    continue
                ex = float(cg.EX[q]); ey = float(cg.EY[q])
                g0 = yf - yb[i, j]
                if g0 <= 0.0:
                    continue
                def _g(s):
                    return (yf - s * ey * dy) - 0.0
                g1 = _g(1.0)
                if g1 >= 0.0:
                    q_frac = 1.0
                else:
                    q_frac = 0.5
                low = 1 if q_frac < 0.5 else 0
                bi, bj = -1, -1
                if low:
                    bi = i + int(cg.EY[q]); bj = j + int(cg.EX[q])
                    if not (0 <= bi < N and 0 <= bj < N and fluid[bi, bj]):
                        low = 0; bi, bj = -1, -1; q_frac = 0.5
                bb_q.append(q); bb_i.append(i); bb_j.append(j)
                bb_qf.append(q_frac); bb_low.append(low)
                bb_bi.append(bi); bb_bj.append(bj)
    geom = {
        "N": N, "n": 0.0, "alpha": 0.0, "Nw": 1, "a": a, "b": b, "H": H, "L": L,
        "dx": dx, "dy": dy, "xc": xc, "yc": yc, "X": X, "Y": Y, "yb": yb,
        "wall_kind": "flat_square",
        "solid": np.ascontiguousarray(solid),
        "fluid": np.ascontiguousarray(fluid),
        "bb_q": np.ascontiguousarray(bb_q, dtype=np.int64),
        "bb_i": np.ascontiguousarray(bb_i, dtype=np.int64),
        "bb_j": np.ascontiguousarray(bb_j, dtype=np.int64),
        "bb_qf": np.ascontiguousarray(bb_qf, dtype=np.float64),
        "bb_low": np.ascontiguousarray(bb_low, dtype=np.int64),
        "bb_bi": np.ascontiguousarray(bb_bi, dtype=np.int64),
        "bb_bj": np.ascontiguousarray(bb_bj, dtype=np.int64),
    }
    geom["n_bb"] = len(bb_q)
    geom["fluid_area"] = float(fluid.sum()) * dx * dy
    return geom


class SquareMRT(CombinedMRT):
    def __init__(self, Re, N, lid="uniform", **kw):
        super().__init__(Re=Re, N=N, n=4.0, alpha=0.0, Nw=1, lid=lid, **kw)
        self.geom = build_flat_geometry(N)
        self.fluid = self.geom["fluid"]
        self.solid = self.geom["solid"]
        self.fluid_area = self.geom["fluid_area"]
        self.L = self.geom["L"]


def main():
    Re = int(sys.argv[1]); N = int(sys.argv[2])
    lid = sys.argv[3] if len(sys.argv) > 3 else "uniform"
    os.makedirs(OUT, exist_ok=True); os.makedirs(LOG, exist_ok=True)
    tag = "Re=%d N=%d lid=%s" % (Re, N, lid)
    t0 = time.time()
    log_path = os.path.join(LOG, "square_Re%d_N%d_%s.log" % (Re, N, lid))
    fh = open(log_path, "w", encoding="utf-8")
    try:
        s = SquareMRT(Re=Re, N=N, lid=lid, max_iter=400000, tol=1e-8,
                      verbose=True, log=fh, case_id="square_Re%d_N%d" % (Re, N))
        s.run()
        # extract velocity
        import lbm_wavy as _lw
        f = s.f if hasattr(s, "f") else None
        ux = getattr(s, "ux", None); uy = getattr(s, "uy", None)
        d = {"Re": Re, "N": N, "lid": lid, "tau": float(s.tau),
             "U_ref": float(s.U_ref), "Ma": float(s.Ma),
             "fluid_area": float(s.fluid_area), "n_bb": int(s.geom["n_bb"])}
        if ux is not None and uy is not None:
            jc = N // 2
            ic = N // 2
            yc = s.geom["yc"]
            xc = s.geom["xc"]
            d["u_centerline"] = {"y": yc.tolist(), "u": np.asarray(ux[:, jc]).tolist()}
            d["v_centerline"] = {"x": xc.tolist(), "v": np.asarray(uy[ic, :]).tolist()}
            d["u_max"] = float(np.max(np.abs([ux[i, jc] for i in range(N)])))
            d["v_min"] = float(np.min([uy[ic, j] for j in range(N)]))
        with open(os.path.join(OUT, "square_Re%d_N%d_%s.json" % (Re, N, lid)), "w") as g:
            json.dump(d, g, indent=2)
        line = "[DONE] %s  %.0fs  u_max=%.4f v_min=%.4f" % (
            tag, time.time() - t0, d.get("u_max", float("nan")),
            d.get("v_min", float("nan")))
    except Exception:
        import traceback
        line = "[FAIL] %s\n%s" % (tag, traceback.format_exc())
    finally:
        fh.close()
    with open(os.path.join(LOG, "square_done.log"), "a", encoding="utf-8") as g:
        g.write(line + "\n")
    print(line)


if __name__ == "__main__":
    main()
