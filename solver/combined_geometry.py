# -*- coding: utf-8 -*-
"""
solver/combined_geometry.py
===========================
Masque solide/fluide + liens Bouzidi pour la geometrie COMBINEE
(fond superelliptique + ondulation), au format attendu par le noyau LBM
Projet A.

Repere du solveur (identique au Projet A) :
    x' in [0, L],  L = 2a
    y' in [0, H],  lid en y' = H   (rangee i = N-1)
    fond :  Y_b(x') = y_b(x' - a; n, alpha, Nw, a, b)  in [0, b]
Cells carrees : dx = dy = L/N = H/N  (boite carree : H = L = 2a).
"""

from __future__ import annotations

import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
_PROV_A = os.path.join(_ROOT, "provenance", "Lid_driven_ondulation")
for _p in (_HERE, _PROV_A):
    if _p not in sys.path:
        sys.path.insert(0, _p)
# NB : _ROOT n'est PAS ajoute a sys.path ici pour ne pas masquer le
# `geometry.py` du Projet A (requis par lbm_wavy). La formule du fond est
# dupliquee ci-dessous (source de verite : geometry/combined_geometry.py).

EX = np.array([0, 1, 0, -1, 0, 1, -1, -1, 1], dtype=np.int64)
EY = np.array([0, 0, 1, 0, -1, 1, 1, -1, -1], dtype=np.int64)
OPP = np.array([0, 3, 4, 1, 2, 7, 8, 5, 6], dtype=np.int64)


def _relief0(x, n, a, b):
    axn = np.clip((np.abs(x) / a) ** n, 0.0, 1.0)
    return b * (1.0 - (1.0 - axn) ** (1.0 / n))


def _eta(x, alpha, Nw, a):
    s = (x + a) / (2.0 * a)
    return (alpha / 2.0) * (1.0 - np.cos(2.0 * np.pi * Nw * s))


def y_b_centered(x, n, alpha, Nw, a, b):
    """Fond combine en repere boite centre x in [-a, a]."""
    x = np.asarray(x, dtype=np.float64)
    return _relief0(x, n, a, b) + _eta(x, alpha, Nw, a)


def build_combined_geometry(N, n, alpha, Nw, a=0.5, b=0.5, H=1.0):
    L = 2.0 * a
    if abs(L - H) > 1e-12:
        raise ValueError("boite non carree : L=%g, H=%g (cells non carrees)" % (L, H))
    dx = L / N
    dy = H / N
    xc = (np.arange(N) + 0.5) * dx
    yc = (np.arange(N) + 0.5) * dy
    X, Y = np.meshgrid(xc, yc)                      # [i, j]

    def _ybfun(xp):
        xp = np.asarray(xp, dtype=np.float64)
        return y_b_centered(xp - a, n, alpha, Nw, a, b)

    yb = _ybfun(X)
    if float(np.max(yb)) >= H:
        raise ValueError("fond au niveau du lid (b trop grand)")
    solid = (Y <= yb).astype(np.uint8)
    fluid = (1 - solid).astype(np.uint8)
    if not fluid.any():
        raise ValueError("geometrie entierement solide")

    bb_q, bb_i, bb_j, bb_qf, bb_low, bb_bi, bb_bj = [], [], [], [], [], [], []

    for i in range(N):
        for j in range(N):
            if not fluid[i, j]:
                continue
            xf = (j + 0.5) * dx
            yf = (i + 0.5) * dy
            for q in range(1, 9):
                pi_ = i - EY[q]
                pj_ = j - EX[q]
                inside = (0 <= pi_ < N) and (0 <= pj_ < N)
                if inside and fluid[pi_, pj_]:
                    continue
                if (pj_ < 0) or (pj_ >= N):
                    continue
                if pi_ >= N:
                    continue
                ex = float(EX[q])
                ey = float(EY[q])
                g0 = yf - float(_ybfun(xf))
                if g0 <= 0.0:
                    continue
                def _g(s):
                    return (yf - s * ey * dy) - float(_ybfun(xf - s * ex * dx))
                g1 = _g(1.0)
                if g1 >= 0.0:
                    q_frac = 1.0
                else:
                    lo, hi = 0.0, 1.0
                    for _ in range(60):
                        mid = 0.5 * (lo + hi)
                        if _g(mid) > 0.0:
                            lo = mid
                        else:
                            hi = mid
                    q_frac = 0.5 * (lo + hi)
                    if q_frac < 1e-8:
                        q_frac = 1e-8
                low = 1 if q_frac < 0.5 else 0
                bi, bj = -1, -1
                if low:
                    bi = i + int(EY[q])
                    bj = j + int(EX[q])
                    if not (0 <= bi < N and 0 <= bj < N and fluid[bi, bj]):
                        low = 0
                        bi, bj = -1, -1
                        q_frac = 0.5
                bb_q.append(q)
                bb_i.append(i)
                bb_j.append(j)
                bb_qf.append(q_frac)
                bb_low.append(low)
                bb_bi.append(bi)
                bb_bj.append(bj)

    geom = {
        "N": N, "n": n, "alpha": alpha, "Nw": Nw,
        "a": a, "b": b, "H": H, "L": L,
        "dx": dx, "dy": dy,
        "xc": xc, "yc": yc, "X": X, "Y": Y, "yb": yb,
        "wall_kind": "combined",
        "solid": np.ascontiguousarray(solid),
        "fluid": np.ascontiguousarray(fluid),
        "bb_q": np.ascontiguousarray(bb_q, dtype=np.int64),
        "bb_i": np.ascontiguousarray(bb_i, dtype=np.int64),
        "bb_j": np.ascontiguousarray(bb_j, dtype=np.int64),
        "bb_qf": np.ascontiguousarray(bb_qf, dtype=np.float64),
        "bb_low": np.ascontiguousarray(bb_low, dtype=np.int64),
        "bb_bi": np.ascontiguousarray(bb_bi, dtype=np.int64),
        "bb_bj": np.ascontiguousarray(bb_bj, dtype=np.int64),
        "n_bb": len(bb_q),
    }
    geom["fluid_area"] = float(fluid.sum()) * dx * dy
    return geom


if __name__ == "__main__":
    for n in (4.0, 8.0):
        for al in (0.0, 0.10):
            g = build_combined_geometry(256, n, al, 1)
            print("n=%g a=%.2f n_bb=%4d area=%.6f" % (n, al, g["n_bb"], g["fluid_area"]))
