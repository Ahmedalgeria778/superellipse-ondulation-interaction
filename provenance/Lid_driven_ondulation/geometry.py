# -*- coding: utf-8 -*-
"""
geometry.py
===========
Geometry of the lid-driven cavity with a wavy lower wall.

Lower wall (plan, section 1.2):

        y_b(x) = (h/2) [ 1 - cos(2 pi x / lambda) ]
               = (h/2) [ 1 - cos(2 pi N_w x / L) ]

with L = H = 1, N_w = L/lambda in {1,2,3,4}.

Domain conventions (consistent with the reference MRT-LBM solver) :
  - N nodes per direction, indices i (vertical, 0 bottom -> N-1 top), j (horizontal).
  - The fluid domain occupies the box [0,1]x[0,1] in dimensionless coordinates.
  - Physical coordinate of the node (i,j) center: x=(j+0.5)/N, y=(i+0.5)/N.
  - The lid (y=1) and the side walls (x=0,1) are treated with
    half-way bounce-back, as in the reference solver.
  - The wavy lower wall is treated with interpolated bounce-back (Bouzidi),
    the fraction q being computed by exact intersection of the link with y_b(x).

Derived dimensionless quantities (plan, section 1.3) :
  alpha = h/H,  Lambda = lambda/L = 1/N_w,
  S_max = pi alpha N_w          (maximum slope),
  K_max = 2 pi^2 alpha N_w^2    (maximum curvature),
  A_f   = L (H - h/2) = 1 - h/2 (mean fluid area).
"""

from __future__ import annotations

import numpy as np

# ----------------------------------------------------------------------
#  D2Q9 constants (identical to the reference solver)
# ----------------------------------------------------------------------
EX = np.array([0, 1, 0, -1, 0, 1, -1, -1, 1], dtype=np.int64)
EY = np.array([0, 0, 1, 0, -1, 1, 1, -1, -1], dtype=np.int64)
OPP = np.array([0, 3, 4, 1, 2, 7, 8, 5, 6], dtype=np.int64)
ENORM = np.sqrt(EX.astype(np.float64) ** 2 + EY.astype(np.float64) ** 2)


def wall_height(x, h_over_H, Nw, L=1.0):
    """Vertical position y_b(x) of the wavy lower wall (plan 1.2)."""
    x = np.asarray(x, dtype=np.float64)
    return 0.5 * h_over_H * (1.0 - np.cos(2.0 * np.pi * Nw * x / L))


def wall_slope(x, h_over_H, Nw, L=1.0):
    """Local slope dy_b/dx."""
    x = np.asarray(x, dtype=np.float64)
    return 0.5 * h_over_H * (2.0 * np.pi * Nw / L) * np.sin(2.0 * np.pi * Nw * x / L)


def derived_params(h_over_H, Nw, L=1.0, H=1.0):
    """Dimensionless geometric quantities of the plan (section 1.3)."""
    alpha = h_over_H
    return {
        "alpha": alpha,
        "Lambda": 1.0 / Nw,
        "S_max": np.pi * alpha * Nw,
        "K_max": 2.0 * np.pi ** 2 * alpha * Nw ** 2,
        "area": L * (H - 0.5 * alpha * H),
    }


def lid_profile_regularized(x, eps=0.02):
    """Quasi-uniform lid profile (corner regularization, plan 3.2)."""
    x = np.asarray(x, dtype=np.float64)
    g = np.ones_like(x)
    lo = x < eps
    hi = x > 1.0 - eps
    tl = x[lo] / eps
    g[lo] = 3.0 * tl ** 2 - 2.0 * tl ** 3
    th = (1.0 - x[hi]) / eps
    g[hi] = 3.0 * th ** 2 - 2.0 * th ** 3
    return g


def build_geometry(N, h_over_H, Nw, L=1.0, H=1.0, flat_height=None):
    """Builds the solid/fluid mask and the wall boundary links.

    By default: wavy wall y_b(x) (parameters h_over_H, Nw).
    If flat_height is provided (EQ volume check, plan 7.2), the lower
    wall is flat at the level y = flat_height.

    Returns a dictionary containing:
      solid, fluid (N,N) uint8 : masks
      xc, yc, X, Y : physical coordinates of the node centers
      yb           : wall height per node
      bb_*         : wall link lists (Bouzidi interpolated bounce-back)
    """
    dx = L / N
    dy = H / N
    xc = (np.arange(N) + 0.5) * dx
    yc = (np.arange(N) + 0.5) * dy
    X, Y = np.meshgrid(xc, yc)                      # [i, j]
    if flat_height is not None:
        yb = np.full_like(X, float(flat_height))
        wall_kind = "flat"

        def _ybfun(x):
            return np.full_like(np.asarray(x, dtype=np.float64), float(flat_height))
    else:
        yb = wall_height(X, h_over_H, Nw, L)
        wall_kind = "wavy"

        def _ybfun(x):
            return wall_height(x, h_over_H, Nw, L)

    solid = (Y <= yb).astype(np.uint8)
    fluid = (1 - solid).astype(np.uint8)
    if not fluid.any():
        raise ValueError("geometry entirely solid (flat_height too large?)")

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
                    continue  # standard streaming
                if (pj_ < 0) or (pj_ >= N):
                    continue  # flat side wall
                if pi_ >= N:
                    continue  # lid
                # --- wavy lower wall: link / y_b(x) intersection ---
                ex = float(EX[q])
                ey = float(EY[q])
                g0 = yf - float(_ybfun(xf))
                if g0 <= 0.0:
                    continue  # fluid node below the wall (inconsistency) -> half-way bounce-back
                # g(s) = y(s) - yb(x(s)),  x(s)=xf - s ex dx, y(s)=yf - s ey dy
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
                        # no node behind: we fall back to pure half-way bounce-back.
                        # We force q_frac=0.5, otherwise the solver would apply the
                        # q>=0.5 branch to a q_frac<0.5 -> extrapolation
                        # (1/(2q), (2q-1)/(2q)) which diverges (EQ case, bottom corners).
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
        "N": N,
        "h_over_H": h_over_H,
        "Nw": Nw,
        "L": L,
        "H": H,
        "dx": dx,
        "dy": dy,
        "xc": xc,
        "yc": yc,
        "X": X,
        "Y": Y,
        "yb": yb,
        "wall_kind": wall_kind,
        "flat_height": flat_height,
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
        "params": derived_params(h_over_H, Nw, L, H),
    }
    geom["fluid_area"] = float(fluid.sum()) * dx * dy
    return geom


def geometry_id(h_over_H, Nw):
    """Deterministic Gnn identifier from (h/H, Nw) according to the plan (section 4)."""
    hs = [0.025, 0.05, 0.075, 0.10]
    ih = int(np.argmin([abs(h_over_H - v) for v in hs]))
    return f"G{ih * 4 + int(Nw):02d}"


if __name__ == "__main__":
    for (h, nw) in [(0.0, 1), (0.05, 2), (0.10, 4)]:
        g = build_geometry(128, h, nw)
        print(f"h/H={h} Nw={nw} : id={geometry_id(h, nw)} n_bb={g['n_bb']} "
              f"area={g['fluid_area']:.4f} S_max={g['params']['S_max']:.3f}")
