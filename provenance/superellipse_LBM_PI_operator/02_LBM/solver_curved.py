# -*- coding: utf-8 -*-
"""
02_LBM/solver_curved.py

Solveur D2Q9-MRT pour geometrie superelliptique (frontiere courbe), base sur :
    - collision MRT identique au solveur de reference valide,
    - streaming PULL,
    - operateur de frontiere Bouzidi LINEAIRE + couvercle tangent mobile.

Le noyau Numba est la traduction exacte de streaming.stream (version numpy) et
de mrt_collision (version numpy), ce qui permet de les tester l'un contre
l'autre.

Sortie : ux, uy normalises par U_ref (comme le solveur de reference).
"""

from __future__ import annotations

import os
import sys
import time

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
for _p in (_HERE, os.path.join(os.path.dirname(_HERE), "01_geometry")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from numba import njit, prange  # noqa: E402

from mrt_collision import (EX, EY, W, OPP, M, M_INV, U_REF, s_diag,  # noqa: E402
                           collide, macroscopic)
from geometry_links import build_link_arrays  # noqa: E402
from streaming import stream as _stream_numpy  # noqa: E402


@njit(parallel=True, cache=True)
def _step(f, M, M_INV, S_diag, EX, EY, W, OPP, fluid, qfrac, uxw, uyw):
    N = f.shape[1]
    rho = np.zeros((N, N))
    ux = np.zeros((N, N))
    uy = np.zeros((N, N))
    f_post = np.zeros_like(f)

    for i in prange(N):
        for j in range(N):
            if fluid[i, j] == 0:
                continue
            r = 0.0
            mx = 0.0
            my = 0.0
            for q in range(9):
                r += f[q, i, j]
                mx += EX[q] * f[q, i, j]
                my += EY[q] * f[q, i, j]
            r = max(r, 1e-10)
            ux_ = mx / r
            uy_ = my / r
            rho[i, j] = r
            ux[i, j] = ux_
            uy[i, j] = uy_
            u2 = ux_ * ux_ + uy_ * uy_
            meq = np.zeros(9)
            meq[0] = r
            meq[1] = -2.0 * r + 3.0 * r * u2
            meq[2] = r - 3.0 * r * u2
            meq[3] = r * ux_
            meq[4] = -r * ux_
            meq[5] = r * uy_
            meq[6] = -r * uy_
            meq[7] = r * (ux_ * ux_ - uy_ * uy_)
            meq[8] = r * ux_ * uy_
            m = np.zeros(9)
            for q in range(9):
                s = 0.0
                for k in range(9):
                    s += M[q, k] * f[k, i, j]
                m[q] = s
            ms = np.zeros(9)
            for q in range(9):
                ms[q] = m[q] - S_diag[q] * (m[q] - meq[q])
            for q in range(9):
                s = 0.0
                for k in range(9):
                    s += M_INV[q, k] * ms[k]
                f_post[q, i, j] = s

    f_new = np.zeros_like(f)
    for i in prange(N):
        for j in range(N):
            if fluid[i, j] == 0:
                continue
            for q in range(9):
                si = i - int(EY[q])
                sj = j - int(EX[q])
                if 0 <= si < N and 0 <= sj < N and fluid[si, sj] == 1:
                    f_new[q, i, j] = f_post[q, si, sj]
                    continue
                o = OPP[q]
                # lien vers le mur en direction o = OPP[q] (source en q solide)
                qf = qfrac[o, i, j]
                if qf <= 0.0:
                    f_new[q, i, j] = f_post[o, i, j]
                    continue
                if qf < 0.5:
                    bi = i + int(EY[q])
                    bj = j + int(EX[q])
                    if 0 <= bi < N and 0 <= bj < N and fluid[bi, bj] == 1:
                        fb = f_post[o, bi, bj]
                    else:
                        fb = f_post[o, i, j]
                    base = 2.0 * qf * f_post[o, i, j] + (1.0 - 2.0 * qf) * fb
                    wscale = 1.0
                else:
                    base = (1.0 / (2.0 * qf)) * f_post[o, i, j] \
                        + (1.0 - 1.0 / (2.0 * qf)) * f_post[q, i, j]
                    wscale = 1.0 / (2.0 * qf)
                base += wscale * 6.0 * W[q] * (EX[q] * uxw[o, i, j] + EY[q] * uyw[o, i, j]) * rho[i, j]
                f_new[q, i, j] = base
    return f_new, rho, ux, uy


class CurvedLBM:
    def __init__(self, geom, Re, N, U_lid=U_REF, max_iter=200000, tol=1e-8,
                 lid_mode="upper", lid_param=None, L_char=None):
        self.geom = geom
        self.Re = Re
        self.N = N
        self.U_lid = U_lid
        self.max_iter = max_iter
        self.tol = tol
        self.L_char = L_char
        self.la = build_link_arrays(geom, N, U_lid=U_lid, lid_mode=lid_mode,
                                    lid_param=lid_param)
        self.S = s_diag(Re, N, L_char=L_char, U_ref=U_lid)
        rho0 = np.ones((N, N))
        self.f = np.ascontiguousarray(W[:, None, None] * rho0[None, :, :])

    def step_numpy(self):
        rho, ux, uy = macroscopic(self.f)
        rho = np.maximum(rho, 1e-10)
        f_post = collide(self.f, self.S, self.la.fluid.astype(bool))
        f_new = _stream_numpy(f_post, self.la, rho)
        self.f = f_new
        return rho, ux, uy

    def step_fast(self):
        f_new, rho, ux, uy = _step(
            self.f, M, M_INV, self.S, EX, EY, W, OPP,
            self.la.fluid, self.la.qfrac, self.la.uxw, self.la.uyw)
        self.f = f_new
        return rho, ux, uy

    def run(self, verbose=True, use_fast=True, logfile=None):
        N = self.N
        fluid = self.la.fluid.astype(bool)
        ux_prev = np.zeros((N, N))
        hist = []
        t0 = time.time()
        step = self.step_fast if use_fast else self.step_numpy
        for it in range(1, self.max_iter + 1):
            rho, ux, uy = step()
            if it % 100 == 0 or it == 1:
                du = np.abs(ux - ux_prev)[fluid].max() if fluid.any() else 0.0
                ux_prev = ux.copy()
                hist.append((it, du))
                if verbose and (it % max(1, self.max_iter // 20) == 0 or it == 1):
                    print(f"  it={it:>7d}  dUmax={du:.3e}  t={time.time() - t0:.1f}s")
                if du < self.tol and it > 500:
                    if verbose:
                        print(f"  converge a it={it} (dUmax={du:.3e})")
                    break
        rho, ux, uy = macroscopic(self.f)
        self.iters = it
        self.residual = hist[-1][1] if hist else float("nan")
        self.history = hist
        self.ux = ux / self.U_lid
        self.uy = uy / self.U_lid
        self.rho = rho
        return self

    def save(self, path):
        np.savez_compressed(
            path, ux=self.ux, uy=self.uy, rho=self.rho, fluid=self.la.fluid,
            Re=self.Re, N=self.N, iters=getattr(self, "iters", -1),
            residual=getattr(self, "residual", float("nan")),
            n=(np.inf if not hasattr(self.geom, "n") else self.geom.n),
            U_lid=self.U_lid, Re_u=1.0 / self.Re,
            L_char=(-1.0 if self.L_char is None else self.L_char),
            nu_lbm=self.U_lid * (self.N if self.L_char is None else self.L_char) / self.Re)


def _main() -> int:
    import argparse

    from superellipse import Superellipse

    parser = argparse.ArgumentParser(description="Solveur LBM courbe (Bouzidi lineaire)")
    parser.add_argument("--n", type=float, default=2.0)
    parser.add_argument("--Re", type=float, default=100.0)
    parser.add_argument("--N", type=int, default=128)
    parser.add_argument("--U", type=float, default=U_REF)
    parser.add_argument("--max-iter", type=int, default=200000)
    parser.add_argument("--tol", type=float, default=1e-8)
    parser.add_argument("--numpy", action="store_true", help="force la version numpy (lente)")
    parser.add_argument("--out", type=str, default=None)
    args = parser.parse_args()

    sol = CurvedLBM(Superellipse(n=args.n), args.Re, args.N, U_lid=args.U,
                    max_iter=args.max_iter, tol=args.tol)
    print(f"geom: {sol.geom.summary()}")
    print(f"solver: Re={args.Re} N={args.N} links={sol.la.n_links} U={args.U}")
    sol.run(verbose=True, use_fast=not args.numpy)
    print(f"iters={sol.iters} residual={sol.residual:.3e}")
    if args.out:
        sol.save(args.out)
        print(f"sauve: {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(_main())
