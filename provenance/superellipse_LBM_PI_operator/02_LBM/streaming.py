# -*- coding: utf-8 -*-
"""
02_LBM/streaming.py

Streaming PULL + application de l'operateur de frontiere (Bouzidi lineaire +
correction de mur mobile). Implementation numpy de reference, lente mais
lisible, utilisee pour tester le noyau Numba du solveur.

    f_new[q,i,j] = f_post[q, i-ey_q, j-ex_q]        si la source est fluide
                 = Bouzidi(qfrac[q,i,j], ...) + C   sinon
"""

from __future__ import annotations

import numpy as np

from mrt_collision import EX, EY, OPP
from bouzidi_linear import bouzidi_linear
from moving_wall import wall_momentum


def stream(f_post: np.ndarray, la, rho: np.ndarray) -> np.ndarray:
    fluid = la.fluid
    qfrac = la.qfrac
    uxw = la.uxw
    uyw = la.uyw
    N = la.N
    f_new = np.zeros_like(f_post)
    for i in range(N):
        for j in range(N):
            if fluid[i, j] == 0:
                continue
            for q in range(9):
                si = i - int(EY[q])
                sj = j - int(EX[q])
                if 0 <= si < N and 0 <= sj < N and fluid[si, sj] == 1:
                    f_new[q, i, j] = f_post[q, si, sj]
                    continue
                o = int(OPP[q])
                # le lien vers le mur est en direction o = OPP[q] (source en q solide)
                qf = qfrac[o, i, j]
                if qf <= 0.0:
                    f_new[q, i, j] = f_post[o, i, j]
                    continue
                if qf < 0.5:
                    bi = i + int(EY[q])
                    bj = j + int(EX[q])
                    if 0 <= bi < N and 0 <= bj < N and fluid[bi, bj] == 1:
                        f_behind = f_post[o, bi, bj]
                    else:
                        f_behind = f_post[o, i, j]
                    base = bouzidi_linear(f_post[o, i, j], 0.0, f_behind, qf)
                else:
                    base = bouzidi_linear(f_post[o, i, j], f_post[q, i, j], 0.0, qf)
                base += wall_momentum(q, uxw[o, i, j], uyw[o, i, j], rho[i, j], qf)
                f_new[q, i, j] = base
    return f_new
