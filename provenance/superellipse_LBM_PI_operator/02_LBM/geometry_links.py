# -*- coding: utf-8 -*-
"""
02_LBM/geometry_links.py

Precalcul unique (geometrie fixe) des liens de frontiere pour le solveur courbe.

Pour chaque noeud fluide x_f et chaque direction q dont le voisin x_f + e_q est
solide, on stocke :
    qfrac[q,i,j] = q = |x_f - x_w| / (|e_q| dx)  in (0,1]  (fraction le long du lien)
    uxw[q,i,j], uyw[q,i,j] = vitesse du mur au point d'intersection x_w

C'est le seul endroit ou l'intersection murale est calculee ; le noyau LBM ne
fait que lire ces tableaux.

q est normalise par la longueur du lien |e_q| (1 cardinal, sqrt(2) diagonal) :
c'est le q de Bouzidi.

Reutilise la construction de liens validee de lbm_geometry_grid.SuperellipseGrid.
"""

from __future__ import annotations

import os
import sys
from dataclasses import dataclass

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
for _p in (_HERE, os.path.join(os.path.dirname(_HERE), "01_geometry")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from lbm_geometry_grid import SuperellipseGrid, U_REF  # noqa: E402
from mrt_collision import EX, EY, OPP  # noqa: E402


@dataclass
class LinkArrays:
    fluid: np.ndarray      # (N,N) uint8, 1 = fluide
    qfrac: np.ndarray      # (9,N,N) float64
    uxw: np.ndarray        # (9,N,N) float64
    uyw: np.ndarray        # (9,N,N) float64
    grid: SuperellipseGrid
    N: int

    @property
    def n_links(self) -> int:
        return int((self.qfrac > 0.0).sum())

    def summary(self) -> str:
        return (
            f"links={self.n_links}  fluid={int(self.fluid.sum())}  "
            f"N={self.N}  " + self.grid.summary()
        )


def build_link_arrays(geom, N: int, U_lid: float = U_REF,
                      lid_mode: str = "upper", lid_param=None) -> LinkArrays:
    grid = SuperellipseGrid(geom, N, U_lid=U_lid, lid_mode=lid_mode, lid_param=lid_param)
    L = grid.links
    fluid = grid.fluid.astype(np.uint8)
    qfrac = np.zeros((9, N, N), dtype=np.float64)
    uxw = np.zeros((9, N, N), dtype=np.float64)
    uyw = np.zeros((9, N, N), dtype=np.float64)
    for k in range(L["q"].size):
        q = int(L["q"][k]); i = int(L["i"][k]); j = int(L["j"][k])
        qfrac[q, i, j] = L["frac"][k]
        uxw[q, i, j] = L["ux"][k]
        uyw[q, i, j] = L["uy"][k]
    return LinkArrays(fluid=fluid, qfrac=qfrac, uxw=uxw, uyw=uyw, grid=grid, N=N)


if __name__ == "__main__":
    import argparse

    from superellipse import Superellipse

    parser = argparse.ArgumentParser(description="Precalcul des liens de frontiere")
    parser.add_argument("--n", type=float, default=2.0)
    parser.add_argument("--N", type=int, default=128)
    parser.add_argument("--U", type=float, default=U_REF)
    args = parser.parse_args()

    la = build_link_arrays(Superellipse(n=args.n), args.N, U_lid=args.U)
    q = la.qfrac[la.qfrac > 0.0]
    print(la.summary())
    print(f"q in [{q.min():.4f}, {q.max():.4f}]  mean={q.mean():.4f}")
