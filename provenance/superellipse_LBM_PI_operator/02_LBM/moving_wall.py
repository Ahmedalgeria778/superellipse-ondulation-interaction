# -*- coding: utf-8 -*-
"""
02_LBM/moving_wall.py

Couvercle tangent mobile : vitesse du mur u_w au point d'intersection et
correction de quantite de mouvement pour le rebond (paroi mobile).

Vitesse tangentielle (convention verrouillee) :
    t = (n_y, -n_x) orientee vers +x,  u_w = U_lid * t
    -> |u_w| = U_lid et u_w . n = 0 (impermeabilite).

Correction de mur mobile ajoutee a la population reflechie d'indice q
(la population qui ARRIVE en x_f, de direction e_q) :

    C_q = 6 w_q rho (e_q . u_w) * s(q_link),    q_link = qfrac du lien vers le mur

Le facteur s(q) est celui de l'operateur lineaire de Bouzidi (verifie sur la
solution exacte de Couette plan ; voir 03_LBM_validation/curved_boundary) :

    s(q) = 1          si q <  1/2
    s(q) = 1/(2q)     si q >= 1/2        (continu en q = 1/2 : s = 1)

A q = 1/2 on retrouve l'injection du solveur de reference valide (couvercle
plat, nord) :
    q = 7 (SW) : e_7 . u_w = -u  ->  C = -6 w rho u
    q = 8 (SE) : e_8 . u_w = +u  ->  C = +6 w rho u
    q = 4 (S)  : e_4 . u_w = 0   ->  C = 0
"""

from __future__ import annotations

import numpy as np

from mrt_collision import EX, EY, W


def wall_scale(qf: float) -> float:
    """Facteur de correction de mur mobile de Bouzidi lineaire."""
    return 1.0 if qf < 0.5 else 1.0 / (2.0 * qf)


def wall_momentum(q: int, ux_w: float, uy_w: float, rho: float = 1.0,
                  qf: float = 0.5) -> float:
    """C_q = 6 w_q rho (e_q . u_w) * s(qf)."""
    return 6.0 * W[q] * rho * (EX[q] * ux_w + EY[q] * uy_w) * wall_scale(qf)


def wall_momentum_array(q: np.ndarray, ux_w, uy_w, rho=1.0, qf=0.5) -> np.ndarray:
    q = np.asarray(q)
    return 6.0 * W[q] * rho * (EX[q] * np.asarray(ux_w) + EY[q] * np.asarray(uy_w)) \
        * wall_scale(qf)


def tangential_velocity(geom, wx, wy, U_lid: float):
    """u_w = U_lid * t aux points (wx, wy) ; t=(n_y,-n_x) orientee vers +x.

    Entrees scalaires ou arrays ; sorties de meme forme.
    """
    wx = np.asarray(wx, dtype=np.float64)
    wy = np.asarray(wy, dtype=np.float64)
    nx, ny = geom.normal(wx, wy)
    nx = np.asarray(nx, dtype=np.float64)
    ny = np.asarray(ny, dtype=np.float64)
    tx, ty = ny, -nx
    return U_lid * tx, U_lid * ty, nx, ny


if __name__ == "__main__":
    rho = 1.0
    u = 0.05
    print("Controle du signe (couvercle plat, u_w=(u,0)) :")
    for q in (4, 7, 8):
        print(f"  q={q}: C={wall_momentum(q, u, 0.0, rho):+.6f}  "
              f"(e_q.u_w={EX[q] * u + EY[q] * 0.0:+.3f})")
