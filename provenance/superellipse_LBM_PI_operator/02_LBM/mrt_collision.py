# -*- coding: utf-8 -*-
"""
02_LBM/mrt_collision.py

Constantes D2Q9 et collision MRT (Diagonale de S), reutilisees telles quelles
depuis le solveur de reference valide (lbm_mrt_reference_square.py).

    U_ref = 0.05, nu_lbm = U_ref * N / Re, tau = 3 nu_lbm + 0.5
    S_diag = [0, 1.4, 1.4, 0, 1.2, 0, 1.2, 1/tau, 1/tau]

Fournit :
    - constantes EX, EY, W, OPP, M, M_INV
    - s_diag(Re, N)
    - macroscopic(f)
    - collide(f, S_diag)  (version numpy, pour tests/unites)

Le noyau Numba utilise par le solveur courbe reprend exactement ces formules.
"""

from __future__ import annotations

import numpy as np

EX = np.array([0, 1, 0, -1, 0, 1, -1, -1, 1], dtype=np.float64)
EY = np.array([0, 0, 1, 0, -1, 1, 1, -1, -1], dtype=np.float64)
W = np.array([4 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 36, 1 / 36, 1 / 36, 1 / 36], dtype=np.float64)
OPP = np.array([0, 3, 4, 1, 2, 7, 8, 5, 6], dtype=np.int32)

M = np.array([
    [1, 1, 1, 1, 1, 1, 1, 1, 1],
    [-4, -1, -1, -1, -1, 2, 2, 2, 2],
    [4, -2, -2, -2, -2, 1, 1, 1, 1],
    [0, 1, 0, -1, 0, 1, -1, -1, 1],
    [0, -2, 0, 2, 0, 1, -1, -1, 1],
    [0, 0, 1, 0, -1, 1, 1, -1, -1],
    [0, 0, -2, 0, 2, 1, 1, -1, -1],
    [0, 1, -1, 1, -1, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 1, -1, 1, -1],
], dtype=np.float64)
M_INV = np.linalg.inv(M)

U_REF = 0.05


def s_diag(Re: float, N: int, L_char: float | None = None,
           U_ref: float = U_REF) -> np.ndarray:
    """Relaxation MRT. Re = U_ref * L_char / nu_lbm.

    L_char=None (defaut, retrocompatible avec le solveur de reference) : la
    longueur caracteristique vaut N (nombre de noeuds). Pour une comparaison
    physique avec la FEM (L_c = 2a), passer L_char = N - 1 (la boite mesure
    N-1 pas de grille).

    U_ref : vitesse caracteristique effective. Doit valoir la vitesse physique
    du couvercle pour que le Re demande soit exactement respecte.
    """
    L = float(N) if L_char is None else float(L_char)
    nu_lbm = U_ref * L / Re
    tau = 3.0 * nu_lbm + 0.5
    omega = 1.0 / tau
    return np.array([0.0, 1.4, 1.4, 0.0, 1.2, 0.0, 1.2, omega, omega], dtype=np.float64)


def macroscopic(f: np.ndarray):
    rho = f.sum(axis=0)
    rho_safe = np.maximum(rho, 1e-10)
    ux = np.tensordot(EX, f, axes=(0, 0)) / rho_safe
    uy = np.tensordot(EY, f, axes=(0, 0)) / rho_safe
    return rho, ux, uy


def collide(f: np.ndarray, S: np.ndarray, fluid: np.ndarray | None = None) -> np.ndarray:
    """Collision MRT node-par-node (numpy). fluid None => tous les noeuds."""
    rho, ux, uy = macroscopic(f)
    rho = np.maximum(rho, 1e-10)
    u2 = ux * ux + uy * uy
    meq = np.empty_like(f)
    meq[0] = rho
    meq[1] = -2.0 * rho + 3.0 * rho * u2
    meq[2] = rho - 3.0 * rho * u2
    meq[3] = rho * ux
    meq[4] = -rho * ux
    meq[5] = rho * uy
    meq[6] = -rho * uy
    meq[7] = rho * (ux * ux - uy * uy)
    meq[8] = rho * ux * uy
    m = np.tensordot(M, f, axes=(1, 0))
    m_star = m - S[:, None, None] * (m - meq)
    f_post = np.tensordot(M_INV, m_star, axes=(1, 0))
    if fluid is not None:
        f_post = np.where(fluid[None, :, :], f_post, 0.0)
    return f_post
