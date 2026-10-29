# -*- coding: utf-8 -*-
"""
solver/mrt_collision.py
D2Q9 MRT collision (moments space), diagonal relaxation matrix S.
Identique aux versions validées A/B.
"""

from __future__ import annotations
import numpy as np

w = np.array([4/9, 1/9,1/9,1/9,1/9, 1/36,1/36,1/36,1/36], dtype=np.float64)
cx = np.array([0,1,0,-1,0, 1,-1,-1,1], dtype=np.int32)
cy = np.array([0,0,1,0,-1, 1,1,-1,-1], dtype=np.int32)
q = 9

# MRT relaxation matrix (S)
S_mrt = np.array([
    0.0,        # m0
    1.4,        # m1 (e_xx)
    1.4,        # m2 (e_yy)
    0.0,        # m3 (j_x)
    1.2,        # m4 (q_x)
    0.0,        # m5 (j_y)
    1.2,        # m6 (q_y)
    1.0/1.3333, # m7 (p_xx) approx? but original uses 1/tau; better use 1/tau coeff form? wait
], dtype=np.float64)

# Better keep standard form: s7=s8=1/tau
# Store tau-based form set at runtime; provide S_diag(tau)
def S_diag(tau: float) -> np.ndarray:
    s = 1.0 / tau
    return np.array([0.0, 1.4, 1.4, 0.0, 1.2, 0.0, 1.2, s, s], dtype=np.float64)

# MRT transformation (orthogonal D2Q9)
M = np.array([
    [1,1,1,1,1, 1,1,1,1],
    [-4,-1,-1,-1,-1, 2,2,2,2],
    [4,-2,-2,-2,-2, 1,1,1,1],
    [0,1,0,-1,0, 1,-1,-1,1],
    [0,-2,0,2,0, 1,-1,-1,1],
    [0,0,1,0,-1, 1,1,-1,-1],
    [0,0,-2,0,2, 1,1,-1,-1],
    [0,1,-1,1,-1, 0,0,0,0],
    [0,0,0,0,0, 1,-1,1,-1],
], dtype=np.float64)
M_inv = np.linalg.inv(M)
