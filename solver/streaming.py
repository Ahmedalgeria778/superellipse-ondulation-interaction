# -*- coding: utf-8 -*-
"""
solver/streaming.py
Pull streaming (D2Q9).
"""
from __future__ import annotations
import numpy as np

cx = np.array([0,1,0,-1,0, 1,-1,-1,1], dtype=np.int32)
cy = np.array([0,0,1,0,-1, 1,1,-1,-1], dtype=np.int32)
q = 9

def pull(f, f_new):
    ny, nx = f.shape[1], f.shape[2]
    for i in range(q):
        f_new[i, 1:ny-1, 1:nx-1] = f[i, 1-cy[i]:ny-1-cy[i], 1-cx[i]:nx-1-cx[i]]
    return f_new
