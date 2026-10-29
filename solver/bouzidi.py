# -*- coding: utf-8 -*-
"""
solver/bouzidi.py
Interpolated bounce-back (Bouzidi et al. 2001) - skeleton (à affiner).
"""
from __future__ import annotations
import numpy as np

cx = np.array([0,1,0,-1,0,1,-1,-1,1], dtype=np.int32)
cy = np.array([0,0,1,0,-1,1,1,-1,-1], dtype=np.int32)
qopp = np.array([0,3,4,1,2,7,8,5,6])
q = 9
