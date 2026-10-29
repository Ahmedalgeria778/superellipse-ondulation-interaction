# -*- coding: utf-8 -*-
"""
geometry/wall_functions.py
==========================
Ondulation eta(x; alpha, Nw, a) s'annulant aux extremites x = ±a.
"""

from __future__ import annotations

import numpy as np


def eta(x, alpha: float = 0.0, Nw: int = 1, a: float = 0.5):
    """Ondulation periodisee sur la largeur 2a, nulle en x = ±a.

    eta(x; alpha,Nw,a) = (alpha/2) [1 - cos(2 pi Nw (x+a)/(2a))]
    """
    x = np.asarray(x, dtype=np.float64)
    s = (x + a) / (2.0 * a)
    return (alpha / 2.0) * (1.0 - np.cos(2.0 * np.pi * Nw * s))


def deta_dx(x, alpha: float = 0.0, Nw: int = 1, a: float = 0.5):
    x = np.asarray(x, dtype=np.float64)
    s = (x + a) / (2.0 * a)
    ds_dx = 1.0 / (2.0 * a)
    return (alpha / 2.0) * (2.0 * np.pi * Nw) * np.sin(2.0 * np.pi * Nw * s) * ds_dx
