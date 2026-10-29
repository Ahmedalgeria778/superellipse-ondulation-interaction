# -*- coding: utf-8 -*-
"""
geometry/combined_geometry.py
=============================
Geometrie combinee canonique (repere boite).

Boite : x in [-a, a], y in [0, H] ; lid en y = H.
Fond  : y_b(x) = relief0(x; n, a, b) + eta(x; alpha, Nw, a)

  - y_b(±a) = b   (epaulement, continuite avec les parois verticales)
  - y_b(0)  = eta(0) = alpha (si Nw=1) : centre le plus bas
"""

from __future__ import annotations

import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

from superellipse import relief0, drelief0_dx          # noqa: E402
from wall_functions import eta, deta_dx                # noqa: E402


def y_b(x, n: float = 4.0, alpha: float = 0.0, Nw: int = 1,
        a: float = 0.5, b: float = 0.5):
    return relief0(x, n, a, b) + eta(x, alpha, Nw, a)


def dy_b_dx(x, n: float = 4.0, alpha: float = 0.0, Nw: int = 1,
            a: float = 0.5, b: float = 0.5):
    return drelief0_dx(x, n, a, b) + deta_dx(x, alpha, Nw, a)


def area_fluid(n: float = 4.0, alpha: float = 0.0, Nw: int = 1,
               a: float = 0.5, b: float = 0.5, H: float = 1.0, Nx: int = 8192):
    x = np.linspace(-a, a, Nx + 1)
    yb = y_b(x, n, alpha, Nw, a, b)
    return float(np.trapezoid(H - yb, x))


def min_yb(n: float = 4.0, alpha: float = 0.0, Nw: int = 1,
           a: float = 0.5, b: float = 0.5, Nx: int = 8192):
    x = np.linspace(-a, a, Nx + 1)
    return float(np.min(y_b(x, n, alpha, Nw, a, b)))


if __name__ == "__main__":
    for n in (4.0, 8.0):
        for al in (0.0, 0.10):
            print("n=%g alpha=%.2f y_b(±a)=b=%.2f y_b_min=%.4f A_f=%.5f"
                  % (n, al, 0.5, min_yb(n, al), area_fluid(n, al)))
