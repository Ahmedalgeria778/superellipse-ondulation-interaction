# -*- coding: utf-8 -*-
"""
geometry/superellipse.py
========================
Fond de base superelliptique en repere "boite" (canonique).

Boite : x in [-a, a], y in [0, H].  Le lid est en y = H.
Fond de base : epaulement y = b aux parois (x = ±a), minimum y = 0 au centre.

    relief0(x; n, a, b) = b [ 1 - (1 - (|x|/a)^n)^(1/n) ]

  - relief0(±a) = b   (epaulement)
  - relief0(0)  = 0   (centre, fond le plus bas)
  - n = 2 : arc parabolique/elliptique ; n -> inf : fond quasi rectangulaire
    avec aretes vives pres des parois.
"""

from __future__ import annotations

import numpy as np


def relief0(x, n: float, a: float = 0.5, b: float = 0.5):
    """Fond de base superelliptique (repere boite)."""
    x = np.asarray(x, dtype=np.float64)
    axn = np.clip((np.abs(x) / a) ** n, 0.0, 1.0)
    return b * (1.0 - (1.0 - axn) ** (1.0 / n))


def drelief0_dx(x, n: float, a: float = 0.5, b: float = 0.5):
    """Derivee analytique de relief0 (0 aux parois/centre pour n>1)."""
    x = np.asarray(x, dtype=np.float64)
    axn = np.clip((np.abs(x) / a) ** n, 0.0, 1.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        term = (1.0 - axn) ** (1.0 / n - 1.0)
        val = b * term * (np.abs(x) / a) ** (n - 1.0) * np.sign(x) / a
    val = np.where(np.isfinite(val), val, 0.0)
    return val


def area_under_relief(n: float, a: float = 0.5, b: float = 0.5, Nx: int = 8192):
    """Integrale de relief0 sur [-a, a] (numerique)."""
    x = np.linspace(-a, a, Nx + 1)
    return float(np.trapezoid(relief0(x, n, a, b), x))


if __name__ == "__main__":
    print(relief0(np.array([-0.5, 0.0, 0.5]), 4.0))
    print(area_under_relief(4.0), area_under_relief(8.0))
