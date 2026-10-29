# -*- coding: utf-8 -*-
"""
02_LBM/bouzidi_linear.py

Operateur de frontiere de Bouzidi (variante LINEAIRE, 2001), en formulation
"node implementation" (collision + advection combinees, populations
post-collision f*).

Convention (fixee une fois pour toutes) :

    x_f   : noeud fluide
    e_o   : direction incidente (du noeud fluide VERS le mur), o = OPP[q]
    e_q   : direction reflechie (du mur vers le fluide), q = OPP[o]
    q     : frac[q,i,j] = |x_f - x_w| / (|e_o| dx)  in (0,1]

Population reflechie arrivant en x_f dans la direction q :

    q < 1/2 :
        f_q(x_f) = 2q f_o*(x_f) + (1-2q) f_o*(x_f - e_o) + C

    q >= 1/2 :
        f_q(x_f) = (1/(2q)) f_o*(x_f) + (1 - 1/(2q)) f_q*(x_f) + C

avec f_o* = post-collision dans la direction incidente,
     f_q* = post-collision dans la direction reflechie (meme noeud),
     x_f - e_o = voisin derriere x_f, cote fluide.

Propriete fondamentale : a q = 1/2, les deux branches redonnent exactement le
half-way bounce-back f_o*(x_f).
"""

from __future__ import annotations

import numpy as np


def coefficients(q: float):
    """Poids (a_inc, a_refl, a_behind) de l'operateur lineaire de Bouzidi."""
    if q >= 0.5:
        return (1.0 / (2.0 * q), 1.0 - 1.0 / (2.0 * q), 0.0)
    return (2.0 * q, 0.0, 1.0 - 2.0 * q)


def bouzidi_linear(f_incident, f_reflected_post, f_behind, q, C=0.0):
    """Operateur Bouzidi lineaire. Entrees scalaires ou array-like."""
    if q >= 0.5:
        return (1.0 / (2.0 * q)) * f_incident + (1.0 - 1.0 / (2.0 * q)) * f_reflected_post + C
    return 2.0 * q * f_incident + (1.0 - 2.0 * q) * f_behind + C


def halfway_bounce_back(f_incident, C=0.0):
    """Half-way bounce-back (limite q = 1/2)."""
    return f_incident + C


if __name__ == "__main__":
    for q in (0.1, 0.25, 0.5, 0.75, 0.99):
        a, b, c = coefficients(q)
        print(f"q={q:<5} a_inc={a:+.4f} a_refl={b:+.4f} a_behind={c:+.4f}  sum={a + b + c:+.4f}")
    q = 0.5
    lhs = bouzidi_linear(2.0, 7.0, 9.0, q)
    print(f"q=0.5 -> {lhs:.6f} (attendu 2.0 = HBB)")
