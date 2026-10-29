# -*- coding: utf-8 -*-
"""
solver/lbm_combined.py
======================
Solveur D2Q9-MRT pour la geometrie COMBINEE (fond superelliptique + ondulation).

Reutilise, sans le modifier, le noyau Numba valide du Projet A
(`mrt_step_wavy` dans provenance/Lid_driven_ondulation/lbm_wavy.py) et sa
classe `WavyMRT`, en remplacant uniquement la geometrie (masque + liens
Bouzidi) par celle construite dans `solver/combined_geometry.py`.

Toutes les CL (lid, parois, fond), la normalisation et le forçage sont
identiques pour tous les cas (n, alpha, Nw).
"""

from __future__ import annotations

import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
_PROV_A = os.path.join(_ROOT, "provenance", "Lid_driven_ondulation")
for _p in (_HERE, _ROOT, _PROV_A):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from lbm_wavy import WavyMRT  # noqa: E402  (noyau Projet A)
from solver.combined_geometry import build_combined_geometry  # noqa: E402


class CombinedMRT(WavyMRT):
    """Cavite boite + fond superelliptique ondule, lid plat horizontal."""

    def __init__(self, Re, N, n=4.0, alpha=0.0, Nw=1, a=0.5, b=0.5, H=1.0,
                 lid="regularized", eps_lid=0.02, U_ref=None, target_tau=None,
                 max_iter=200000, tol=1e-9, verbose=True, log=None,
                 case_id=None):
        # Geometrie provisoire (sera remplacee juste apres).
        super().__init__(Re=Re, N=N, h_over_H=0.0, Nw=1, lid=lid,
                         eps_lid=eps_lid, U_ref=U_ref, target_tau=target_tau,
                         max_iter=max_iter, tol=tol, verbose=verbose, log=log)
        self.case_id = case_id
        self.n = float(n)
        self.alpha = float(alpha)
        self.Nw = int(Nw)
        self.a = float(a)
        self.b = float(b)
        self.H = float(H)
        self.h_over_H = float(alpha)  # etiquetage diagnostics (=> alpha)

        self.geom = build_combined_geometry(self.N, self.n, self.alpha,
                                            self.Nw, self.a, self.b, self.H)
        self.fluid = self.geom["fluid"]
        self.solid = self.geom["solid"]
        self.fluid_area = self.geom["fluid_area"]
        self.L = self.geom["L"]

    def diagnostics(self):
        from diagnostics import compute_diagnostics
        return compute_diagnostics(self)


if __name__ == "__main__":
    s = CombinedMRT(Re=400, N=64, n=4.0, alpha=0.10, Nw=1, max_iter=2000,
                    verbose=False)
    print("case:", s.n, s.alpha, s.Nw, "area=%.5f" % s.fluid_area,
          "n_bb=%d" % s.geom["n_bb"], "tau=%.4f" % s.tau,
          "U_ref=%.5f" % s.U_ref)
