# -*- coding: utf-8 -*-
"""
solver/checks.py
Controles qualite d'un cas convergé : bilan mecanique, Mach, masse,
localisation du tourbillon principal.

Note importante : la cavité entraînée n'est PAS symetrique gauche-droite.
La reflexion x -> -x transforme la vitesse du lid (u_lid, 0) en
(-u_lid, 0) : les conditions aux limites ne sont donc pas invariantes et
le tourbillon principal est décalé du côté aval (sens du lid). Le benchmark
de Ghia et al. (1982) confirme ce décalage (Re=400 -> centre ~ (0.5547,
0.6055)). On ne teste donc pas de symetrie : on suit la position du
tourbillon principal comme indicateur quantitatif.
"""

from __future__ import annotations


def vortex_center(diag):
    """Centre (x, y) du tourbillon principal, ou (nan, nan) si absent."""
    mv = diag.get("main_vortex")
    if not mv:
        return float("nan"), float("nan")
    return float(mv["x"]), float(mv["y"])


def balance(diag):
    """Residu du bilan mecanique r = (P - D)/D."""
    P = diag["P_star"]
    D = diag["D_star"]
    return float((P - D) / D) if D else float("nan")


def summarize(solver, diag):
    vx, vy = vortex_center(diag)
    return {
        "case_id": getattr(solver, "case_id", None),
        "n": solver.n, "alpha": solver.alpha, "Nw": solver.Nw,
        "Re": solver.Re, "N": solver.N,
        "status": solver.status, "iters": solver.iters,
        "residual": float(solver.last_residual),
        "Ma": float(solver.Ma),
        "fluid_area": float(solver.fluid_area),
        "D_star": diag["D_star"], "P_star": diag["P_star"], "K_star": diag["K_star"],
        "Z": diag.get("Z"),
        "balance": balance(diag),
        "vortex_x": vx, "vortex_y": vy,
        "n_vortices": diag.get("n_vortices"),
    }
