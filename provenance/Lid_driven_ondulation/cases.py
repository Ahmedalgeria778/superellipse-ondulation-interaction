# -*- coding: utf-8 -*-
"""
cases.py
========
Case matrix of the plan :

  - main      : 48 cases (G01..G16) x Re in {100,400,1000}, regularized lid ;
  - flatref   : F00-R{Re} flat wall (h=0) regularized lid, dP* reference ;
  - eq        : EQ05/EQ10 flat wall of equal area (volume check, plan 7.2) ;
  - l1        : G04,G06,G13 x Re with uniform lid (L1 check, plan 7.3) ;
  - validation: V01..V03 flat wall uniform lid (solver validation) ;
  - fem       : 15 cases for the Taylor-Hood finite element solver (plan 8).

A case is a serializable dict ; N (resolution) is set at run time.
"""

from __future__ import annotations

# (h/H, Nw) -> G01..G16 (order: increasing h, then increasing Nw)
GEOMETRIES = [
    (0.025, 1), (0.025, 2), (0.025, 3), (0.025, 4),
    (0.050, 1), (0.050, 2), (0.050, 3), (0.050, 4),
    (0.075, 1), (0.075, 2), (0.075, 3), (0.075, 4),
    (0.100, 1), (0.100, 2), (0.100, 3), (0.100, 4),
]  # G01 .. G16

RES = [100, 400, 1000]

GID = {g: f"G{idx + 1:02d}" for idx, g in enumerate(GEOMETRIES)}


def _case(cid, group, h, nw, Re, lid="regularized", eps=0.02, flat_height=None):
    return {
        "id": cid, "group": group, "h_over_H": h, "Nw": nw, "Re": Re,
        "lid": lid, "eps_lid": eps, "flat_height": flat_height,
    }


def main_cases(lid="regularized"):
    out = []
    for (h, nw) in GEOMETRIES:
        gid = GID[(h, nw)]
        for Re in RES:
            out.append(_case(f"{gid}-R{Re}", "main", h, nw, Re, lid))
    return out


def flatref_cases(lid="regularized"):
    return [_case(f"F00-R{Re}", "flatref", 0.0, 1, Re, lid) for Re in RES]


def eq_cases():
    """Flat wall with the same area as G06 (h/H=0.05 -> y_b=h/2=0.025) and G13/G16
    (h/H=0.10 -> y_b=0.05). They are compared to the wavy average."""
    out = []
    for Re in RES:
        out.append(_case(f"EQ05-R{Re}", "eq", 0.05, 2, Re, "regularized",
                         flat_height=0.025))
        out.append(_case(f"EQ10-R{Re}", "eq", 0.10, 4, Re, "regularized",
                         flat_height=0.05))
    return out


def l1_cases():
    """Sensitivity to lid regularization (uniform lid)."""
    out = []
    for (h, nw) in [(0.025, 4), (0.05, 2), (0.10, 1)]:   # G04, G06, G13
        gid = GID[(h, nw)]
        for Re in RES:
            out.append(_case(f"L1-{gid}-R{Re}", "l1", h, nw, Re, "uniform"))
    return out


def validation_cases():
    """Validation flat wall, uniform lid (to be compared with the reference)."""
    return [_case(f"V{i + 1:02d}-R{Re}", "validation", 0.0, 1, Re, "uniform")
            for i, Re in enumerate(RES)]


def fem_cases():
    """15 FEM cases: reduced h/H x Nw grid, Re in {100,400}."""
    out = []
    for (h, nw) in [(0.025, 2), (0.05, 1), (0.05, 2), (0.075, 2), (0.10, 2),
                    (0.05, 3), (0.10, 3)]:
        for Re in (100, 400):
            out.append(_case(f"FEM-{GID[(h, nw)]}-R{Re}", "fem", h, nw, Re,
                             "regularized"))
    return out[:15]


GROUPS = {
    "main": main_cases,
    "flatref": flatref_cases,
    "eq": eq_cases,
    "l1": l1_cases,
    "validation": validation_cases,
    "fem": fem_cases,
}


def build(groups=("main",), lid="regularized"):
    out = []
    for g in groups:
        fn = GROUPS[g]
        try:
            out += fn(lid=lid) if g in ("main", "flatref") else fn()
        except TypeError:
            out += fn()
    return out


if __name__ == "__main__":
    for g in GROUPS:
        c = GROUPS[g]()
        print(f"[{g:10s}] {len(c):3d} cases  ex: {c[0]['id']} ... {c[-1]['id']}")
