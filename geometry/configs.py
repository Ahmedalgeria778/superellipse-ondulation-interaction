# -*- coding: utf-8 -*-
"""
geometry/configs.py
Configurations pilotes validées (Option A).
"""

from __future__ import annotations
from dataclasses import dataclass


@dataclass(frozen=True)
class CaseConfig:
    name: str
    case_id: str
    n: float
    alpha: float
    Nw: int
    Re: float
    N: int
    a: float = 0.5
    b: float = 0.5
    H: float = 1.0
    U_ref: float = 1.0


PILOT_4CASES = [
    CaseConfig("A", "se_n4_w0p00_Nw1_Re400_N256", n=4.0, alpha=0.00, Nw=1, Re=400, N=256),
    CaseConfig("B", "se_n4_w0p10_Nw1_Re400_N256", n=4.0, alpha=0.10, Nw=1, Re=400, N=256),
    CaseConfig("C", "se_n8_w0p00_Nw1_Re400_N256", n=8.0, alpha=0.00, Nw=1, Re=400, N=256),
    CaseConfig("D", "se_n8_w0p10_Nw1_Re400_N256", n=8.0, alpha=0.10, Nw=1, Re=400, N=256),
]


def nu_lbm_like(U_ref=1.0, H=1.0, Re=400.0):
    return U_ref * H / Re
