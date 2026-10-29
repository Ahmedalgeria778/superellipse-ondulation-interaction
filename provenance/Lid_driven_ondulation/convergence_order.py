# -*- coding: utf-8 -*-
"""Observed convergence order + Richardson extrapolation (non-uniform levels).

Three levels Q(h1,h2,h3) with h = 1/N (h1 > h2 > h3):
    (Q1 - Q2)/(Q2 - Q3) = (h1^p - h2^p)/(h2^p - h3^p)     -> p observed (solved by Brent)
    C = (Q2 - Q3)/(h2^p - h3^p)
    Q_inf = Q3 - C * h3^p                                   -> limit h -> 0

We report the residual gap at N_max as a % of Q_inf: it is the
discretization error remaining at N_max. For the differential effect Delta_Q,
the same treatment gives the uncertainty on the geometric effect.

Output: results/analysis/convergence_order.csv
"""

from __future__ import annotations

import csv
import json
import os
import sys

from scipy.optimize import brentq

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

RES = os.path.join("results", "convergence", "resolution.csv")
TABLE = os.path.join("results", "scalar_table.json")
OUT = os.path.join("results", "analysis", "convergence_order.csv")

# (series, quantity); "Delta_Gxx-Ryyy" = relative effect vs F00-Ryyy at the same N
SERIES = [
    ("F00-R100", "P_star"), ("F00-R100", "D_star"),
    ("F00-R400", "P_star"), ("F00-R400", "D_star"),
    ("F00-R1000", "P_star"), ("F00-R1000", "D_star"),
    ("G06-R400", "P_star"), ("G06-R400", "D_star"),
    ("G12-R400", "P_star"), ("G12-R400", "D_star"),
    ("Delta_G06-R400", "P_star"), ("Delta_G12-R400", "P_star"),
]


def load():
    rows = {}
    with open(TABLE, encoding="utf-8") as fh:
        for r in json.load(fh):
            rows[(r["case"], int(r["N"]))] = {k: float(r[k])
                                              for k in ("P_star", "D_star", "K_star")}
    if os.path.exists(RES):
        with open(RES, newline="", encoding="utf-8") as fh:
            for r in csv.DictReader(fh):
                if r.get("variant", "consistent") != "consistent":
                    continue
                rows[(r["case"], int(r["N"]))] = {k: float(r[k])
                                                  for k in ("P_star", "D_star")}
    return rows


def solve_order(h, Q):
    target = (Q[0] - Q[1]) / (Q[1] - Q[2])

    def f(p):
        return (h[0] ** p - h[1] ** p) / (h[1] ** p - h[2] ** p) - target

    try:
        return float(brentq(f, 0.05, 6.0))
    except ValueError:
        return float("nan")


def main():
    rows = load()
    out = []
    print(f"{'series':17} {'Q':>8} {'N':>16} {'p':>6} {'Q_inf':>10} {'resid at Nmax':>13}")
    for name, key in SERIES:
        if name.startswith("Delta_"):
            case = name.split("_", 1)[1]
            Ns = sorted({n for (c, n) in rows if c == case}
                        & {n for (c, n) in rows if c == "F00-" + case.split("-", 1)[1]})
            vals = [(n, 100.0 * (rows[(case, n)][key]
                                  - rows[("F00-" + case.split("-", 1)[1], n)][key])
                     / rows[("F00-" + case.split("-", 1)[1], n)][key]) for n in Ns]
            label = name.replace("Delta_", "d")
        else:
            vals = sorted((n, rows[(name, n)][key]) for (c, n) in rows if c == name)
            label, key = name, key

        if len(vals) < 3:
            print(f"{label:17} {key:>8} {'/'.join(str(n) for n, _ in vals):>16}"
                  "   (n<3 levels, no extrapolation)")
            continue

        # Richardson analysis triplet: {128,256,512} (refinement ratio 2)
        # when available, otherwise the three largest grids.
        dvals = dict(vals)
        prefer = [128, 256, 512]
        if all(n in dvals for n in prefer):
            Ns3 = prefer
            Q = [dvals[n] for n in prefer]
        else:
            Ns3 = [n for n, _ in vals][:3]
            Q = [v for _, v in vals][:3]
        h = [1.0 / n for n in Ns3]
        p = solve_order(h, Q)
        if not (p == p):
            print(f"{label:17} {key:>8} {'/'.join(map(str, Ns3)):>16}      p unavailable")
            continue
        C = (Q[1] - Q[2]) / (h[1] ** p - h[2] ** p)
        Qinf = Q[2] - C * h[2] ** p
        resid = 100.0 * (Q[2] - Qinf) / Qinf
        out.append({"series": label, "quantity": key,
                    "N_used": "/".join(map(str, Ns3)), "N_last": Ns3[-1],
                    "p_observed": p, "Q_last": Q[2], "Q_inf": Qinf,
                    "resid_last_pct": resid})
        print(f"{label:17} {key:>8} {'/'.join(map(str, Ns3)):>16} {p:6.3f} "
              f"{Qinf:10.4f} {resid:12.2f}%")

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0].keys()))
        w.writeheader()
        w.writerows(out)
    print(f"\n-> {OUT}")


if __name__ == "__main__":
    main()
