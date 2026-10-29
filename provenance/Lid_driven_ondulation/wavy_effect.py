# -*- coding: utf-8 -*-
"""Relative geometric effect of the wavy bottom vs the flat bottom, at fixed N.

    Delta_Q(N) = (Q_wavy(N) - Q_flat(N)) / Q_flat(N) * 100 %

Central question (proof 1 + 3 of the demonstration): is Delta_Q(N) stable
as N grows? If so, the effect is measured, not a residual of discretization
error. If the amplitude or the sign changes with N, the conclusion is suspended.

Data sources:
  - results/scalar_table.json          (campaign, mainly N=256)
  - results/convergence/resolution.csv (convergence study, N=128/256/384 + 512)

Priority given to the 'consistent' variant of resolution.csv.
Output: results/analysis/wavy_effect.csv
"""

from __future__ import annotations

import csv
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

RES = os.path.join("results", "convergence", "resolution.csv")
TABLE = os.path.join("results", "scalar_table.json")
OUT = os.path.join("results", "analysis", "wavy_effect.csv")

# (wavy case, flat case) to compare at the same N and same Re
PAIRS = [
    ("G06-R400", "F00-R400"),
    ("G12-R400", "F00-R400"),
    ("G06-R100", "F00-R100"),
    ("G12-R100", "F00-R100"),
    ("G06-R1000", "F00-R1000"),
    ("G12-R1000", "F00-R1000"),
    ("G13-R400", "F00-R400"),
    ("G13-R1000", "F00-R1000"),
]
QUANT = ("P_star", "D_star", "K_star")


def load():
    rows = {}
    with open(TABLE, encoding="utf-8") as fh:
        for r in json.load(fh):
            rows[(r["case"], int(r["N"]))] = r
    if os.path.exists(RES):
        with open(RES, newline="", encoding="utf-8") as fh:
            for r in csv.DictReader(fh):
                if r.get("variant", "consistent") != "consistent":
                    continue
                rows[(r["case"], int(r["N"]))] = r
    return rows


def main():
    rows = load()
    out = []
    print(f"{'wavy case':12} {'N':>4} {'dP*%':>8} {'dD*%':>8} {'dK*%':>8} {'nv':>6}")
    for wavy, flat in PAIRS:
        Ns = sorted({n for (c, n) in rows if c == wavy} &
                    {n for (c, n) in rows if c == flat})
        for N in Ns:
            a, b = rows[(wavy, N)], rows[(flat, N)]
            rec = {"case_wavy": wavy, "case_flat": flat, "N": N,
                   "Re": int(a["Re"]), "n_vortices_wavy": a.get("n_vortices"),
                   "n_vortices_flat": b.get("n_vortices")}
            for q in QUANT:
                va, vb = float(a[q]), float(b[q])
                rec[f"{q}_wavy"] = va
                rec[f"{q}_flat"] = vb
                rec[f"d{q}_pct"] = 100.0 * (va - vb) / vb
            out.append(rec)
            print(f"{wavy:12} {N:4d} {rec['dP_star_pct']:8.3f} {rec['dD_star_pct']:8.3f} "
                  f"{rec['dK_star_pct']:8.3f} {rec['n_vortices_flat']}/{rec['n_vortices_wavy']:>2}")
        if not Ns:
            print(f"{wavy:12}    -   (no reference at the same N)")

    # stability of the effect in N (proof 1)
    print("\n--- stability of dP* in N (same case) ---")
    for wavy, _ in PAIRS:
        vals = [(r["N"], r["dP_star_pct"]) for r in out if r["case_wavy"] == wavy]
        vals.sort()
        if len(vals) >= 2:
            ns = "/".join(str(n) for n, _ in vals)
            vs = "/".join(f"{v:.3f}" for _, v in vals)
            spread = max(v for _, v in vals) - min(v for _, v in vals)
            print(f"  {wavy:12} N={ns:16} dP*={vs:24} amplitude={spread:.3f} pp")

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0].keys()))
        w.writeheader()
        w.writerows(out)
    print(f"\n-> {OUT}")


if __name__ == "__main__":
    main()
