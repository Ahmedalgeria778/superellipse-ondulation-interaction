"""Collect the derived quantities cited in the revised manuscript.

The script reads the derived scalar files produced by the LBM solver; it does not
re-run any simulation.  Its only output is revision_data.json, which is used to
check every number written in the revised text.
"""
from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCALARS = ROOT / "results" / "scalars"
OUT = Path(__file__).resolve().parent / "revision_data.json"

KEEP = (
    "case_id", "Re", "N", "n", "alpha", "Nw", "status", "iters", "residual",
    "Ma", "fluid_area", "D_star", "P_star", "K_star", "balance", "n_vortices",
)


def key(case: dict) -> tuple:
    return (case["Re"], case["N"], case["alpha"], case["Nw"], case["n"])


def main() -> None:
    cases = []
    for path in sorted(SCALARS.glob("*.json")):
        try:
            data = json.loads(path.read_text(encoding="utf-8"))
        except (ValueError, OSError):
            continue
        if not all(k in data for k in ("Re", "N", "n", "alpha", "D_star")):
            continue
        cases.append({k: data[k] for k in KEEP if k in data})

    index = {key(c): c for c in cases if c.get("status") == "converged"}

    contrasts = {}
    for (reyn, size, alpha, nw, _), case in index.items():
        for left, right in ((2, 4), (4, 6), (6, 8)):
            a = index.get((reyn, size, alpha, nw, float(left)))
            b = index.get((reyn, size, alpha, nw, float(right)))
            if a is None or b is None:
                continue
            contrasts[(reyn, size, alpha, nw, left, right)] = {
                "I_D_pct": 100.0 * (b["D_star"] - a["D_star"]) / a["D_star"],
                "I_P_pct": 100.0 * (b["P_star"] - a["P_star"]) / a["P_star"],
                "area_pct": 100.0 * (b["fluid_area"] - a["fluid_area"]) / a["fluid_area"],
                "delta_IP_over_ID": (
                    (b["P_star"] - b["D_star"]) - (a["P_star"] - a["D_star"])
                ) / (b["D_star"] - a["D_star"]),
            }

    spans = {}
    by_transition = defaultdict(list)
    for (reyn, size, alpha, nw, left, right), value in contrasts.items():
        if alpha == 0.0 and nw == 1:
            by_transition[(reyn, left, right)].append((size, value["I_D_pct"]))
    for k, pairs in by_transition.items():
        pairs.sort()
        values = [v for _, v in pairs]
        sub = [v for n, v in pairs if n <= 384]
        spans[k] = {
            "resolutions": [n for n, _ in pairs],
            "values_pct": values,
            "full_span_pp": max(values) - min(values),
            "span_128_384_pp": max(sub) - min(sub) if sub else None,
        }

    new = [c for c in cases if c["N"] == 768 and c["Re"] in (400, 1000)]
    payload = {
        "n_cases": len(cases),
        "n_converged": len(index),
        "new_N768_cases": sorted(new, key=key),
        "contrasts": {"/".join(map(str, k)): v for k, v in contrasts.items()},
        "spans": {"/".join(map(str, k)): v for k, v in spans.items()},
    }
    OUT.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"wrote {OUT} ({len(cases)} scalar files, {len(index)} converged)")


if __name__ == "__main__":
    main()
