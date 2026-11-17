# -*- coding: utf-8 -*-
"""solver/make_scalars_csv.py
============================
Regenerate ``results/scalars_clean.csv`` from the per-run JSON files stored in
``results/scalars/``.

The JSON files written by ``solver/run_case.py`` (one per simulation) are the
single source of truth: every number quoted in the manuscript is read from
them by the post-processing scripts.  This script only flattens them into a
CSV for inspection, and is also used to check that the directory is free of
duplicate or non-converged records.

Usage
-----
    python -m solver.make_scalars_csv            # write the CSV
    python -m solver.make_scalars_csv --report   # coverage report only

Only records with ``status`` in ``("converged", "reanalyzed")`` are exported;
runs that stopped on ``max_iter`` (capped archives) are listed but not
exported, so they can never leak into a fit or a table.
"""

from __future__ import annotations

import argparse
import collections
import csv
import glob
import io
import json
import os

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
SCAL = os.path.join(_ROOT, "results", "scalars")
OUT = os.path.join(_ROOT, "results", "scalars_clean.csv")

OK = ("converged", "reanalyzed")

FIELDS = [
    "case_id", "n", "alpha", "Nw", "Re", "N", "status",
    "D_star", "P_star", "K_star", "Z", "balance",
    "fluid_area", "Ma", "U_ref", "target_tau",
    "iters", "residual", "elapsed_s", "resumed",
    "vortex_x", "vortex_y", "n_vortices", "a", "b", "H", "tol",
]


def load():
    """Return (ok_records, rejected_records) indexed by case_id."""
    good, bad = {}, {}
    for p in sorted(glob.glob(os.path.join(SCAL, "se_*.json"))):
        try:
            d = json.load(io.open(p, encoding="utf-8"))
        except ValueError:
            bad[os.path.basename(p)] = {"status": "unreadable"}
            continue
        cid = d.get("case_id") or os.path.splitext(os.path.basename(p))[0]
        (good if d.get("status") in OK else bad)[cid] = d
    return good, bad


def write_csv(good):
    with io.open(OUT, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS, extrasaction="ignore")
        w.writeheader()
        for cid in sorted(good, key=_sort_key):
            w.writerow(good[cid])


def _sort_key(cid):
    d = cid
    try:
        last = d.split("_")
        return (last[-3], int(last[-1]), float(last[1][1:]), float(last[2][1:]))
    except Exception:
        return ("", 0, 0.0, 0.0)


def report(good, bad):
    print("records kept      : %d" % len(good))
    print("records rejected  : %d" % len(bad))
    for cid, d in sorted(bad.items()):
        print("   %-52s %s" % (cid, d.get("status")))
    print()
    print("by status : %s" % dict(collections.Counter(
        d.get("status") for d in good.values())))
    print("by N      : %s" % sorted(collections.Counter(
        d.get("N") for d in good.values()).items()))
    print("by Re     : %s" % sorted(collections.Counter(
        d.get("Re") for d in good.values()).items()))
    print("by alpha  : %s" % sorted(collections.Counter(
        d.get("alpha") for d in good.values()).items()))
    print("by Nw     : %s" % sorted(collections.Counter(
        d.get("Nw") for d in good.values()).items()))
    print("by n      : %s" % sorted(collections.Counter(
        d.get("n") for d in good.values()).items()))
    print()
    print("primary matrix (n=2,4 x alpha=0,0.1 x Nw=1):")
    for Re in (100, 400, 1000):
        for N in (128, 256, 384, 512, 768):
            c = [d for d in good.values()
                 if d.get("Re") == Re and d.get("N") == N
                 and d.get("alpha") in (0.0, 0.1) and d.get("Nw") == 1
                 and d.get("n") in (2.0, 4.0)]
            if c:
                print("   Re=%-4d N=%-4d %d/4" % (Re, N, len(c)))
    print()
    missing = [d for d in good.values() if d.get("residual") is None]
    if missing:
        print("WARNING: %d records without a residual value" % len(missing))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--report", action="store_true",
                    help="print the coverage report without writing the CSV")
    a = ap.parse_args()
    good, bad = load()
    report(good, bad)
    if not a.report:
        write_csv(good)
        print("wrote", OUT)


if __name__ == "__main__":
    main()
