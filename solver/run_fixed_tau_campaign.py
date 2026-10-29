# -*- coding: utf-8 -*-
"""Regenerate the resolution/contrast study at a fixed target tau per Re.

Convention (TARGET_TAU in provenance/Lid_driven_ondulation/lbm_wavy.py):
    Re=100 -> tau=0.65, Re=400 -> tau=0.60, Re=1000 -> tau=0.54
with U = (tau-1/2)/3 * Re / N  (diffusive: U ~ 1/N, tau fixed).

This driver only runs the cases that are missing or inconsistent with the
target tau, so it is resumable.  It writes fields+scalars exactly like
run_case.py so the rest of the pipeline is unchanged.
"""
from __future__ import annotations

import json
import os
import sys
import time
import traceback

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
for _p in (_HERE, _ROOT, os.path.join(_ROOT, "provenance", "Lid_driven_ondulation")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from solver.run_case import run_one, case_id, tol_suffix  # noqa: E402

FIELDS_DIR = os.path.join(_ROOT, "results", "fields")
TARGET = {100: 0.65, 400: 0.60, 1000: 0.54}

# (Re, N, n, alpha)
JOBS = []
for _Re in (100, 400, 1000):
    for _N in (128, 256, 384):
        for _n in (2, 4):
            for _a in (0.0, 0.10):
                JOBS.append((_Re, _N, _n, _a))


def need_run(Re, N, n, alpha, tol=1e-8):
    import numpy as np
    fn = os.path.join(FIELDS_DIR, case_id(n, alpha, 1, Re, N, tol) + ".npz")
    if not os.path.exists(fn):
        return True, "missing"
    try:
        d = dict(np.load(fn))
        tau = float(d["tau"])
    except Exception:
        return True, "unreadable"
    if abs(tau - TARGET[Re]) > 1e-3:
        return True, "tau=%.4f != %.2f" % (tau, TARGET[Re])
    return False, "ok"


def main():
    log = open(os.path.join(_ROOT, "results", "logs",
                            "fixed_tau_campaign.log"), "w", encoding="utf-8")
    todo = []
    for (Re, N, n, a) in JOBS:
        do, why = need_run(Re, N, n, a)
        line = "Re=%d N=%d n=%d a=%.2f : %s (%s)" % (Re, N, n, a,
                                                     "RUN" if do else "skip", why)
        print(line)
        log.write(line + "\n")
        log.flush()
        if do:
            todo.append((Re, N, n, a))
    print("total to run:", len(todo))
    log.write("total to run: %d\n" % len(todo))
    log.flush()

    for (Re, N, n, a) in todo:
        t0 = time.time()
        tag = "Re=%d N=%d n=%d a=%.2f" % (Re, N, n, a)
        try:
            summ = run_one(n, a, 1, Re, N, tol=1e-8, max_iter=400000,
                           verbose=True, log_path=os.path.join(
                               _ROOT, "results", "logs",
                               case_id(n, a, 1, Re, N, 1e-8) + ".log"))
            line = "[DONE] %s  D*=%.5f  %.0fs" % (tag, summ["D_star"],
                                                  time.time() - t0)
        except Exception:
            line = "[FAIL] %s\n%s" % (tag, traceback.format_exc())
        print(line)
        log.write(line + "\n")
        log.flush()
    log.write("campaign finished\n")
    log.close()


if __name__ == "__main__":
    main()
