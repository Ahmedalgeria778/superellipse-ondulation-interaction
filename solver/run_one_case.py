# -*- coding: utf-8 -*-
"""Run a single (Re, N, n, alpha) fixed-tau case. Used for parallel dispatch.

Usage:
    python run_one_case.py <Re> <N> <n> <alpha>
"""
from __future__ import annotations

import os
import sys
import time
import traceback

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
for _p in (_HERE, _ROOT, os.path.join(_ROOT, "provenance", "Lid_driven_ondulation")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from solver.run_case import run_one, case_id  # noqa: E402

LOG_DIR = os.path.join(_ROOT, "results", "logs")

def main():
    Re = int(sys.argv[1]); N = int(sys.argv[2])
    n = float(sys.argv[3]); a = float(sys.argv[4])
    tag = "Re=%d N=%d n=%g a=%.2f" % (Re, N, n, a)
    t0 = time.time()
    line = None
    try:
        summ = run_one(n, a, 1, Re, N, tol=1e-8, max_iter=400000, verbose=True,
                       log_path=os.path.join(
                           LOG_DIR, case_id(n, a, 1, Re, N, 1e-8) + ".log"))
        line = "[DONE] %s  D*=%.5f  %.0fs" % (tag, summ["D_star"], time.time() - t0)
    except Exception:
        line = "[FAIL] %s\n%s" % (tag, traceback.format_exc())
    with open(os.path.join(LOG_DIR, "parallel_done.log"), "a", encoding="utf-8") as fh:
        fh.write(line + "\n")

if __name__ == "__main__":
    main()
