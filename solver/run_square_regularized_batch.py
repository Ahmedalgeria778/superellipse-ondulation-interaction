# -*- coding: utf-8 -*-
"""
solver/run_square_regularized_batch.py
======================================
Launches the missing flat-square, REGULARIZED-lid Ghia-benchmark cases so that
the effect of the regularized lid alone can be isolated on a strictly flat
geometry (A_f = 1), at both N=128 and N=256.

Inventory BEFORE this batch (results/ghia/square_Re*_N*_{lid}.json):
    uniform     : Re=100/400/1000  x  N=128, N=256   (complete)
    regularized : Re=400 x N=256                     (present)
Missing (this batch):
    regularized : Re=100 N=128, Re=100 N=256,
                  Re=400 N=128,
                  Re=1000 N=128, Re=1000 N=256

Each case is run synchronously via solver.run_square_ghia.main()'s worker
(imported) to avoid spawning child processes.  Per-case logs go to
results/logs/square_Re{Re}_N{N}_regularized.log ; a batch log records
start/end and elapsed time.

Usage: python -m solver.run_square_regularized_batch
"""

from __future__ import annotations

import os
import sys
import time

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
for _p in (_HERE, _ROOT, os.path.join(_ROOT, "provenance", "Lid_driven_ondulation")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

# (Re, N) list of MISSING regularized cases
JOBS = [
    (400, 128),
    (100, 128),
    (100, 256),
    (1000, 128),
    (1000, 256),
]

GHIA = os.path.join(_ROOT, "results", "ghia")
LOGS = os.path.join(_ROOT, "results", "logs")
BATCH_LOG = os.path.join(LOGS, "square_regularized_batch.log")


def main():
    os.makedirs(GHIA, exist_ok=True)
    os.makedirs(LOGS, exist_ok=True)
    sys.argv = [sys.argv[0], "0", "0", "regularized"]  # placeholder
    import run_square_ghia as RS

    with open(BATCH_LOG, "a", encoding="utf-8") as bl:
        bl.write("=== batch start %s ===\n" % time.strftime("%Y-%m-%d %H:%M:%S"))
        for Re, N in JOBS:
            out = os.path.join(GHIA, "square_Re%d_N%d_regularized.json" % (Re, N))
            if os.path.exists(out):
                line = "[SKIP] Re=%d N=%d regularized already present" % (Re, N)
                print(line)
                bl.write(line + "\n")
                continue
            t0 = time.time()
            print("[START] Re=%d N=%d regularized" % (Re, N))
            RS.OUT = GHIA
            RS.LOG = LOGS
            sys.argv = ["run_square_ghia", str(Re), str(N), "regularized"]
            try:
                RS.main()
                line = "[OK]   Re=%d N=%d regularized  %.0fs" % (
                    Re, N, time.time() - t0)
            except Exception as exc:  # noqa: BLE001
                line = "[FAIL] Re=%d N=%d regularized  %s: %s" % (
                    Re, N, type(exc).__name__, exc)
            print(line)
            bl.write(line + "\n")
            bl.flush()
        bl.write("=== batch end %s ===\n" % time.strftime("%Y-%m-%d %H:%M:%S"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
