# -*- coding: utf-8 -*-
"""Orchestrateur serialise par priorite.

Attend la fin de P1.1 + P1.2 + P4.1 ( Marqueurs de fin dans leurs logs ),
puis lance en serie :
  1. P3.2  (4 runs LBM a =0.05/0.2, N=384,  2 workers)
  2. P3.1  (4 runs LBM N_w=2,        N=384,  2 workers)
  3. P2.1  (6 runs LBM N=512,             3 workers)
  4. P2.2  (verifie les sweeps restants,  2 workers)  [rattrapage si tue]

Usage :
  python solver/orchestrator_P2_P3.py
"""
from __future__ import annotations

import os
import re
import subprocess
import sys
import time

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
LOG_DIR = os.path.join(_ROOT, "results", "logs")

PY = sys.executable

# log -> regex de fin de campagne
DONE_MARKERS = {
    "campaign_P1.1.log":   r"\[campagne\] terminee",
    "campaign_P4.1.log":   r"\[campagne\] terminee",
}
# P1.2 : driver v1 tue puis remplacé par v2 -> fin = log v2 OU v1 (au cas où)
DONE_MARKERS_ALT = {
    "campaign_P1.2_fem_v2.log": r"CAMPAGNE FEM P1.2 TERMINEE",
    "campaign_P1.2_fem.log":    r"TOUS LES RUNS SONT TERMINES",
}


def _marker_present(logs_dir, name, pat):
    path = os.path.join(logs_dir, name)
    if not os.path.exists(path):
        return False
    try:
        with open(path, encoding="utf-8", errors="replace") as fh:
            return bool(re.search(pat, fh.read()))
    except OSError:
        return False


def all_done(logs_dir, markers, markers_alt):
    for name, pat in markers.items():
        if not _marker_present(logs_dir, name, pat):
            return False
    # au moins un des marqueurs alternatifs (P1.2 v1 ou v2)
    if not any(_marker_present(logs_dir, n, p) for n, p in markers_alt.items()):
        return False
    return True


def run_phase(name, cmd):
    log_path = os.path.join(LOG_DIR, "orchestrated_%s.log" % name)
    print("[orchestrator] LANCE %s -> %s" % (name, log_path), flush=True)
    t0 = time.time()
    with open(log_path, "w", encoding="utf-8") as log:
        rc = subprocess.call(cmd, cwd=_ROOT, stdout=log, stderr=subprocess.STDOUT)
    print("[orchestrator] %s termine rc=%d (%.1f min)"
          % (name, rc, (time.time() - t0) / 60.0), flush=True)
    return rc


def main():
    env_threads = None  # herite de l'environnement du lanceur

    print("[orchestrator] attente de P1.1 + P1.2 + P4.1...", flush=True)
    n_check = 0
    while not all_done(LOG_DIR, DONE_MARKERS, DONE_MARKERS_ALT):
        time.sleep(300)
        n_check += 1
        if n_check % 12 == 0:  # every ~1h
            print("[orchestrator] toujours en attente (%d verifs)" % n_check,
                  flush=True)

    # petite patience pour que les descripteurs de fichiers se liberent
    time.sleep(30)

    phases = [
        # P3.1/P3.2 lancés manuellement le 26/10 -> seul P2.1 reste à orchestrer
        ("P2.1_N512", [PY, "-m", "solver.run_campaign",
                       "--N", "512", "--Re-list", "100,400,1000",
                       "--n-list", "2,4", "--alpha-list", "0",
                       "--tol", "1e-8", "--workers", "3",
                       "--max_iter", "4000000"]),
    ]

    rc_total = 0
    for name, cmd in phases:
        rc_total |= run_phase(name, cmd)

    print("[orchestrator] PHASES P3/P2.1 TERMINEES rc_total=%d" % rc_total,
          flush=True)
    return rc_total


if __name__ == "__main__":
    os.environ.setdefault("NUMBA_NUM_THREADS", "2")
    os.environ.setdefault("OMP_NUM_THREADS", "2")
    sys.exit(main())
