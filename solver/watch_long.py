# -*- coding: utf-8 -*-
"""solver/watch_long.py — superviseur de run long (LBM, N=768).

But
---
Lancer un cas LBM en mode << aucun plafond >> (`tol=1e-8`, `max_iter=1e8`
dans `run_one_case_long.py`), et le relancer automatiquement s'il meurt
(coupure, kill, erreur transitoire).  Le superviseur s'arrete des que le
fichier scalaire du cas porte `status == "converged"`.

Chaque tentative ecrase `results/logs/<case>.log` (run_case ouvre le log en
"w") : le superviseur archive donc le log de la tentative precedente dans
`<case>.history.log` avant de relancer, pour que la chronologie survive.

Usage
-----
    python solver/watch_long.py <Re> <N> <n> <alpha>

Exemple :
    python solver/watch_long.py 400 768 2 0.0
"""
from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import time

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
LOGS = os.path.join(_ROOT, "results", "logs")
SCAL = os.path.join(_ROOT, "results", "scalars")
PY = sys.executable

RETRY_SLEEP = 120          # s entre deux tentatives
MAX_ATTEMPTS = 50          # garde-fou, jamais atteint en pratique


def case_id(Re, N, n, alpha):
    return "se_n%g_w%.2f_Nw1_Re%g_N%d_tol1e-08" % (n, alpha, Re, N)


def log(path, msg):
    line = "[%s] %s" % (time.strftime("%Y-%m-%d %H:%M:%S"), msg)
    with open(path, "a", encoding="utf-8") as fh:
        fh.write(line + "\n")
    print(line, flush=True)


def converged(cid):
    p = os.path.join(SCAL, cid + ".json")
    if not os.path.exists(p):
        return False
    try:
        with open(p, encoding="utf-8") as fh:
            return json.load(fh).get("status") == "converged"
    except (OSError, ValueError):
        return False


def rotate(path, history):
    if not os.path.exists(path):
        return
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as fh:
            body = fh.read()
    except OSError:
        return
    if not body.strip():
        return
    with open(history, "a", encoding="utf-8") as fh:
        fh.write("\n===== tentative precedente %s =====\n"
                 % time.strftime("%Y-%m-%d %H:%M:%S"))
        fh.write(body)


def main(argv):
    if len(argv) != 5:
        raise SystemExit("usage: watch_long.py <Re> <N> <n> <alpha>")
    Re, N, n, alpha = int(argv[1]), int(argv[2]), float(argv[3]), float(argv[4])
    cid = case_id(Re, N, n, alpha)
    os.makedirs(LOGS, exist_ok=True)
    wlog = os.path.join(LOGS, "watch_Re%d_N%d_n%g.log" % (Re, N, n))
    runlog = os.path.join(LOGS, cid + ".log")
    history = os.path.join(LOGS, cid + ".history.log")
    stdout = os.path.join(LOGS, cid + ".stdout.log")

    log(wlog, "superviseur demarre cible=%s (max_iter=1e8, tol=1e-8)" % cid)
    for attempt in range(1, MAX_ATTEMPTS + 1):
        if converged(cid):
            log(wlog, "JSON converge present : arret")
            return 0
        rotate(runlog, history)
        cmd = [PY, "-X", "utf8",
               os.path.join(_HERE, "run_one_case_long.py"),
               str(Re), str(N), str(n), str(alpha)]
        log(wlog, "tentative %d : %s" % (attempt, " ".join(cmd[2:])))
        with open(stdout, "a", encoding="utf-8") as out:
            out.write("\n===== tentative %d %s =====\n"
                      % (attempt, time.strftime("%Y-%m-%d %H:%M:%S")))
            out.flush()
            rc = subprocess.call(cmd, cwd=_ROOT, stdout=out,
                                 stderr=subprocess.STDOUT)
        log(wlog, "tentative %d : code de sortie %s" % (attempt, rc))
        if converged(cid):
            log(wlog, "converge : arret")
            return 0
        if attempt < MAX_ATTEMPTS:
            log(wlog, "relance dans %d s" % RETRY_SLEEP)
            time.sleep(RETRY_SLEEP)
    log(wlog, "abandon apres %d tentatives" % MAX_ATTEMPTS)
    return 1


if __name__ == "__main__":
    os.environ.setdefault("NUMBA_NUM_THREADS", "2")
    os.environ.setdefault("OMP_NUM_THREADS", "2")
    sys.exit(main(sys.argv))
