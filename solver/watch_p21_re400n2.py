# -*- coding: utf-8 -*-
"""Relance le 4e cas P2.1 (Re=400, n=2, N=512) des qu'un slot se libere.

Le pilote de campagne N=512 a ete arrete volontairement apres le lancement des
quatre cas utiles : les deux cas Re=100 en N=512 ne doivent pas tourner
(Ma_local different, aucun interet pour le span de (2,4) a ce Reynolds nombre).
Les trois cas N=512 deja en vol ecrivent leurs JSON normalement.

Usage : python solver/watch_p21_re400n2.py
"""
from __future__ import annotations

import json
import os
import re
import subprocess
import sys
import time

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
SCAL = os.path.join(_ROOT, "results", "scalars")
LOG = os.path.join(_ROOT, "results", "logs", "watch_p21_re400n2.log")
PY = sys.executable

TARGET = "se_n2_w0.00_Nw1_Re400_N512_tol1e-08"
MAX_LIVE = 8


def log(msg):
    stamp = time.strftime("%H:%M:%S")
    line = "[%s] %s" % (stamp, msg)
    with open(LOG, "a", encoding="utf-8") as fh:
        fh.write(line + "\n")
    print(line, flush=True)


def live_cases():
    out = subprocess.run(
        ["powershell", "-NoProfile", "-Command",
         "Get-CimInstance Win32_Process -Filter \"Name='python.exe'\" | "
         "Select-Object -ExpandProperty CommandLine"],
        capture_output=True, text=True).stdout
    n = len(re.findall(r"-m solver\.run_case", out.replace("\r", "").replace("\n", " ")) )
    return n // 2 if n else 0


def target_done():
    for suf in ("_tol1e-08", ""):
        p = os.path.join(SCAL, TARGET + suf + ".json")
        if os.path.exists(p):
            try:
                with open(p) as fh:
                    return json.load(fh).get("status") == "converged"
            except (OSError, ValueError):
                return False
    return False


def launch():
    logpath = os.path.join(_ROOT, "results", "logs", TARGET + ".log")
    cmd = [PY, "-m", "solver.run_case", "--n", "2", "--alpha", "0", "--Nw", "1",
           "--Re", "400", "--N", "512", "--tol", "1e-08", "--max_iter", "4000000",
           "--log", logpath]
    flags = 0x00000008 | 0x00000200  # DETACHED_PROCESS | NEW_PROCESS_GROUP
    with open(logpath, "w", encoding="utf-8") as fh:
        subprocess.Popen(cmd, cwd=_ROOT, stdout=fh, stderr=subprocess.STDOUT,
                         creationflags=flags)
    log("LANCE %s" % TARGET)


def main():
    log("watchdog P2.1 demarre ( cible=%s )" % TARGET)
    while True:
        if target_done():
            log("JSON converge deja present : rien a lancer, sort")
            return 0
        n = live_cases()
        if n < MAX_LIVE:
            log("slots libres (%d cas en vol) : lancement" % n)
            launch()
            return 0
        log("machine saturee (%d cas en vol) : reessai dans 300 s" % n)
        time.sleep(300)


if __name__ == "__main__":
    os.environ.setdefault("NUMBA_NUM_THREADS", "2")
    os.environ.setdefault("OMP_NUM_THREADS", "2")
    sys.exit(main())
