# -*- coding: utf-8 -*-
"""
solver/resume_all.py
====================
Relance avec reprise sur checkpoint tous les cas de P2.1 et P3.2 qui ne sont
pas encore convergés. run_case recharge automatiquement results/ckpt/<cid>.npz.

Lancement par liste d'arguments (pas de shell) : un chemin contenant des espaces
ne peut donc plus etre tronque, contrairement a un --log colle dans une chaine.

Usage :
  python -m solver.resume_all --dry-run
  python -m solver.resume_all
"""

from __future__ import annotations

import glob
import io
import json
import os
import subprocess
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
LOGS = os.path.join(_ROOT, "results", "logs")
SCAL = os.path.join(_ROOT, "results", "scalars")

PY = sys.executable
FLAGS = 0x00000008 | 0x00000200  # DETACHED_PROCESS | NEW_PROCESS_GROUP

# (n, alpha, Nw, Re, N, max_iter)
CASES = [
    # P2.1 priorite 1 : Re=1000, n=6 et n=8
    (6, 0.00, 1, 1000, 512, 4_000_000),
    (8, 0.00, 1, 1000, 512, 4_000_000),
    # P2.1 priorite 2
    (2, 0.00, 1, 1000, 512, 4_000_000),
    (4, 0.00, 1, 1000, 512, 4_000_000),
    (4, 0.00, 1, 400, 512, 4_000_000),
    (2, 0.00, 1, 400, 512, 4_000_000),
    # P3.2 restants
    (2, 0.05, 1, 1000, 384, 2_000_000),
    (2, 0.20, 1, 1000, 384, 2_000_000),
]


def cid(n, alpha, Nw, Re, N, tol=1e-8):
    return "se_n%g_w%.2f_Nw%d_Re%d_N%d_tol1e-08" % (n, alpha, Nw, Re, N)


def converged(name):
    p = os.path.join(SCAL, name + ".json")
    if not os.path.exists(p):
        return False
    try:
        return json.load(io.open(p, encoding="utf-8")).get("status") == "converged"
    except ValueError:
        return False


def live():
    q = ("Get-CimInstance Win32_Process -Filter \"Name='python.exe'\" | "
         "Select-Object -ExpandProperty CommandLine")
    out = subprocess.run(["powershell", "-NoProfile", "-Command", q],
                         capture_output=True, text=True).stdout
    return out.replace("\r", " ").replace("\n", " ")


def main():
    dry = "--dry-run" in sys.argv
    running = live()
    for (n, alpha, Nw, Re, N, mi) in CASES:
        name = cid(n, alpha, Nw, Re, N)
        if converged(name):
            print("[skip] %s deja converge" % name)
            continue
        if name in running:
            print("[skip] %s deja en vol" % name)
            continue
        logp = os.path.join(LOGS, name + ".log")
        cmd = [PY, "-m", "solver.run_case", "--n", str(n), "--alpha", str(alpha),
               "--Nw", str(Nw), "--Re", str(Re), "--N", str(N), "--tol", "1e-08",
               "--max_iter", str(mi), "--log", logp]
        if "python.exe" not in PY:
            print("[ERREUR] interpreteur non venv : %s" % PY)
            return 2
        if dry:
            print("[dry ] %s" % name)
            continue
        # stdout vers le neant : run_case ecrit deja dans --log, deux handles sur
        # le meme fichier se recouvriraient et effaceraient la ligne [ckpt].
        subprocess.Popen(cmd, cwd=_ROOT, stdout=subprocess.DEVNULL,
                         stderr=subprocess.DEVNULL, creationflags=FLAGS)
        print("[lanc] %s" % name)
    return 0


if __name__ == "__main__":
    sys.exit(main())
