# -*- coding: utf-8 -*-
"""Campagne FEM P2.2 — Factorial (tau,U) a Re=100 et Re=1000.

16 runs LBM (2 Re x 2 N x 4 U x 2 n) pour montrer que la near-cancellation
de l'erreur commune dans la difference appariee n'est pas specifique a Re=400.
"""
import subprocess, sys, os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PY = sys.executable
DRV = os.path.join(ROOT, "solver", "run_ma_sweep.py")

# Configuration: Re x N x n x U-list
configs = [
    # Re=100, N=128
    (100, 128, 2, "0.013021,0.026042,0.052083,0.104167"),
    (100, 128, 4, "0.013021,0.026042,0.052083,0.104167"),
    # Re=100, N=256
    (100, 256, 2, "0.013021,0.026042,0.052083,0.104167"),
    (100, 256, 4, "0.013021,0.026042,0.052083,0.104167"),
    # Re=1000, N=128
    (1000, 128, 2, "0.034722,0.069444,0.138889,0.277778"),
    (1000, 128, 4, "00.034722,0.069444,0.138889,0.277778"),
    # Re=1000, N=256
    (1000, 256, 2, "0.034722,0.069444,0.138889,0.277778"),
    (1000, 256, 4, "0.034722,0.069444,0.138889,0.277778"),
]

print(f"[P2.2] Lancement de {len(configs)} sweeps factoriaux", flush=True)
for i, (Re, N, n, U_list) in enumerate(configs, 1):
    cmd = [PY, "-m", "solver.run_ma_sweep", "--Re", str(Re), "--N", str(N),
           "--n", str(n), "--alpha", "0.0", "--U-list", U_list,
           "--workers", "2", "--max-iter", "600000"]
    print(f"[{i}/{len(configs)}] Re={Re}, N={N}, n={n}, U-list={U_list}", flush=True)
    r = subprocess.run(cmd, cwd=ROOT)
    if r.returncode != 0:
        print(f"  ERREUR rc={r.returncode}", flush=True)
    else:
        print(f"  OK", flush=True)

print("[P2.2] TOUS LES SWEES SONT TERMINES", flush=True)
