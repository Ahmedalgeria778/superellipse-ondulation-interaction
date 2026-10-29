# -*- coding: utf-8 -*-
"""Campagne FEM P1.2 — Re=100 et Re=1000, n=2 et n=4, nx=64/128/256.

12 runs FEM Taylor-Hood P2/P1 pour fermer la limitation
"FEM verification only at Re=400".
"""
import subprocess, sys, os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PY = sys.executable
DRV = os.path.join(ROOT, "solver", "run_fem_xcheck.py")

# Configuration: Re x n x nx
jobs = [
    # Re=100
    (100, 2, 64), (100, 2, 128), (100, 2, 256),
    (100, 4, 64), (100, 4, 128), (100, 4, 256),
    # Re=1000
    (1000, 2, 64), (1000, 2, 128), (1000, 2, 256),
    (1000, 4, 64), (1000, 4, 128), (1000, 4, 256),
]

print(f"[P1.2] Lancement de {len(jobs)} runs FEM", flush=True)
for i, (Re, n, nx) in enumerate(jobs, 1):
    cmd = [PY, DRV, "--Re", str(Re), "--n", str(n), "--alpha", "0.0",
           "--nx", str(nx), "--ny", str(nx), "--N", "128"]
    print(f"[{i}/{len(jobs)}] Re={Re}, n={n}, nx={nx}", flush=True)
    r = subprocess.run(cmd, cwd=ROOT)
    if r.returncode != 0:
        print(f"  ERREUR rc={r.returncode}", flush=True)
    else:
        print(f"  OK", flush=True)

print("[P1.2] TOUS LES RUNS SONT TERMINES", flush=True)
