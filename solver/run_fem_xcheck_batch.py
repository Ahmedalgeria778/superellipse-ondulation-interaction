# -*- coding: utf-8 -*-
"""Run the FEM superellipse cross-checks needed for the PoF revision.

Cases:
  Re=400, n=2 & n=4, nx=ny=256   (mesh-convergence of the FEM cross-check)
  Re=1000, n=2 & n=4, nx=ny=128  (second Reynolds number)
  Re=400, n=2 & n=4, nx=ny=128   (re-run: adds native P*/K*)
Each run also records native D*, P*, K*.
"""
import subprocess, sys, os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PY = sys.executable
DRV = os.path.join(ROOT, "solver", "run_fem_xcheck.py")

jobs = [
    (400, 2, 128), (400, 4, 128),
    (1000, 2, 128), (1000, 4, 128),
    (400, 2, 256), (400, 4, 256),
]

for Re, n, mesh in jobs:
    cmd = [PY, DRV, "--Re", str(Re), "--n", str(n), "--alpha", "0.0",
           "--nx", str(mesh), "--ny", str(mesh), "--N", "128"]
    print("RUN", cmd, flush=True)
    r = subprocess.run(cmd, cwd=ROOT)
    print("DONE rc=%d" % r.returncode, flush=True)
print("ALL DONE", flush=True)
