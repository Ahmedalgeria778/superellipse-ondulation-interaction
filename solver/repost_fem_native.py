# -*- coding: utf-8 -*-
"""solver/repost_fem_native.py
============================
Re-post-process a saved FEM native field (_natfield.npz) with the native
diagnostics, independently of the original run.  Used to audit the FEM
cross-check: it recomputes D*, P*, K* (native) exactly as
run_fem_xcheck.py does, from the saved mesh and solution vector.

Usage:
    python solver/repost_fem_native.py <path_to_*_natfield.npz> [...]
"""
from __future__ import annotations

import argparse
import glob
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
_FEM = os.path.join(_ROOT, "fem")
_PROV = os.path.join(_ROOT, "provenance", "Lid_driven_ondulation")
for _p in (_HERE, _FEM, _PROV):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import run_fem_xcheck as RX          # native fns + install_superellipse


def repost(path):
    z = np.load(path, allow_pickle=True)
    mesh = {"nodes": z["nodes"], "elems": z["elems"]}
    U = z["U"]
    M = int(z["M"])
    n = float(z["n"])
    alpha = float(z["alpha"])
    eps = float(z["eps_lid"])
    RX.install_superellipse(n, alpha, 1)
    D = RX.fem_native_dissipation(mesh, U, M)
    P = RX.fem_native_power(mesh, U, M, eps=eps)
    K = RX.fem_native_kinetic(mesh, U, M)
    area = RX.fem_mesh_area(mesh)
    print("%s" % os.path.basename(path))
    print("  n=%g alpha=%g nx=%d ny=%d Re=%g" %
          (n, alpha, int(z["nx"]), int(z["ny"]), float(z["Re"])))
    print("  A_mesh=%.8f  K*=%.8f  P*(nat)=%.8f  D*(nat)=%.8f"
          % (area, K, P, D))
    print("  native balance P/D-1 = %+.4f%%" % (100.0 * (P - D) / D))
    return D, P, K, area


def main(argv=None):
    ap = argparse.ArgumentParser(description="Re-post-process FEM native field")
    ap.add_argument("paths", nargs="+", help="*_natfield.npz paths or globs")
    args = ap.parse_args(argv)
    files = []
    for p in args.paths:
        matches = glob.glob(p)
        files.extend(matches if matches else [p])
    for f in files:
        if os.path.exists(f):
            repost(f)
            print()
        else:
            print("MISSING: %s" % f)


if __name__ == "__main__":
    main()
