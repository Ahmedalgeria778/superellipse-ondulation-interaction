# -*- coding: utf-8 -*-
"""
solver/run_case.py
==================
Execute un cas combine (superellipse + ondulation), sauvegarde :
  results/scalars/<case_id>.json   (diagnostics + controles)
  results/fields/<case_id>.npz     (ux, uy, rho, psi + masque)

Usage :
  python solver/run_case.py A                 # cas pilote A (via cases.py)
  python solver/run_case.py --n 4 --alpha 0.10 --Nw 1 --Re 400 --N 256
  python solver/run_case.py --all             # les 4 cas pilotes
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
_PROV_A = os.path.join(_ROOT, "provenance", "Lid_driven_ondulation")
for _p in (_HERE, _PROV_A):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from solver.lbm_combined import CombinedMRT  # noqa: E402
from solver import checks                     # noqa: E402

SCALARS_DIR = os.path.join(_ROOT, "results", "scalars")
FIELDS_DIR = os.path.join(_ROOT, "results", "fields")


def tol_suffix(tol):
    return "" if float(tol) == 1e-9 else "_tol%.0e" % float(tol)


def case_id(n, alpha, Nw, Re, N, tol=1e-9):
    return ("se_n%g_w%.2f_Nw%d_Re%d_N%d" % (n, alpha, Nw, int(Re), int(N))
            + tol_suffix(tol))


def run_one(n, alpha, Nw, Re, N, a=0.5, b=0.5, H=1.0, max_iter=300000,
            tol=1e-9, lid="regularized", verbose=True, save=True,
            log_path=None, ckpt=True, ckpt_every=200000, U_ref=None):
    os.makedirs(SCALARS_DIR, exist_ok=True)
    os.makedirs(FIELDS_DIR, exist_ok=True)
    cid = case_id(n, alpha, Nw, Re, N, tol)
    ck = None
    ckpt_path = None
    resumed = False
    if ckpt:
        import ckpt as _ck
        ck = _ck
        ckpt_path = os.path.join(ck.CKPT_DIR, cid + ".npz")
    t0 = time.time()
    log_fh = open(log_path, "w", encoding="utf-8") if log_path else None
    try:
        s = CombinedMRT(Re=Re, N=N, n=n, alpha=alpha, Nw=Nw, a=a, b=b, H=H,
                        lid=lid, max_iter=max_iter, tol=tol, verbose=verbose,
                        log=log_fh, case_id=cid, U_ref=U_ref)
        if ck is not None:
            ck.install(every=ckpt_every)
            f0 = ck.load(ckpt_path, cid)
            if f0 is not None:
                s.f_init = f0
                resumed = True
            ck.activate(cid, ckpt_path)
        try:
            s.run()
        finally:
            if ck is not None:
                ck.deactivate()
    finally:
        if log_fh:
            log_fh.close()
    diag = s.diagnostics()
    summ = checks.summarize(s, diag)
    summ["elapsed_s"] = time.time() - t0
    summ["a"], summ["b"], summ["H"] = a, b, H
    summ["tol"] = float(tol)
    summ["U_ref"] = float(s.U_ref)
    summ["target_tau"] = float(getattr(s, "target_tau", s.tau))
    summ["resumed"] = bool(resumed)
    print("[run_case] %s  D*=%.5f P*=%.5f K*=%.6f bal=%.3f%%  it=%d res=%.2e  %.1fs"
          % (cid, summ["D_star"], summ["P_star"], summ["K_star"],
             100 * summ["balance"], summ["iters"], summ["residual"], summ["elapsed_s"]))
    if save:
        with open(os.path.join(SCALARS_DIR, cid + ".json"), "w", encoding="utf-8") as fh:
            json.dump(summ, fh, indent=2, ensure_ascii=False)
        np.savez_compressed(
            os.path.join(FIELDS_DIR, cid + ".npz"),
            ux=s.ux, uy=s.uy, rho=s.rho,
            fluid=s.fluid, xc=s.geom["xc"], yc=s.geom["yc"],
            n=n, alpha=alpha, Nw=Nw, Re=Re, N=N, U_ref=s.U_ref, tau=s.tau,
        )
    if ck is not None and summ.get("status") == "converged":
        ck.clear(ckpt_path)
    return summ


def run_pilot(max_iter=300000, N=256, alpha_list=(0.0, 0.10)):
    from solver.cases import pilot_cases
    out = []
    for c in pilot_cases():
        if c.alpha not in alpha_list:
            continue
        out.append(run_one(c.n, c.alpha, c.Nw, c.Re, N, a=c.a, b=c.b, H=c.H,
                           max_iter=max_iter))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("case", nargs="?", help="nom du cas pilote (A,B,C,D) ou case_id")
    ap.add_argument("--n", type=float)
    ap.add_argument("--alpha", type=float, default=0.0)
    ap.add_argument("--Nw", type=int, default=1)
    ap.add_argument("--Re", type=float, default=400)
    ap.add_argument("--N", type=int, default=256)
    ap.add_argument("--max_iter", type=int, default=300000)
    ap.add_argument("--tol", type=float, default=1e-9)
    ap.add_argument("--U", type=float, default=None,
                    help="vitesse de maquette imposee (sinon tau du solveur) ; "
                         "U=(40/3)/N donne Ma=0.1804, constant en Re")
    ap.add_argument("--log", type=str, default=None)
    ap.add_argument("--all", action="store_true")
    args = ap.parse_args()

    if args.all:
        run_pilot(max_iter=args.max_iter, N=args.N)
        return 0
    if args.case and args.n is None:
        from solver.cases import by_name
        c = by_name(args.case)
        run_one(c.n, c.alpha, c.Nw, c.Re, args.N or c.N, a=c.a, b=c.b, H=c.H,
                max_iter=args.max_iter, tol=args.tol, log_path=args.log,
                U_ref=args.U)
        return 0
    if args.n is None:
        ap.error("precisez un cas (A/B/C/D) ou --n")
    run_one(args.n, args.alpha, args.Nw, args.Re, args.N,
            max_iter=args.max_iter, tol=args.tol, log_path=args.log,
            U_ref=args.U)
    return 0


if __name__ == "__main__":
    sys.exit(main())
