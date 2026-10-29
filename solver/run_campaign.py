# -*- coding: utf-8 -*-
"""
solver/run_campaign.py
======================
Lanceur de campagne LBM en parallele borne (pool de sous-processus).

Un cas par processus (`solver.run_case`), `--workers` cas simultanes. Les
cas deja convergents a la bonne tolerance sont ignores (reuse). Les cas
les plus couteux (Re eleve) sont lances en premier pour equilibrer le pool.

Usage
-----
  python solver/run_campaign.py --N 256 --tol 1e-8 --workers 6
  python solver/run_campaign.py --N 256 --Re-list 100,400 --workers 8
  python solver/run_campaign.py --N 512 --Re-list 1000 --workers 4
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
for _p in (_HERE, _ROOT):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from solver.run_case import case_id, SCALARS_DIR          # noqa: E402

LOG_DIR = os.path.join(_ROOT, "results", "logs")


def status_of(cid, tol):
    """Un cas est rejoue s'il n'existe pas deja une solution acceptable.

    Reutilise tout enregistrement `converged` OU `reanalyzed` dont le
    residu stocke est deja <= tol. Les fichiers `reanalyzed` herites de
    l'ancienne campagne n'ont pas de residu : ils sont consideres valides
    (pas de relance d'un cas deja calcule).
    """
    path = os.path.join(SCALARS_DIR, cid + ".json")
    if not os.path.exists(path):
        return None, path
    try:
        with open(path, encoding="utf-8") as fh:
            rec = json.load(fh)
    except Exception:
        return None, path
    st = rec.get("status")
    if st in ("converged", "reanalyzed"):
        res = rec.get("residual")
        if res is None or float(res) <= float(tol):
            return rec, path
    return None, path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int, default=256)
    ap.add_argument("--Re-list", type=str, default="100,400,1000")
    ap.add_argument("--n-list", type=str, default="2,4,6,8")
    ap.add_argument("--alpha-list", type=str, default="0,0.10")
    ap.add_argument("--Nw", type=int, default=1)
    ap.add_argument("--tol", type=float, default=1e-8)
    ap.add_argument("--U", type=float, default=None,
                    help="vitesse de maquette imposee a tous les cas ; "
                         "U=(40/3)/N => Ma=0.1804 constant en Re. SANS cette "
                         "option le solveur retombe sur DEFAULT_TAU=0.60 pour "
                         "tout Re hors {100,400,1000}, ce qui fait croitre Ma "
                         "lineairement avec Re (invalide pour une etude en Re).")
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--max_iter", type=int, default=1200000)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    Res = [float(x) for x in args.Re_list.split(",")]
    ns = [float(x) for x in args.n_list.split(",")]
    alphas = [float(x) for x in args.alpha_list.split(",")]

    cases = [(n, al, Re) for Re in Res for n in ns for al in alphas]
    cases.sort(key=lambda c: (-c[2], -c[0], c[1]))

    os.makedirs(LOG_DIR, exist_ok=True)
    pending, done = [], []
    for n, al, Re in cases:
        cid = case_id(n, al, args.Nw, Re, args.N, args.tol)
        rec, _ = status_of(cid, args.tol)
        if rec is None and args.tol != 1e-9:
            # cas deja calcule a la tolerance par defaut (suffixe absent)
            rec, _ = status_of(case_id(n, al, args.Nw, Re, args.N, 1e-9),
                               args.tol)
        (done if rec else pending).append((n, al, Re, cid))

    print("[campagne] N=%d tol=%g workers=%d | %d deja convergents, "
          "%d a lancer" % (args.N, args.tol, args.workers, len(done),
                           len(pending)), flush=True)
    for n, al, Re, cid in pending:
        print("  a lancer : %-42s Re=%g n=%g alpha=%.2f" % (cid, Re, n, al),
              flush=True)
    if args.dry_run or not pending:
        return 0

    exe = sys.executable
    t0 = time.time()
    failures = []

    def run_one_case(item):
        n, al, Re, cid = item
        cmd = [exe, "-m", "solver.run_case", "--n", "%g" % n,
               "--alpha", "%g" % al, "--Nw", "%d" % args.Nw,
               "--Re", "%g" % Re, "--N", "%d" % args.N,
               "--tol", "%g" % args.tol,
               "--max_iter", "%d" % args.max_iter,
               "--log", os.path.join(LOG_DIR, cid + ".log")]
        if args.U is not None:
            cmd += ["--U", "%.12f" % args.U]
        t = time.time()
        rc = subprocess.call(cmd, cwd=_ROOT)
        return item, rc, time.time() - t

    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        futs = [ex.submit(run_one_case, it) for it in pending]
        for fu in as_completed(futs):
            (n, al, Re, cid), rc, dt = fu.result()
            ok = False
            if rc == 0:
                rec, _ = status_of(cid, args.tol)
                ok = rec is not None
            tag = "OK " if ok else "ECHEC"
            if not ok:
                failures.append(cid)
            print("[%s] %-42s Re=%g n=%g alpha=%.2f  %6.1f min  rc=%d"
                  % (tag, cid, Re, n, al, dt / 60.0, rc), flush=True)

    print("\n[campagne] terminee en %.1f min ; %d echec(s)"
          % ((time.time() - t0) / 60.0, len(failures)), flush=True)
    for cid in failures:
        print("  ECHEC %s" % cid)
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
