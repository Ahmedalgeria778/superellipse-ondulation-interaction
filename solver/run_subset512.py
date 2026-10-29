# -*- coding: utf-8 -*-
"""
solver/run_subset512.py
=======================
Campagne N=512 CIBLEE de 16 cas (recommandation du relecteur) plutot que les
152 cas de la campagne complete.

Selection (see docs/10_audit_contrastes_bilan.md, section 3.2) :

  Re    n            alpha          cas   question visee
  100   2, 4         0, 0.10         4   reduire le budget ~28 % sur I(2,4)
  400   2, 4, 6, 8   0, 0.10         8   tester les 3 intervalles 2->4/4->6/6->8
  1000  2, 4         0, 0.10         4   verifier I(2,4) a grand Reynolds

Conventions conservees a l'identique :
  * tolerance unique            : tol = 1e-8
  * loi de vitesse              : U = (40/3)/N  (SANS --U ; TARGET_TAU le
    fournit deja pour Re dans {100, 400, 1000})  ->  Ma = 0.0451 a N=512,
    tau = 0.90 / 0.60 / 0.54 pour Re = 100 / 400 / 1000
  * checkpoint                  : repris automatiquement par solver.run_case
    (ck.load -> s.f_init, sum[\"resumed\"]=True), jamais recalcule a zero.

Ordre de lancement : le checkpoint existant
  results/ckpt/se_n2_w0.00_Nw1_Re400_N512_tol1e-08.npz  (it = 200 000)
est soumis EN PREMIER, puis les cas les plus couteux d'abord (Re decroissant,
n decroissant) pour equilibrer le pool.

Usage :
  python -m solver.run_subset512 --dry-run
  python -m solver.run_subset512 --workers 8 --max_iter 4000000
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
from solver.run_campaign import status_of, LOG_DIR        # noqa: E402

# --- la selection des 16 cas ------------------------------------------------
CHECKPOINT_CASE = (400.0, 2.0, 0.0)          # repris en premier

SUBSET = (
    [(100.0, n, a) for n in (2.0, 4.0) for a in (0.0, 0.10)]
    + [(400.0, n, a) for n in (2.0, 4.0, 6.0, 8.0) for a in (0.0, 0.10)]
    + [(1000.0, n, a) for n in (2.0, 4.0) for a in (0.0, 0.10)]
)
assert len(SUBSET) == 16, len(SUBSET)
assert len(set(SUBSET)) == 16, "doublon dans la selection"


def ordered_cases():
    """Checkpoint d'abord, puis Re decroissant puis n decroissant."""
    rest = [c for c in SUBSET if c != CHECKPOINT_CASE]
    rest.sort(key=lambda c: (-c[0], -c[1], c[2]))
    return [CHECKPOINT_CASE] + rest


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int, default=512)
    ap.add_argument("--tol", type=float, default=1e-8)
    ap.add_argument("--workers", type=int, default=8)
    # estimation : 122k/270k/622k iterations a N=256 pour Re=100/400/1000,
    # x ~2.8 pour doubler N  ->  ~0.35M / 0.76M / 1.75M a N=512.
    ap.add_argument("--max_iter", type=int, default=4000000)
    ap.add_argument("--Nw", type=int, default=1)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    os.makedirs(LOG_DIR, exist_ok=True)
    ordered = ordered_cases()
    ck_cid = case_id(CHECKPOINT_CASE[1], CHECKPOINT_CASE[2], args.Nw,
                     CHECKPOINT_CASE[0], args.N, args.tol)

    pending, done = [], []
    for Re, n, al in ordered:
        cid = case_id(n, al, args.Nw, Re, args.N, args.tol)
        rec, _ = status_of(cid, args.tol)
        (done if rec else pending).append((Re, n, al, cid, rec))

    print("[subset512] N=%d tol=%g workers=%d max_iter=%d | %d cas, "
          "%d deja convergents, %d a lancer"
          % (args.N, args.tol, args.workers, args.max_iter,
             len(ordered), len(done), len(pending)), flush=True)
    for i, (Re, n, al) in enumerate(ordered, 1):
        cid = case_id(n, al, args.Nw, Re, args.N, args.tol)
        rec, _ = status_of(cid, args.tol)
        if rec:
            tag = "fait "
        elif cid == ck_cid:
            tag = "ckpt "
        else:
            tag = "lancer"
        print("  %2d. [%s] %-44s Re=%-5g n=%g alpha=%.2f"
              % (i, tag, cid, Re, n, al), flush=True)
    if args.dry_run or not pending:
        return 0

    exe = sys.executable
    t0 = time.time()
    failures = []

    def run_one(item):
        Re, n, al, cid, _ = item
        cmd = [exe, "-m", "solver.run_case", "--n", "%g" % n,
               "--alpha", "%g" % al, "--Nw", "%d" % args.Nw,
               "--Re", "%g" % Re, "--N", "%d" % args.N,
               "--tol", "%g" % args.tol,
               "--max_iter", "%d" % args.max_iter,
               "--log", os.path.join(LOG_DIR, cid + ".log")]
        t = time.time()
        rc = subprocess.call(cmd, cwd=_ROOT)
        return item, rc, time.time() - t

    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        futs = [ex.submit(run_one, it) for it in pending]
        for fu in as_completed(futs):
            (Re, n, al, cid, _), rc, dt = fu.result()
            rec, _ = status_of(cid, args.tol)
            ok = rec is not None
            if not ok:
                failures.append(cid)
            resumed = (rec or {}).get("resumed")
            print("[%s] %-44s Re=%-5g n=%g alpha=%.2f  %7.1f min  rc=%d "
                  "resum=%s iters=%s"
                  % ("OK  " if ok else "ECHEC", cid, Re, n, al, dt / 60.0,
                     rc, resumed, (rec or {}).get("iters")), flush=True)

    print("\n[subset512] termine en %.1f min ; %d echec(s)"
          % ((time.time() - t0) / 60.0, len(failures)), flush=True)
    for cid in failures:
        print("  ECHEC %s" % cid)
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
