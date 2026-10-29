# -*- coding: utf-8 -*-
"""
solver/run_pilot.py
===================
Execute la campagne de CETTE etude (geometrie combinee) a une resolution fixe,
en reutilisant les cas deja calcules, puis calcule la difference de differences

    I_Q = [Q(n_b, alpha) - Q(n_b, 0)] - [Q(n_a, alpha) - Q(n_a, 0)]

pour toutes les paires (n_a < n_b) et chaque alpha, sur Q in {D_star, P_star, K_star}.

Aucune etude de convergence ici : la resolution est une hypothese de travail
(validee dans les projets d'origine).

Usage :
  python solver/run_pilot.py --N 128 --Re 400 --n-list 2,4,6,8 --alpha-list 0,0.10
  python solver/run_pilot.py --N 128 --Re 400 --n-list 4,8 --alpha-list 0,0.10 --reuse
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
_PROV_A = os.path.join(_ROOT, "provenance", "Lid_driven_ondulation")
for _p in (_HERE, _PROV_A):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from solver.run_case import run_one, case_id, tol_suffix, SCALARS_DIR  # noqa: E402

KEYS = ["D_star", "P_star", "K_star", "balance", "status", "iters",
        "residual", "Ma", "fluid_area", "vortex_x", "vortex_y", "n_vortices"]


def get_or_run(n, alpha, Nw, Re, N, a, b, H, max_iter, reuse, tol=1e-9,
               log_path=None):
    cid = case_id(n, alpha, Nw, Re, N, tol)
    path = os.path.join(SCALARS_DIR, cid + ".json")
    if reuse and os.path.exists(path):
        with open(path, encoding="utf-8") as fh:
            cached = json.load(fh)
        if (cached.get("status") == "converged"
                and float(cached.get("tol", 1e-9)) == float(tol)):
            print("[reuse] %s (converged, tol=%g)" % (cid, tol))
            return cached
        print("[redo ] %s (cache status=%s tol=%s)" %
              (cid, cached.get("status"), cached.get("tol", 1e-9)))
    else:
        print("[new  ] %s" % cid)
    return run_one(n, alpha, Nw, Re, N, a=a, b=b, H=H, max_iter=max_iter,
                   tol=tol, log_path=log_path)


def interaction(res, n_a, n_b, alpha):
    out = {}
    for key in ("D_star", "P_star", "K_star"):
        da = res[(n_a, alpha)][key] - res[(n_a, 0.0)][key]
        db = res[(n_b, alpha)][key] - res[(n_b, 0.0)][key]
        base_a = res[(n_a, 0.0)][key]
        base_b = res[(n_b, 0.0)][key]
        out[key] = {
            "delta_a": da, "delta_b": db, "I": db - da,
            "I_rel_a_pct": 100.0 * (db - da) / base_a if base_a else None,
        }
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int, default=128)
    ap.add_argument("--Re", type=float, default=400)
    ap.add_argument("--Nw", type=int, default=1)
    ap.add_argument("--n-list", type=str, default="4,8")
    ap.add_argument("--alpha-list", type=str, default="0,0.10")
    ap.add_argument("--max_iter", type=int, default=200000)
    ap.add_argument("--tol", type=float, default=1e-9)
    ap.add_argument("--reuse", action="store_true")
    ap.add_argument("--a", type=float, default=0.5)
    ap.add_argument("--b", type=float, default=0.5)
    ap.add_argument("--H", type=float, default=1.0)
    args = ap.parse_args()

    n_list = [float(x) for x in args.n_list.split(",")]
    alpha_list = [float(x) for x in args.alpha_list.split(",")]

    res = {}
    for n in n_list:
        for al in alpha_list:
            res[(n, al)] = get_or_run(n, al, args.Nw, args.Re, args.N, args.a,
                                      args.b, args.H, args.max_iter, args.reuse,
                                      tol=args.tol)

    suf = tol_suffix(args.tol)
    os.makedirs(SCALARS_DIR, exist_ok=True)
    tab = os.path.join(SCALARS_DIR,
                       "pilot_Re%g_N%d%s.csv" % (args.Re, args.N, suf))
    with open(tab, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["n", "alpha"] + KEYS)
        for (n, al), s in sorted(res.items()):
            w.writerow([n, al] + [s.get(k) for k in KEYS])

    interaction_all = {"Re": args.Re, "N": args.N, "Nw": args.Nw, "tol": args.tol,
                       "n_list": n_list, "alpha_list": alpha_list, "pairs": {}}
    print("\n=== Interactions I_Q (Re=%g, N=%d, Nw=%d) ===" % (args.Re, args.N, args.Nw))
    for i in range(len(n_list)):
        for j in range(i + 1, len(n_list)):
            na, nb = n_list[i], n_list[j]
            for al in alpha_list:
                if al == 0.0:
                    continue
                iq = interaction(res, na, nb, al)
                interaction_all["pairs"]["n%d_n%d_a%.3f" % (na, nb, al)] = iq
                print("  n_a=%g n_b=%g alpha=%.2f" % (na, nb, al))
                for key in ("D_star", "P_star", "K_star"):
                    v = iq[key]
                    print("     %-7s dQ(na)=%+.5f dQ(nb)=%+.5f  I=%+.6f  (I/Dna=%.3f%%)"
                          % (key, v["delta_a"], v["delta_b"], v["I"],
                             v["I_rel_a_pct"] or 0.0))

    iq_file = os.path.join(SCALARS_DIR,
                           "interaction_Re%g_N%d%s.json" % (args.Re, args.N, suf))
    with open(iq_file, "w", encoding="utf-8") as fh:
        json.dump(interaction_all, fh, indent=2, ensure_ascii=False)
    print("\n[ok] table   = %s" % tab)
    print("[ok] I_Q     = %s" % iq_file)
    return 0


if __name__ == "__main__":
    sys.exit(main())
