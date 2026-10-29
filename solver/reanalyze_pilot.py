# -*- coding: utf-8 -*-
"""
solver/reanalyze_pilot.py
=========================
Recalcule les diagnostics (D*, P*, K*, Z, balance, tourbillons) a partir des
champs sauvegardes (results/fields/<cid>.npz), sans relancer le solveur.

Utile pour rafraichir les scalaires apres un changement de post-traitement
(ex. remplacement de l'indicateur de symetrie -- sans objet pour la cavite
entrainee -- par la position du tourbillon principal).

Usage :
  python solver/reanalyze_pilot.py --Re 400 --N 128 --n-list 2,4,6,8 --alpha-list 0,0.10
  python solver/reanalyze_pilot.py --Re 400 --N 128          # tous les champs trouves
"""

from __future__ import annotations

import argparse
import csv
import glob
import json
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
_PROV_A = os.path.join(_ROOT, "provenance", "Lid_driven_ondulation")
for _p in (_HERE, _PROV_A):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from solver.combined_geometry import build_combined_geometry  # noqa: E402
from solver.run_case import SCALARS_DIR, FIELDS_DIR, case_id   # noqa: E402
from solver.run_pilot import interaction, KEYS                 # noqa: E402

from diagnostics import compute_diagnostics_fields  # noqa: E402  (Projet A)


def lid_profile_regularized(x, eps=0.02):
    x = np.asarray(x, dtype=np.float64)
    g = np.ones_like(x)
    lo = x < eps
    hi = x > 1.0 - eps
    tl = x[lo] / eps
    g[lo] = 3.0 * tl ** 2 - 2.0 * tl ** 3
    th = (1.0 - x[hi]) / eps
    g[hi] = 3.0 * th ** 2 - 2.0 * th ** 3
    return g


def reanalyze_one(path, a=0.5, b=0.5, H=1.0, eps_lid=0.02):
    d = np.load(path)
    N = int(d["N"])
    n = float(d["n"])
    alpha = float(d["alpha"])
    Nw = int(d["Nw"])
    Re = float(d["Re"])
    ux, uy, rho = d["ux"], d["uy"], d["rho"]
    U_ref = float(d["U_ref"])
    tau = float(d["tau"])
    cid = os.path.splitext(os.path.basename(path))[0]

    geom = build_combined_geometry(N, n, alpha, Nw, a, b, H)
    U_lid_norm = lid_profile_regularized(geom["xc"], eps_lid)
    diag = compute_diagnostics_fields(
        ux, uy, geom, U_lid_norm,
        case_id=cid, Re=Re, N=N, h_over_H=alpha, Nw=Nw,
        lid="regularized", eps_lid=eps_lid, U_ref=U_ref, tau=tau,
        Ma=float(np.sqrt(3.0) * U_ref), rho=rho)

    mv = diag.get("main_vortex") or {}
    return {
        "case_id": cid, "n": n, "alpha": alpha, "Nw": Nw, "Re": Re, "N": N,
        "status": "reanalyzed", "iters": None, "residual": None,
        "Ma": float(np.sqrt(3.0) * U_ref),
        "fluid_area": diag["fluid_area"],
        "D_star": diag["D_star"], "P_star": diag["P_star"], "K_star": diag["K_star"],
        "Z": diag["Z"], "balance": diag["balance"],
        "vortex_x": float(mv.get("x", np.nan)), "vortex_y": float(mv.get("y", np.nan)),
        "n_vortices": diag["n_vortices"],
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--Re", type=float, default=400)
    ap.add_argument("--N", type=int, default=128)
    ap.add_argument("--Nw", type=int, default=1)
    ap.add_argument("--n-list", type=str, default=None)
    ap.add_argument("--alpha-list", type=str, default=None)
    ap.add_argument("--a", type=float, default=0.5)
    ap.add_argument("--b", type=float, default=0.5)
    ap.add_argument("--H", type=float, default=1.0)
    args = ap.parse_args()

    n_filt = set(float(x) for x in args.n_list.split(",")) if args.n_list else None
    a_filt = set(float(x) for x in args.alpha_list.split(",")) if args.alpha_list else None

    paths = sorted(glob.glob(os.path.join(FIELDS_DIR, "se_n*_N%d.npz" % args.N)))
    res = {}
    for p in paths:
        base = os.path.basename(p)
        if "_Re%d_" % int(args.Re) not in base:
            continue
        d0 = np.load(p)
        n = float(d0["n"]); al = float(d0["alpha"])
        if n_filt is not None and n not in n_filt:
            continue
        if a_filt is not None and al not in a_filt:
            continue
        summ = reanalyze_one(p, a=args.a, b=args.b, H=args.H)
        with open(os.path.join(SCALARS_DIR, summ["case_id"] + ".json"), "w",
                  encoding="utf-8") as fh:
            json.dump(summ, fh, indent=2, ensure_ascii=False)
        res[(summ["n"], summ["alpha"])] = summ
        print("[re] %-28s D*=%.5f P*=%.5f vort=(%.3f,%.3f) nv=%d"
              % (summ["case_id"], summ["D_star"], summ["P_star"],
                 summ["vortex_x"], summ["vortex_y"], summ["n_vortices"]))

    n_list = sorted({k[0] for k in res})
    alpha_list = sorted({k[1] for k in res})

    tab = os.path.join(SCALARS_DIR, "pilot_Re%g_N%d.csv" % (args.Re, args.N))
    with open(tab, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["n", "alpha"] + KEYS)
        for (n, al), s in sorted(res.items()):
            w.writerow([n, al] + [s.get(k) for k in KEYS])

    inter = {"Re": args.Re, "N": args.N, "Nw": args.Nw,
             "n_list": n_list, "alpha_list": alpha_list, "pairs": {}}
    for i in range(len(n_list)):
        for j in range(i + 1, len(n_list)):
            na, nb = n_list[i], n_list[j]
            for al in alpha_list:
                if al == 0.0:
                    continue
                if (na, al) in res and (na, 0.0) in res and \
                   (nb, al) in res and (nb, 0.0) in res:
                    inter["pairs"]["n%d_n%d_a%.3f" % (na, nb, al)] = \
                        interaction(res, na, nb, al)
    iq_file = os.path.join(SCALARS_DIR, "interaction_Re%g_N%d.json" % (args.Re, args.N))
    with open(iq_file, "w", encoding="utf-8") as fh:
        json.dump(inter, fh, indent=2, ensure_ascii=False)
    print("\n[ok] table = %s" % tab)
    print("[ok] I_Q   = %s" % iq_file)
    return 0


if __name__ == "__main__":
    sys.exit(main())
