# -*- coding: utf-8 -*-
"""
analyze.py
==========
Comparative analysis of the campaign (read-only from results/scalar_table.json).

Produces in results/analysis/ :

  normalize.csv : for each WAVY case, normalization by the flat reference
                  F00 with the same lid (regularized) and same Re/N :
                    K*_n = K*/K*_F00 , P*_n = P*/P*_F00 , D*_n = D*/D*_F00 ,
                    K*_a = K*/A_f   , D*_a = D*/A_f   (intensities / fluid area)
                  -> P* is NOT divided by A_f (lid power, lid area is
                     constant) ; we report P*_n = P*/P*_F00.

  lid_reg.csv   : effect of lid regularization, same flat wall :
                  F00 (regularized) vs V0x (uniform) at given Re/N.
                    Delta P* = P*_reg/P*_unif - 1 , Delta D* = D*_reg/D*_unif - 1

  topology.csv  : vortex count occurrences as a function of (h, Nw, Re).

Usage :
  python analyze.py [--outdir results] [--N 256]
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def _load_table(outdir):
    with open(os.path.join(outdir, "scalar_table.json"), "r", encoding="utf-8") as fh:
        return json.load(fh)


def _idx(rows):
    return {(r["case"], int(r["N"])): r for r in rows}


def _af(r):
    return float(r["fluid_area"]) if r.get("fluid_area") else (
        1.0 - 0.5 * float(r.get("h_over_H") or 0.0))


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="results")
    ap.add_argument("--N", type=int, default=None,
                    help="restrict to a single resolution (default: all)")
    args = ap.parse_args(argv)

    rows = _load_table(args.outdir)
    by = _idx(rows)
    outdir = os.path.join(args.outdir, "analysis")
    os.makedirs(outdir, exist_ok=True)

    # ---- 1) normalization by F00 ------------------------------------
    norm = []
    for r in rows:
        if r.get("group") != "main":
            continue
        Re, N = int(r["Re"]), int(r["N"])
        if args.N and N != args.N:
            continue
        ref = by.get((f"F00-R{Re}", N))
        if ref is None or ref.get("status") != "converged":
            continue
        af = _af(r)
        k, p, d = float(r["K_star"]), float(r["P_star"]), float(r["D_star"])
        k0, p0, d0 = float(ref["K_star"]), float(ref["P_star"]), float(ref["D_star"])
        norm.append({
            "case": r["case"], "Re": Re, "N": N,
            "h_over_H": r["h_over_H"], "Nw": r["Nw"], "A_f": af,
            "K_star": k, "P_star": p, "D_star": d,
            "K_per_A": k / af, "D_per_A": d / af,
            "K_norm_F00": k / k0, "P_norm_F00": p / p0, "D_norm_F00": d / d0,
            "dP_pct_F00": 100.0 * (p / p0 - 1.0),
            "dD_pct_F00": 100.0 * (d / d0 - 1.0),
            "dK_pct_F00": 100.0 * (k / k0 - 1.0),
            "balance_pct": 100.0 * float(r.get("balance") or 0.0),
            "n_vortices": r.get("n_vortices"),
            "checks_ok": r.get("checks_ok"),
        })
    norm.sort(key=lambda x: (x["Re"], x["h_over_H"], x["Nw"]))
    _write(os.path.join(outdir, "normalize.csv"), norm)

    # ---- 2) effect of lid regularization ----------------------------
    lid = []
    vmap = {100: "V01-R100", 400: "V02-R400", 1000: "V03-R1000"}
    for Re, vid in vmap.items():
        for N in sorted({int(r["N"]) for r in rows}):
            if args.N and N != args.N:
                continue
            reg = by.get((f"F00-R{Re}", N))
            uni = by.get((vid, N))
            if not reg or not uni:
                continue
            if reg.get("status") != "converged" or uni.get("status") != "converged":
                continue
            p_r, p_u = float(reg["P_star"]), float(uni["P_star"])
            d_r, d_u = float(reg["D_star"]), float(uni["D_star"])
            k_r, k_u = float(reg["K_star"]), float(uni["K_star"])
            lid.append({
                "Re": Re, "N": N,
                "P_reg": p_r, "P_unif": p_u, "dP_reg_pct": 100.0 * (p_r / p_u - 1.0),
                "D_reg": d_r, "D_unif": d_u, "dD_reg_pct": 100.0 * (d_r / d_u - 1.0),
                "K_reg": k_r, "K_unif": k_u, "dK_reg_pct": 100.0 * (k_r / k_u - 1.0),
                "bal_reg_pct": 100.0 * float(reg.get("balance") or 0.0),
                "bal_unif_pct": 100.0 * float(uni.get("balance") or 0.0),
            })
    _write(os.path.join(outdir, "lid_reg.csv"), lid)

    # ---- 3) topology ------------------------------------------------
    topo = []
    for r in rows:
        if r.get("group") != "main":
            continue
        if args.N and int(r["N"]) != args.N:
            continue
        h = float(r["h_over_H"]); nw = int(r["Nw"])
        topo.append({
            "case": r["case"], "Re": int(r["Re"]), "N": int(r["N"]),
            "h_over_H": h, "Nw": nw, "hNw": h * nw,
            "n_vortices": r.get("n_vortices"),
            "main_x": (r.get("main_vortex") or {}).get("x"),
            "main_y": (r.get("main_vortex") or {}).get("y"),
            "main_psi": (r.get("main_vortex") or {}).get("psi"),
            "checks_ok": r.get("checks_ok"),
        })
    topo.sort(key=lambda x: (x["Re"], x["hNw"], x["Nw"]))
    _write(os.path.join(outdir, "topology.csv"), topo)

    # ---- console summary ---------------------------------------------
    print(f"analysis summary (outdir={outdir})")
    print(f"  normalize.csv : {len(norm)} wavy cases normalized by F00")
    for Re in (100, 400, 1000):
        sub = [x for x in norm if x["Re"] == Re]
        if sub:
            dP = [x["dP_pct_F00"] for x in sub]
            dD = [x["dD_pct_F00"] for x in sub]
            print(f"    Re={Re:4d}: dP*/F00 in [{min(dP):+.2f},{max(dP):+.2f}]%  "
                  f"dD*/F00 in [{min(dD):+.2f},{max(dD):+.2f}]%")
    print(f"  lid_reg.csv   : {len(lid)} regularized/uniform comparison points")
    for x in lid:
        print(f"    Re={x['Re']:4d} N={x['N']}: dP*={x['dP_reg_pct']:+.2f}%  "
              f"dD*={x['dD_reg_pct']:+.2f}%  (bal {x['bal_unif_pct']:+.1f} -> "
              f"{x['bal_reg_pct']:+.1f}%)")
    print(f"  topology.csv  : {len(topo)} cases")
    return 0


def _write(path, rows):
    if not rows:
        with open(path, "w", newline="", encoding="utf-8") as fh:
            fh.write("")
        return
    keys = list(rows[0].keys())
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow(r)


if __name__ == "__main__":
    raise SystemExit(main())
