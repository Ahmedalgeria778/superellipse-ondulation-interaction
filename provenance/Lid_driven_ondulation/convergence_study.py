# -*- coding: utf-8 -*-
"""
convergence_study.py
====================
Targeted resolution convergence study (does not replace the campaign).

Two cases in two geometric regimes:
  - G06-R400: moderate corrugation (h/H=0.05, Nw=2, h*Nw=0.10);
  - G12-R400: strong corrugation   (h/H=0.075, Nw=4, h*Nw=0.30).

Resolutions: N in {128, 256, 384}.

Two variants to SEPARATE grid error and compressibility effect:

  A) 'consistent' : campaign scale -- target tau fixed per Re
     (TARGET_TAU), so U_ref ~ 1/N and Ma ~ 1/N (Ma decreases upon refinement).
  B) 'ma_fixed'   : Mach control -- U_ref frozen at the campaign value at
     N=256 (so Ma constant), Re frozen; tau = 3*U_ref*N/Re + 1/2 varies with N.

Output: results/convergence/resolution.csv (one row per case x variant x N)
and observed order estimation (Richardson at 3 levels) for K*, P*, D*.

Usage:
  python convergence_study.py [--N 128 256 384] [--max_iter 1000000] [--tol 1e-8]
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from diagnostics import compute_diagnostics
from checks import run_checks
from lbm_wavy import WavyMRT, TARGET_TAU

# (case, h/H, Nw, Re)
STUDY_CASES = [
    ("G06-R400", 0.05, 2, 400),
    ("G12-R400", 0.075, 4, 400),
]
REF_N = 256  # reference resolution for the Ma value in variant B


def uref_campaign(Re, N):
    """Campaign U_ref (fixed target tau) at resolution N."""
    tau = TARGET_TAU.get(Re, 0.60)
    return (tau - 0.5) / 3.0 * Re / N


def run_one(case, h, nw, Re, N, variant, outdir, max_iter, tol,
            lid="regularized", label=None, log=None):
    """One run of the convergence study.

    lid   : lid condition (standard 'regularized';
            'smoothstep_mean1' = profile renormalized to mean 1.0, protocol B)
    label : variant tag written to the CSV (default = variant)
    log   : optional progress stream passed to the solver (default None)
    """
    label = label or variant
    if variant == "consistent":
        kw = dict(target_tau=TARGET_TAU.get(Re, 0.60))
    else:  # ma_fixed : U_ref frozen -> Ma constant
        kw = dict(U_ref=uref_campaign(Re, REF_N))
    t0 = time.time()
    s = WavyMRT(Re, N, h, nw, lid=lid, eps_lid=0.02, max_iter=max_iter, tol=tol,
                verbose=True, log=log, **kw)
    s.case_id = f"{case}-{label}"
    s.run()
    diag = compute_diagnostics(s)
    diag["status"] = s.status
    diag["case"] = f"{case}-{label}"
    diag["N"] = N
    try:
        diag["checks"] = run_checks(s, diag)
        diag["checks_ok"] = bool(diag["checks"]["checks_ok"])
    except Exception:                                     # noqa: BLE001
        diag["checks_ok"] = None
    # lightweight save (traceability)
    od = os.path.join(outdir, label, f"{case}_N{N}")
    os.makedirs(od, exist_ok=True)
    from diagnostics import streamfunction
    np.savez_compressed(os.path.join(od, "fields.npz"),
                        ux=s.ux, uy=s.uy, rho=s.rho,
                        psi=streamfunction(s.ux, s.geom),
                        xc=s.geom["xc"], yc=s.geom["yc"], yb=s.geom["yb"],
                        fluid=s.geom["fluid"], N=N, Re=Re,
                        h_over_H=h, Nw=nw, U_ref=s.U_ref, tau=s.tau,
                        Ma=s.Ma, status=s.status, iters=s.iters,
                        last_residual=s.last_residual)
    row = {
        "case": case, "variant": label, "Re": Re, "N": N,
        "h_over_H": h, "Nw": nw, "hNw": h * nw,
        "status": s.status, "checks_ok": diag.get("checks_ok"),
        "Ma": float(s.Ma), "tau": float(s.tau), "U_ref": float(s.U_ref),
        "iters": int(s.iters), "resid": float(s.last_residual),
        "elapsed_s": round(time.time() - t0, 1),
        "K_star": float(diag["K_star"]), "P_star": float(diag["P_star"]),
        "D_star": float(diag["D_star"]), "balance_pct": 100.0 * float(diag["balance"]),
        "n_vortices": int(diag["n_vortices"]),
        "main_x": (diag["main_vortex"] or {}).get("x"),
        "main_y": (diag["main_vortex"] or {}).get("y"),
        "main_psi": (diag["main_vortex"] or {}).get("psi"),
    }
    print(f"  {case} {label:18s} N={N:4d}: Ma={s.Ma:.4f} tau={s.tau:.4f} "
          f"K*={row['K_star']:.6f} P*={row['P_star']:.4f} D*={row['D_star']:.4f} "
          f"bal={row['balance_pct']:+.1f}% nv={row['n_vortices']} "
          f"({row['elapsed_s']}s, {s.status})")
    return row


def _order(N1, Q1, N2, Q2, N3, Q3):
    """Observed order p such that Q_N = Q_exact - C N^-p (3 levels)."""
    d12, d23 = Q1 - Q2, Q2 - Q3
    if abs(d23) < 1e-300 or d12 * d23 <= 0:
        return float("nan")
    target = d12 / d23

    def f(p):
        return (N1 ** -p - N2 ** -p) / (N2 ** -p - N3 ** -p) - target

    lo, hi = 1e-3, 8.0
    flo, fhi = f(lo), f(hi)
    if flo * fhi > 0:
        return float("nan")
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        if f(lo) * f(mid) <= 0:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


CSV_KEYS = ["case", "variant", "Re", "N", "h_over_H", "Nw", "hNw", "status",
            "checks_ok", "Ma", "tau", "U_ref", "iters", "resid", "elapsed_s",
            "K_star", "P_star", "D_star", "balance_pct", "n_vortices",
            "main_x", "main_y", "main_psi"]


def write_csv(rows, out):
    with open(out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=CSV_KEYS, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int, nargs="+", default=[128, 256, 384])
    ap.add_argument("--outdir", default=os.path.join("results", "convergence"))
    ap.add_argument("--max_iter", type=int, default=1000000)
    ap.add_argument("--tol", type=float, default=1e-8)
    args = ap.parse_args(argv)

    os.makedirs(args.outdir, exist_ok=True)
    out = os.path.join(args.outdir, "resolution.csv")
    rows = []
    print(f"=== Convergence study (N={args.N}) ===")
    for (case, h, nw, Re) in STUDY_CASES:
        for variant in ("consistent", "ma_fixed"):
            for N in sorted(args.N):
                try:
                    row = run_one(case, h, nw, Re, N, variant,
                                  args.outdir, args.max_iter, args.tol)
                except Exception as exc:                  # noqa: BLE001
                    print(f"  FAILURE {case} {variant} N={N}: {exc}")
                    row = {"case": case, "variant": variant, "N": N,
                           "status": "error", "error": str(exc)}
                rows.append(row)
                write_csv(rows, out)                       # incremental flush
                print(f"  [checkpoint] {out}: {len(rows)} rows")
    print(f"\n-> {out}")

    # observed order ('consistent' variant, 3 levels)
    print("\n=== observed order (Richardson 3 levels, consistent variant) ===")
    for (case, _, _, _) in STUDY_CASES:
        sub = [r for r in rows if r["case"] == case
               and r["variant"] == "consistent" and r.get("status") == "converged"]
        ns = sorted(int(r["N"]) for r in sub)
        if len(ns) < 3:
            print(f"  {case}: insufficient resolutions ({ns})")
            continue
        r1 = next(r for r in sub if int(r["N"]) == ns[0])
        r2 = next(r for r in sub if int(r["N"]) == ns[1])
        r3 = next(r for r in sub if int(r["N"]) == ns[2])
        for q in ("K_star", "P_star", "D_star"):
            p = _order(ns[0], r1[q], ns[1], r2[q], ns[2], r3[q])
            print(f"  {case} {q:7s}: p ~ {p:.2f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
