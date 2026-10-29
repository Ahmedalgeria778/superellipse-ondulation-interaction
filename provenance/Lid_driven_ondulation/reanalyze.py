# -*- coding: utf-8 -*-
"""
reanalyze.py
============
Re-analysis of all already computed LBM/FEM fields (results/<id>/<id>_N<N>.npz) :
recomputes K*, P*, D*, balance, vortex topology with the current
settings and rebuilds results/scalar_table.json consistently.

Usage :
  python reanalyze.py [--outdir results]
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from geometry import build_geometry, lid_profile_regularized
from diagnostics import compute_diagnostics_fields, streamfunction

VALID_DIRNAME = "invalid"  # statuses != converged stored in results/invalid/


def analyze_npz(path):
    d = np.load(path)
    cid = str(d["case_id"])
    N = int(d["N"])
    Re = int(d["Re"])
    status = str(d["status"]) if "status" in d else None
    meta = {}
    base = os.path.splitext(path)[0]
    jpath = base.replace("_invalid", "") + ".json"
    if not os.path.exists(jpath):
        jpath = base + ".json"
    if os.path.exists(jpath):
        with open(jpath, "r", encoding="utf-8") as fh:
            meta = json.load(fh)
    if status is None:
        status = meta.get("status")
    if status is None:
        status = "converged"  # files predating the status check
    h = float(d["h_over_H"]) if "h_over_H" in d else float(meta.get("h_over_H", 0.0))
    nw = int(d["Nw"]) if "Nw" in d else int(meta.get("Nw", 1))
    lid = str(d["lid"]) if "lid" in d else str(meta.get("lid", "regularized"))
    eps = float(d["eps_lid"]) if "eps_lid" in d else float(meta.get("eps_lid", 0.02))
    if "flat_height" in d:
        fh = float(d["flat_height"])
        fh = None if fh == 0.0 else fh
    else:
        fh = meta.get("flat_height")
    U_ref = float(d["U_ref"])
    if "tau" in d:
        tau = float(d["tau"])
        Ma = float(d["U_ref"]) * np.sqrt(3.0)
    else:
        tau, Ma = None, None

    g = build_geometry(N, h, nw, flat_height=fh)
    uid = lid_profile_regularized(g["xc"], eps) if lid == "regularized" else np.ones(N)
    diag = compute_diagnostics_fields(
        d["ux"], d["uy"], {
            "N": N, "dx": g["dx"], "fluid": g["fluid"], "solid": g["solid"],
            "xc": g["xc"], "yc": g["yc"], "X": g["X"], "Y": g["Y"],
            "yb": g["yb"]}, uid,
        case_id=cid, Re=Re, N=N, h_over_H=h, Nw=nw, lid=lid, eps_lid=eps,
        U_ref=U_ref, tau=tau, Ma=Ma, rho=d["rho"] if "rho" in d else None,
        mass_drift=float(d["mass_drift"]),
        last_residual=float(d["last_residual"]))
    diag["iters"] = int(d["iters"])
    diag["flat_height"] = fh
    diag["status"] = status
    diag["group"] = meta.get("group") or (
        "fem" if cid.upper().startswith("FEM") else "main")
    for key in ("tol", "max_iter", "target_tau"):
        if key in d:
            diag[key] = d[key].item() if hasattr(d[key], "item") else d[key]
    # self-validation (same criteria as run_case)
    try:
        from checks import run_checks

        class _G:
            pass
        s = _G()
        s.geom = {"N": N, "h_over_H": h, "Nw": nw, "dx": g["dx"],
                  "fluid": g["fluid"], "flat_height": fh}
        diag["checks"] = run_checks(s, diag)
        diag["checks_ok"] = bool(diag["checks"]["checks_ok"])
    except Exception:                                   # noqa: BLE001
        diag["checks_ok"] = None
    # Preserve any D*_fneq columns (present in the run meta ;
    # not recomputable here since f is not in the main npz).
    for k in ("D_star_fneq", "D_star_fneq_dev", "D_star_fneq_ratio",
              "D_star_fneq_error"):
        if k in meta:
            diag[k] = meta[k]
    return diag


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="results")
    args = ap.parse_args(argv)

    npz_files = sorted(glob.glob(os.path.join(args.outdir, "*", "*_N*.npz")))
    table = []
    n_skip = 0
    for p in npz_files:
        try:
            diag = analyze_npz(p)
        except Exception as exc:                       # noqa: BLE001
            print(f"[skip] {p} : {exc}")
            continue
        if diag.get("status") != "converged":
            n_skip += 1
            print(f"[invalid] {p} : status '{diag.get('status')}' -> excluded from table")
            continue
        group = diag.get("group") or (
            "fem" if os.path.basename(p).startswith("fem") else "LBM")
        row = {k: diag.get(k) for k in (
            "case", "status", "Re", "N", "h_over_H", "Nw", "lid", "flat_height",
            "target_tau", "U_ref", "tau", "Ma", "Ma_max_local", "max_drho", "Ma_ok",
            "tol", "max_iter", "iters", "fluid_area", "mass_drift",
            "last_residual", "K_star", "P_star", "D_star", "Z", "balance",
            "n_vortices", "main_vortex", "checks_ok",
            "D_star_fneq", "D_star_fneq_dev", "D_star_fneq_ratio",
            "D_star_fneq_error")}
        row["group"] = group
        table.append(row)
        print(f"[ok] {row['case']} (N={row['N']}) K*={row['K_star']:.6f} "
              f"bal={100 * row['balance']:+.1f}% vort={row['n_vortices']}")

    out = os.path.join(args.outdir, "scalar_table.json")
    _locked_write(out, table, timeout=180.0, poll=0.25)
    print(f"table rebuilt : {len(table)} cases ({n_skip} excluded), "
          f"-> {out}")
    return 0


def _locked_write(path, table, timeout=180.0, poll=0.25):
    """Atomic write of the table protected by a lock (coexists with
    run_case.save_table_row if a campaign is still running)."""
    lockpath = path + ".lock"
    t0 = time.time()
    while True:
        try:
            fd = os.open(lockpath, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
            os.close(fd)
            break
        except FileExistsError:
            if time.time() - t0 > timeout:
                try:
                    os.remove(lockpath)
                except OSError:
                    pass
                continue
            time.sleep(poll)
    try:
        tmp = path + f".tmp{os.getpid()}"
        with open(tmp, "w", encoding="utf-8") as fh:
            json.dump(table, fh, indent=1, default=str)
        os.replace(tmp, path)
    finally:
        try:
            os.remove(lockpath)
        except OSError:
            pass


if __name__ == "__main__":
    main()