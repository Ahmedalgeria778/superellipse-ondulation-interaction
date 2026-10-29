# -*- coding: utf-8 -*-
"""
run_case.py
===========
Executor for a single case. Builds the WavyMRT solver, drives it to
convergence, computes the diagnostics (topology, energy, power, dissipation) and
writes to results/ :

  results/<id>/<id>_N<N>.npz   : fields u, v, rho, psi, omega + meta
  results/scalar_table.json    : table keyed by case id (CONVERGED ONLY)

Integrity rules (plan 6.4 / requirement) :
  - each record carries the identifier, (Re, h/H, Nw, N), the
    relaxation parameters (tau, omega, U_ref, nu, Ma, target_tau), the
    criterion (tol, max_iter) and the final STATUS ;
  - status 'converged'   -> results/<id>/ + row in scalar_table.json ;
  - status 'max_iter'    -> results/invalid/<id>_N<N>.npz (NO row) ;
  - status 'interrupted' / 'error' -> no save, message in the log.

Usage :
  python run_case.py G06-R400 --N 256 --outdir results
  python run_case.py ALL-MAIN --N 256
  python run_case.py V03-R100 --N 128
"""

from __future__ import annotations

import argparse
import datetime
import json
import os
import sys
import time
import traceback

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from cases import main_cases, flatref_cases, eq_cases, l1_cases
from cases import validation_cases, GROUPS
from diagnostics import compute_diagnostics, streamfunction
from checks import run_checks
from lbm_wavy import WavyMRT

# Default resolution and iteration budget
N_DEFAULT = 256
MAX_ITER = int(4e5)
TOL = 1e-8


def run_case(case, N=N_DEFAULT, outdir="results", max_iter=MAX_ITER, tol=TOL,
             save=True, log=None, case_id_override=None):
    """Run a case (serializable dict) and return its diagnostics.

    The diag['status'] field is 'converged', 'max_iter', 'interrupted' or
    'error'. Only 'converged' produces a row in scalar_table.json.
    """
    cid = case_id_override or case["id"]
    t0 = time.time()
    if log:
        log.write(f"== {cid} (N={N}) Re={case['Re']} h/H={case['h_over_H']} "
                  f"Nw={case['Nw']} lid={case['lid']}\n")
        log.flush()

    solver = None
    status = None
    try:
        solver = WavyMRT(
            case["Re"], N, case["h_over_H"], case["Nw"],
            lid=case["lid"], eps_lid=case["eps_lid"],
            flat_height=case["flat_height"],
            max_iter=max_iter, tol=tol, verbose=True, log=log)
        solver.case_id = cid
        solver.run()
    except KeyboardInterrupt:
        status = "interrupted"
        if log:
            log.write(f"   INTERRUPTED {cid} (N={N}) at {time.time()-t0:.0f}s\n")
            log.flush()
    except Exception:
        status = "error"
        if log:
            log.write(f"   ERROR {cid} (N={N}) at {time.time()-t0:.0f}s :\n")
            traceback.print_exc(file=log)
            log.flush()
    else:
        status = solver.status

    diag = {}
    if solver is not None and status in ("converged", "max_iter"):
        diag = compute_diagnostics(solver)
        diag["status"] = status
    else:
        diag["status"] = status
        diag["case"] = cid
        diag["N"] = N
        if log:
            log.write(f"   {status.upper()} -> not saved to the table.\n")
            log.flush()
        return diag

    diag["case"] = cid
    diag["group"] = case["group"]
    diag["N"] = N
    diag["flat_height"] = case.get("flat_height")
    diag["elapsed_s"] = time.time() - t0
    diag["iters"] = solver.iters
    diag["tol"] = tol
    diag["max_iter"] = max_iter
    diag["target_tau"] = solver.target_tau
    diag["omega"] = solver.omega
    diag["nu"] = solver.nu
    diag["finished_iso"] = datetime.datetime.now().isoformat(timespec="seconds")

    # --- case self-validation (requirement: every case self-checks) -----
    try:
        diag["checks"] = run_checks(solver, diag)
        diag["checks_ok"] = bool(diag["checks"]["checks_ok"])
    except Exception as exc:                             # noqa: BLE001
        diag["checks_ok"] = False
        if log:
            log.write(f"   ERROR checks {cid} : {exc}\n")
            log.flush()

    if diag.get("Ma_ok") is False:
        msg = (f"   MACH WARNING {cid} (N={N}) : Ma_max_local="
               f"{diag.get('Ma_max_local'):.4f} > threshold {diag.get('Ma_threshold')} "
               f"(max|rho-1|={diag.get('max_drho')})\n")
        if log:
            log.write(msg)
            log.flush()

    if save:
        if status == "converged":
            save_results(cid, N, solver, diag, outdir)
            save_table_row(diag, outdir)
        else:
            save_results(cid, N, solver, diag, outdir, invalid=True)
    return diag


# ----------------------------------------------------------------------
def save_results(cid, N, solver, diag, outdir="results", invalid=False):
    """Write the .npz field file plus the case .json diagnostics."""
    od = os.path.join(outdir, "invalid" if invalid else cid)
    os.makedirs(od, exist_ok=True)
    geom = solver.geom
    psi = streamfunction(solver.ux, geom)
    np.savez_compressed(
        os.path.join(od, f"{cid}_N{N}.npz"),
        ux=solver.ux, uy=solver.uy, rho=solver.rho, psi=psi,
        xc=geom["xc"], yc=geom["yc"], yb=geom["yb"],
        fluid=geom["fluid"], solid=geom["solid"],
        N=N, Re=solver.Re, case_id=cid,
        h_over_H=solver.h_over_H, Nw=solver.Nw, lid=solver.lid,
        flat_height=0.0 if solver.flat_height is None else solver.flat_height,
        eps_lid=solver.eps_lid,
        nu=solver.nu, tau=solver.tau, target_tau=solver.target_tau,
        U_ref=solver.U_ref,
        box=geom["L"], mass_drift=solver.mass_drift,
        last_residual=solver.last_residual, iters=solver.iters,
        status=diag["status"], tol=diag["tol"], max_iter=diag["max_iter"],
        Ma=solver.Ma, Ma_max_local=diag.get("Ma_max_local"),
        max_drho=diag.get("max_drho"), Ma_ok=diag.get("Ma_ok"),
        conv_hist=np.asarray(solver.conv_hist, dtype=float)
        if len(solver.conv_hist) else np.empty((0, 3)))

    # meta json (scalar quantities + vortex positions)
    iname = f"{cid}_N{N}.json" if not invalid else f"{cid}_N{N}_invalid.json"
    meta = jsonable(diag)
    with open(os.path.join(od, iname), "w", encoding="utf-8") as fh:
        json.dump(meta, fh, indent=1)
    # OPTIONAL save of the final distributions (for later D*_fneq).
    # Enabled by LBM_SAVE_F=1 ; separate file so as not to weigh down the main npz.
    if os.environ.get("LBM_SAVE_F", "0") == "1":
        ffin = getattr(solver, "f_final", None)
        if ffin is not None:
            # name WITHOUT "_N" so it is not caught by the reanalyze glob
            np.savez_compressed(os.path.join(od, f"{cid}_{N}_fdist.npz"), f=ffin)
    return od


def jsonable(diag):
    """Copy of diag without the numpy arrays (JSON-safe)."""
    out = {}
    for k, v in diag.items():
        if isinstance(v, (np.integer,)):
            out[k] = int(v)
        elif isinstance(v, (np.floating,)):
            out[k] = float(v)
        elif isinstance(v, (np.bool_,)):
            out[k] = bool(v)
        elif isinstance(v, dict):
            out[k] = {kk: (float(vv) if isinstance(vv, (np.floating,)) else
                           int(vv) if isinstance(vv, (np.integer,)) else vv)
                      for kk, vv in v.items()}
        elif isinstance(v, list):
            out[k] = [dict((kk, (float(vv) if isinstance(vv, (np.floating,))
                                else int(vv) if isinstance(vv, (np.integer,))
                                else vv))
                           for kk, vv in e.items()) if isinstance(e, dict) else e
                      for e in v]
        elif isinstance(v, tuple):
            out[k] = list(v)
        else:
            out[k] = v
    return out


TABLE_KEYS = ["case", "group", "Re", "N", "h_over_H", "Nw", "lid", "flat_height",
              "status", "target_tau", "tau", "omega", "nu", "U_ref", "Ma",
              "Ma_max_local", "max_drho", "Ma_ok",
              "tol", "max_iter", "iters", "elapsed_s", "last_residual",
              "fluid_area", "mass_drift",
              "K_star", "P_star", "D_star", "Z", "balance",
              "n_vortices", "main_vortex", "finished_iso", "checks_ok",
              "D_star_fneq", "D_star_fneq_dev", "D_star_fneq_ratio",
              "D_star_fneq_error"]


def save_table_row(diag, outdir="results", timeout=120.0, poll=0.2):
    """Add the diag row to results/scalar_table.json (merge by
    (case, N) pair : each run is identified by its resolution).

    The update is protected by a lock (lockfile + atomic write)
    so that several campaigns launched in parallel do not corrupt
    the table for each other (concurrent read-modify-write)."""
    path = os.path.join(outdir, "scalar_table.json")
    lockpath = path + ".lock"
    row = {k: diag.get(k) for k in TABLE_KEYS}
    key = f"{row['case']}#N{row['N']}"

    t0 = time.time()
    while True:
        try:
            fd = os.open(lockpath, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
            os.close(fd)
            break
        except FileExistsError:
            if time.time() - t0 > timeout:      # orphan lock : break it
                try:
                    os.remove(lockpath)
                except OSError:
                    pass
                continue
            time.sleep(poll)
    try:
        table = []
        if os.path.exists(path):
            with open(path, "r", encoding="utf-8") as fh:
                table = json.load(fh)
        table = [r for r in table if f"{r.get('case')}#N{r.get('N')}" != key]
        table.append(row)
        tmp = path + f".tmp{os.getpid()}"
        with open(tmp, "w", encoding="utf-8") as fh:
            json.dump(table, fh, indent=1, default=str)
        os.replace(tmp, path)
    finally:
        try:
            os.remove(lockpath)
        except OSError:
            pass


# ----------------------------------------------------------------------
def _dispatch(spec):
    sid = spec.upper()
    if sid.startswith("ALL-MAIN"):
        return main_cases()
    if sid.startswith("FLATREF"):
        return flatref_cases()
    if sid.startswith("EQ"):
        return eq_cases()
    if sid.startswith("L1"):
        return l1_cases()
    if sid.startswith("VALIDATION") or sid.startswith("V"):
        return validation_cases()
    if sid in GROUPS:
        return GROUPS[sid]()
    # direct identifier of an already known case
    for c in (main_cases() + flatref_cases() + eq_cases() + l1_cases()
              + validation_cases()):
        if c["id"].upper() == sid or c["id"].upper().rstrip("-R") == sid:
            return [c]
    return None


def main(argv=None):
    ap = argparse.ArgumentParser(description="Run a single corrugated-wall LBM case")
    ap.add_argument("spec", help="case id (e.g. G06-R400, V03-R100, ALL-MAIN)")
    ap.add_argument("--N", type=int, default=N_DEFAULT)
    ap.add_argument("--outdir", default="results")
    ap.add_argument("--max_iter", type=int, default=MAX_ITER)
    ap.add_argument("--tol", type=float, default=TOL)
    ap.add_argument("--skip-done", action="store_true",
                    help="do not recompute a case that has already converged (results/<id>/...npz)")
    args = ap.parse_args(argv)

    cases = _dispatch(args.spec)
    if not cases:
        ap.error(f"unknown spec : {args.spec}")
    os.makedirs(args.outdir, exist_ok=True)
    logf = open(os.path.join(args.outdir, "run.log"), "a", encoding="utf-8")
    n_ok = 0
    n_skip = 0
    for c in cases:
        if args.skip_done and _already_converged(c["id"], args.N, args.outdir):
            logf.write(f"SKIP {c['id']} (N={args.N}) already converged\n")
            logf.flush()
            n_skip += 1
            continue
        d = run_case(c, N=args.N, outdir=args.outdir,
                     max_iter=args.max_iter, tol=args.tol, log=logf)
        st = d.get("status")
        if st == "converged":
            n_ok += 1
            logf.write(f"DONE {c['id']} K*={d['K_star']:.5f} P*={d['P_star']:.4f} "
                       f"D*={d['D_star']:.4f} vort={d['n_vortices']} "
                       f"checks_ok={d.get('checks_ok')}\n")
        else:
            logf.write(f"{st.upper()} {c['id']} (N={args.N}) -> NOT SAVED\n")
        logf.flush()
    logf.close()
    print(f"done. {n_ok} converged, {n_skip} skipped. table : "
          f"{os.path.join(args.outdir, 'scalar_table.json')}")


def _already_converged(cid, N, outdir):
    """True if a converged field already exists for (cid, N)."""
    npz = os.path.join(outdir, cid, f"{cid}_N{N}.npz")
    if not os.path.exists(npz):
        return False
    try:
        d = np.load(npz)
        return str(d["status"]) == "converged" if "status" in d.files else True
    except Exception:                                   # noqa: BLE001
        return False


if __name__ == "__main__":
    main()