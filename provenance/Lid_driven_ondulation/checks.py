# -*- coding: utf-8 -*-
"""
checks.py
=========
Self-validation of each case of the campaign (requirement: each case must be
verified automatically before entering the table).

Criteria (all reported ; `checks_ok` = AND of the blocking criteria) :

  1. area     : discrete fluid area ~ analytic (1 - h/2, or 1 - flat_height) ;
  2. corners  : y_b(0) = y_b(1) = 0 (wavy wall, zero phase) ;
  3. mask     : every column contains fluid, lid (top row)
                entirely fluid, fluid region connected ;
  4. detection: the detector is immune to corners -- every vortex
                detected is at a distance >= margin from any wall ;
  5. main vortex : present, negative sign, significant amplitude ;
  6. balance  : |P* - D*| / D* under control (report only, known FD bias) ;
  7. convergence : status 'converged' and residual <= tol ;
  8. Mach     : Ma_max_local <= threshold (report) ;
  9. mass     : small relative mass drift.
"""

from __future__ import annotations

import numpy as np
from scipy import ndimage

TOL_AREA = 2e-2       # 2 % (mask staircase at low resolution)
TOL_BALANCE = 0.15    # 15 % (systematic wall FD bias)
TOL_MASS = 1e-3
MAIN_AMP_MIN = 0.02


def _corner_values(geom):
    """y_b extrapolated to the two bottom corners (continuous)."""
    h = geom.get("h_over_H", 0.0)
    Nw = geom.get("Nw", 1)
    flat = geom.get("flat_height")
    if flat is not None:
        return float(flat), float(flat)
    if h == 0.0:
        return 0.0, 0.0
    from geometry import wall_height
    return float(wall_height(np.array([0.0]), h, Nw)[0]), \
        float(wall_height(np.array([1.0]), h, Nw)[0])


def run_checks(solver, diag):
    """Returns a dict of checks for a case (bools + values)."""
    geom = solver.geom
    N = geom["N"]
    dx = geom["dx"]
    fluid = geom["fluid"].astype(bool)

    checks = {}

    # 1) fluid area ----------------------------------------------------
    h = float(geom.get("h_over_H", 0.0))
    flat = geom.get("flat_height")
    area_an = (1.0 - float(flat)) if flat is not None else (1.0 - 0.5 * h)
    area = float(diag.get("fluid_area", geom.get("fluid_area")))
    checks["area_analytic"] = area_an
    checks["area_meas"] = area
    checks["area_err_pct"] = 100.0 * (area / area_an - 1.0) if area_an > 0 else np.nan
    checks["area_ok"] = bool(abs(area / area_an - 1.0) <= TOL_AREA)

    # 2) corners ----------------------------------------------------------
    # y_b at the two corners must equal the analytic value of the geometry:
    # 0 for a zero-phase wavy wall (or a flat bottom at h=0); equal to
    # flat_height for a raised flat bottom (EQ case). The two corners must
    # also match (continuity of the wall).
    c0, c1 = _corner_values(geom)
    c_exp = 0.0 if flat is None else float(flat)
    checks["corner_yb_left"] = c0
    checks["corner_yb_right"] = c1
    checks["corner_expected"] = c_exp
    checks["corners_ok"] = bool(
        max(abs(c0 - c_exp), abs(c1 - c_exp), abs(c0 - c1)) < 1e-9)

    # 3) mask ---------------------------------------------------------
    col_has_fluid = bool(np.all(fluid.any(axis=0)))
    top_fluid_frac = float(fluid[N - 1, :].mean())
    lab, nlab = ndimage.label(fluid)
    largest = int(np.bincount(lab.ravel())[1:].max()) if nlab else 0
    connected = bool(largest == int(fluid.sum()))
    checks["mask_cols_ok"] = col_has_fluid
    checks["mask_top_fluid_frac"] = top_fluid_frac
    checks["mask_connected"] = connected
    checks["mask_ok"] = bool(col_has_fluid and connected and top_fluid_frac > 0.99)

    # 4) detection immune to corners -------------------------------------
    vort = diag.get("vortices") or []
    margin = 0.02
    worst = min((v.get("wall_dist", 1.0) for v in vort), default=1.0)
    checks["n_vortices"] = int(diag.get("n_vortices", len(vort)))
    checks["vortex_min_wall_dist"] = float(worst)
    checks["detection_ok"] = bool(worst >= margin - 1e-9)

    # 5) main vortex ------------------------------------------
    main = diag.get("main_vortex")
    checks["main_present"] = main is not None
    checks["main_sign"] = int(main["sign"]) if main else 0
    checks["main_amp"] = float(main["amp"]) if main else 0.0
    checks["main_ok"] = bool(main is not None and main["sign"] == -1
                             and main["amp"] >= MAIN_AMP_MIN)

    # 6) energy balance ---------------------------------------------
    bal = diag.get("balance", np.nan)
    checks["balance_pct"] = 100.0 * bal if bal == bal else np.nan
    checks["balance_ok"] = bool(bal == bal and abs(bal) <= TOL_BALANCE)

    # 7) convergence ----------------------------------------------------
    checks["status"] = diag.get("status")
    checks["last_residual"] = float(diag.get("last_residual", np.inf)
                                    if diag.get("last_residual") is not None
                                    else np.inf)
    checks["tol"] = float(diag.get("tol", 1e-8))
    checks["converged_ok"] = bool(
        diag.get("status") == "converged"
        and checks["last_residual"] <= checks["tol"])

    # 8) Mach -----------------------------------------------------------
    mma = diag.get("Ma_max_local")
    checks["Ma_max_local"] = float(mma) if mma is not None else np.nan
    checks["Ma_ok"] = bool(diag.get("Ma_ok", False))

    # 9) mass ----------------------------------------------------------
    md = diag.get("mass_drift")
    checks["mass_drift"] = float(md) if md is not None else np.nan
    checks["mass_ok"] = bool(md is not None and md <= TOL_MASS)

    blocking = ["area_ok", "corners_ok", "mask_ok", "detection_ok",
                "main_ok", "converged_ok"]
    checks["checks_ok"] = bool(all(checks[k] for k in blocking))
    return checks
