# -*- coding: utf-8 -*-
"""Four-band dA/dB partition used in the revised manuscript.

Same method as extract_dAdB_multires_v2.py, extended to the N=768 fields now
available at Re=400 and Re=1000.  The N=384 rows are checked against
results/partition_bands_N384.json before any new row is accepted.
"""
from __future__ import annotations

import json
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
for _p in (_ROOT, os.path.join(_ROOT, "solver")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import numpy as np
from solver.combined_geometry import build_combined_geometry
from solver.referee_mech_analysis import diss_density

FIELDS = os.path.join(_ROOT, "results", "fields")
REF = os.path.join(_ROOT, "results", "partition_bands_N384.json")
OUT = os.path.join(_HERE, "revision_dAdB.json")
A, B = 0.5, 0.5
LID_FRAC = 0.10
BOTTOM_TH = 0.10
SIDE_TH = 0.10
RES = [128, 256, 384, 512, 768]
REYS = [100, 400, 1000]


def field_path(n, alpha, Re, N):
    for suf in ("_tol1e-08", ""):
        fn = os.path.join(FIELDS, "se_n%d_w%.2f_Nw1_Re%d_N%d%s.npz"
                          % (n, alpha, Re, N, suf))
        if os.path.exists(fn):
            return fn
    return None


def phi_and_geom(n, alpha, Re, N):
    p = field_path(n, alpha, Re, N)
    if p is None:
        return None, None
    d = np.load(p)
    g = build_combined_geometry(N, n, alpha, 1, a=A, b=B)
    phi, fluid = diss_density(d, g)
    return (phi, fluid, g), p


def classify_common(g, f2, f4):
    X, Y, yb = g["X"], g["Y"], g["yb"]
    common = f2 & f4
    lid = (Y >= 1.0 - LID_FRAC) & common
    bottom = ((Y - yb) < BOTTOM_TH) & common & ~lid
    side = ((X < SIDE_TH) | (X > 1.0 - SIDE_TH)) & common & ~lid & ~bottom
    core = common & ~lid & ~bottom & ~side
    return common, {"lid": lid, "bottom": bottom, "side": side, "core": core}


def partition(N, Re, alpha=0.0):
    a2, a4 = phi_and_geom(2, alpha, Re, N), phi_and_geom(4, alpha, Re, N)
    if a2[0] is None or a4[0] is None:
        return None
    (phi2, f2, g2), _ = a2
    (phi4, f4, _), _ = a4
    dx = g2["dx"]
    D2 = float(np.sum(phi2[f2]) * dx * dx)
    common, bands = classify_common(g2, f2, f4)
    added = f4 & (~f2)
    dphi = phi4 - phi2
    band_frac = {name: float(np.sum(dphi[m]) * dx * dx) / D2
                 for name, m in bands.items()}
    dA = float(np.sum(dphi[common]) * dx * dx) / D2
    dB = float(np.sum(phi4[added]) * dx * dx) / D2
    return {"N": N, "Re": Re, "band_frac": band_frac, "dA_frac": dA,
            "dB_frac": dB, "Delta_frac": dA + dB,
            "C": sum(abs(v) for v in band_frac.values()) / abs(dA) if dA else float("nan"),
            "ID_pct": 100.0 * (dA + dB)}


def main():
    ref = json.load(open(REF, encoding="utf-8"))
    ok = True
    rows = []
    print("=== validation N=384 against partition_bands_N384.json ===")
    for Re in REYS:
        r = partition(384, Re)
        rr = ref[str(Re)]
        maxdiff = max(abs(r["band_frac"][b] - rr["bands"][b]) for b in rr["bands"])
        print("Re=%4d maxAbsDiff=%.2e" % (Re, maxdiff))
        ok &= maxdiff <= 2e-6
    if not ok:
        raise SystemExit("validation failed; no new row accepted")

    print("=== partition, all available resolutions ===")
    for Re in REYS:
        for N in RES:
            if N == 512 and Re == 100:
                continue
            r = partition(N, Re)
            if r is None:
                print("N=%d Re=%d : MISSING field" % (N, Re))
                continue
            rows.append(r)
            b = r["band_frac"]
            print("N=%3d Re=%4d | lid %+.2f bottom %+.2f side %+.2f core %+.2f"
                  " | dA %+.2f dB %+.2f D %+.2f C %.1f I_D=%.4f%%"
                  % (N, Re, b["lid"]*1e4, b["bottom"]*1e4, b["side"]*1e4,
                     b["core"]*1e4, r["dA_frac"]*1e4, r["dB_frac"]*1e4,
                     r["Delta_frac"]*1e4, r["C"], r["ID_pct"]))
    with open(OUT, "w", encoding="utf-8") as fh:
        json.dump(rows, fh, indent=2)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
