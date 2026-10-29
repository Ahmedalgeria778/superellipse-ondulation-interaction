# -*- coding: utf-8 -*-
"""
extract_dAdB_multires_v2.py
===========================
Four-band dA/dB partition of the (n=4)-(n=2) dissipation contrast at alpha=0,
for N in {128, 256, 384, 512}, Re in {100, 400, 1000}.

REUSES the authoritative machinery from solver/referee_mech_analysis.py
(build_combined_geometry + diss_density, which use the `fluid` mask and the
one-sided wall/lid gradient correction). This guarantees the recomputed bands
are consistent with partition_bands_N384.json (the source of Table SM.XIII),
so the N=256/512 rows can be added to Table SM.XV without mixing methods.

Output: results/dAdB_multires_v2.json
Prints a validation block comparing N=384 to partition_bands_N384.json.
"""
from __future__ import annotations
import json, os, sys

_HERE = os.path.dirname(os.path.abspath(__file__))
for _p in (_HERE, os.path.join(_HERE, "solver")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import numpy as np
from solver.combined_geometry import build_combined_geometry
from solver.referee_mech_analysis import diss_density

FIELDS = os.path.join(_HERE, "results", "fields")
A, B = 0.5, 0.5
LID_FRAC = 0.10      # top 10%  -> Y >= 1 - 0.10
BOTTOM_TH = 0.10     # Y < yb + 0.10
SIDE_TH = 0.10       # X < 0.10 or X > 0.90
RES = [128, 256, 384, 512]
RES_ALL = [100, 400, 1000]


def field_path(n, alpha, Re, N):
    for suf in ("_tol1e-08", ""):
        fn = os.path.join(
            FIELDS, "se_n%d_w%.2f_Nw1_Re%d_N%d%s.npz" % (n, alpha, Re, N, suf))
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
    """Return dict of boolean masks for the four bands within the COMMON region."""
    X, Y, yb = g["X"], g["Y"], g["yb"]
    common = f2 & f4
    lid = (Y >= 1.0 - LID_FRAC) & common
    near_bot = ((Y - yb) < BOTTOM_TH) & common
    bottom = near_bot & ~lid
    side = ((X < SIDE_TH) | (X > 1.0 - SIDE_TH)) & common & ~lid & ~bottom
    core = common & ~lid & ~bottom & ~side
    return common, {"lid": lid, "bottom": bottom, "side": side, "core": core}


def partition(N, Re, alpha=0.0):
    a2 = phi_and_geom(2, alpha, Re, N)
    a4 = phi_and_geom(4, alpha, Re, N)
    if a2[0] is None or a4[0] is None:
        return None
    (phi2, f2, g2), _ = a2
    (phi4, f4, _), _ = a4
    dx = g2["dx"]
    # reference absolute dissipation of n=2 (nu-free phi units; nu cancels in ratio)
    D2 = float(np.sum(phi2[f2]) * dx * dx)
    D4 = float(np.sum(phi4[f4]) * dx * dx)
    common, bands = classify_common(g2, f2, f4)
    added = f4 & (~f2)
    dphi = phi4 - phi2
    band_frac = {}
    for name, m in bands.items():
        band_frac[name] = float(np.sum(dphi[m]) * dx * dx) / D2
    dA_frac = float(np.sum(dphi[common]) * dx * dx) / D2
    dB_frac = float(np.sum(phi4[added]) * dx * dx) / D2
    Delta_frac = dA_frac + dB_frac
    C = (sum(abs(v) for v in band_frac.values()) / abs(dA_frac)) if dA_frac else float("nan")
    return {"N": N, "Re": Re, "D2": D2, "D4": D4,
            "band_frac": band_frac, "dA_frac": dA_frac, "dB_frac": dB_frac,
            "Delta_frac": Delta_frac, "C": C,
            "ID_pct": 100.0 * Delta_frac,
            "A_added_frac": float(added.sum()) * dx * dx / float(g2["fluid_area"] + added.sum() * dx * dx)}


def main():
    out = []
    print("=== validation: recomputed N=384 vs partition_bands_N384.json ===")
    ref = json.load(open(os.path.join(_HERE, "results", "partition_bands_N384.json")))
    ok = True
    for Re in RES_ALL:
        r = partition(384, Re)
        out.append(r)
        rr = ref[str(Re)]
        maxdiff = max(abs(r["band_frac"][b] - rr["bands"][b]) for b in rr["bands"])
        print("Re=%4d  bands lid %+.6f/%+.6f bottom %+.6f/%+.6f side %+.6f/%+.6f core %+.6f/%+.6f  dA %+.6f/%+.6f dB %+.6f/%+.6f | maxAbsDiff=%.2e"
              % (Re, r["band_frac"]["lid"], rr["bands"]["lid"], r["band_frac"]["bottom"], rr["bands"]["bottom"],
                 r["band_frac"]["side"], rr["bands"]["side"], r["band_frac"]["core"], rr["bands"]["core"],
                 r["dA_frac"], rr["dA"], r["dB_frac"], rr["dB"], maxdiff))
        if maxdiff > 2e-6:
            ok = False
    print("VALIDATION", "PASS (reproducible to 2e-6)" if ok else "FAIL -> do NOT mix into SM.XV")

    if not ok:
        print("Aborting: method does not reproduce the source table.")
        return

    print("\n=== all resolutions ===")
    for Re in RES_ALL:
        for N in RES:
            if N == 512 and Re == 100:
                continue
            if N == 384:
                continue  # already have
            r = partition(N, Re)
            if r is None:
                print("N=%d Re=%d : MISSING field" % (N, Re))
                continue
            out.append(r)
            b = r["band_frac"]
            print("N=%3d Re=%4d | lid %+.2f bottom %+.2f side %+.2f core %+.2f | dA %+.2f dB %+.2f D %+.2f C %.1f  (x1e4)  I_D=%.4f%%"
                  % (N, Re, b["lid"]*1e4, b["bottom"]*1e4, b["side"]*1e4, b["core"]*1e4,
                     r["dA_frac"]*1e4, r["dB_frac"]*1e4, r["Delta_frac"]*1e4, r["C"], r["ID_pct"]))

    with open(os.path.join(_HERE, "results", "dAdB_multires_v2.json"), "w") as fh:
        json.dump(out, fh, indent=2)
    print("\nwrote results/dAdB_multires_v2.json")


if __name__ == "__main__":
    main()
