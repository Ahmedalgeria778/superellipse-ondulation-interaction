# -*- coding: utf-8 -*-
"""
06_physics/analyse_etape2_3_coeur.py

Etape 2 : bande a epaisseur PHYSIQUE fixe. Au lieu d'exclure des disques en
cellules, on retire la couche d(x,Gamma) < d0 avec d0 in {0.05, 0.1, 0.2}
(longueur de couche limite ~ L/sqrt(Re) ~ 0.06 a Re=1000). On calcule
D_core(d0) et A_core(d0) — tous deux sur le MEME masque discret — sur les
30 cas (10 n x Re=300/500/1000) a N=128 et N=256.

Critere : variation 128->256 de D_core/A_core < ~3 % => quantite convergee ;
on recalcule alors n*(d0) aux deux resolutions (accord a +-0.15 requis).

Etape 3 : incertitude sur n*. Pour chaque Re et chaque variante, ecart-type de
R_D(n) = (D/A)(n)/min_n(D/A) sur la grille echantillonnee, compare a la marge
du creux (0.1-0.4 %), plus argmin a N=128 et N=256.

Point 4 : D/A total sur n in {3,4,5,6,8}, Re = 300/500/1000, aux deux N ;
rapport D/A(n=3)/D/A(n=4) a Re=300.

Sorties : 06_physics/results_coeur/
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np
from scipy import ndimage

_HERE = Path(__file__).resolve().parent
_ROOT = _HERE.parent
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))
from etude_convergence_paroi import gradient_valid_mask  # noqa: E402

ALL190 = _ROOT / "05_dataset" / "data_all190"
CONV256 = _ROOT / "05_dataset" / "data_conv256"
OUT_DIR = _HERE / "results_coeur"

NS_GRID = (1.25, 1.5, 1.75, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0)
RES = (300.0, 500.0, 1000.0)
D0 = (0.05, 0.1, 0.2)
JUNCTIONS = ((1.0, 0.0), (-1.0, 0.0))


def case_id(n, Re, N):
    return f"se_n{str(round(n, 2)).replace('.', 'p')}_Re{int(round(Re))}_N{N}"


def load(n, Re, N):
    d = ALL190 if N == 128 else CONV256
    p = d / f"{case_id(n, Re, N)}.npz"
    if not p.exists() or "_ckpt" in p.name:
        return None
    with np.load(p) as z:
        ux = np.asarray(z["ux"], float); uy = np.asarray(z["uy"], float)
        fluid = np.asarray(z["fluid"]).astype(bool)
        xx = np.asarray(z["X"], float); yy = np.asarray(z["Y"], float)
    return ux, uy, fluid, xx, yy


def phi_and_dist(ux, uy, fluid, xx, Re):
    dx = float(abs(xx[0, 1] - xx[0, 0]))
    valid = gradient_valid_mask(fluid)
    gx = np.full(ux.shape, np.nan); gy = np.full(ux.shape, np.nan)
    vx = np.full(ux.shape, np.nan); vy = np.full(ux.shape, np.nan)
    gx[1:-1, 1:-1] = (ux[1:-1, 2:] - ux[1:-1, :-2]) / (2 * dx)
    vx[1:-1, 1:-1] = (uy[1:-1, 2:] - uy[1:-1, :-2]) / (2 * dx)
    gy[1:-1, 1:-1] = (ux[2:, 1:-1] - ux[:-2, 1:-1]) / (2 * dx)
    vy[1:-1, 1:-1] = (uy[2:, 1:-1] - uy[:-2, 1:-1]) / (2 * dx)
    phi = np.where(valid, (2 * gx ** 2 + 2 * vy ** 2 + (gy + vx) ** 2) / Re, 0.0)
    dist = ndimage.distance_transform_edt(fluid) * dx
    return phi, valid, dist, dx


def variant_values(phi, valid, dist, xx, yy, cell):
    out = {}
    keep = valid
    out["total"] = (2.0 * float(np.sum(phi[keep])) * cell,
                    float(np.count_nonzero(keep)) * cell)
    keepj = valid.copy()
    for (xj, yj) in JUNCTIONS:
        keepj &= np.hypot(xx - xj, yy - yj) > 0.25
    out["r0.25"] = (2.0 * float(np.sum(phi[keepj])) * cell,
                    float(np.count_nonzero(keepj)) * cell)
    for d0 in D0:
        k = valid & (dist >= d0)
        out[f"d0={d0:g}"] = (2.0 * float(np.sum(phi[k])) * cell,
                            float(np.count_nonzero(k)) * cell)
    return out


def vertex(nv, yv):
    if len(nv) < 3:
        return float("nan")
    c = np.polyfit(np.asarray(nv, float), np.asarray(yv, float), 2)
    if c[0] == 0:
        return float("nan")
    return float(-c[1] / (2 * c[0]))


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    rows = []
    for n in NS_GRID:
        for Re in RES:
            for N in (128, 256):
                got = load(n, Re, N)
                if got is None:
                    print(f"[coeur] MANQUE n={n} Re={Re} N={N}")
                    continue
                ux, uy, fluid, xx, yy = got
                phi, valid, dist, dx = phi_and_dist(ux, uy, fluid, xx, Re)
                vv = variant_values(phi, valid, dist, xx, yy, dx * dx)
                for name, (D, A) in vv.items():
                    rows.append({"n": n, "Re": Re, "N": N, "variant": name,
                                 "D": D, "A": A,
                                 "D_over_A": D / A if A else float("nan")})
                print(f"[coeur] n={n} Re={Re:g} N={N} OK", flush=True)

    with (OUT_DIR / "etape2_core_30cas.csv").open("w", newline="",
                                                  encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)

    variants = ["total", "r0.25"] + [f"d0={d:g}" for d in D0]
    summary = {"variants": variants, "variation_DoverA_128_256": {},
               "nstar": {}, "point4_DoverA": {}}

    # --- Etape 2 : variation 128 -> 256 de D_core/A_core -----------------
    lines = ["# Etape 2 — coeur a epaisseur physique fixe (D_core/A_core)", "",
             "| variante | mediane | min | max | n cas |", "|---|---|---|---|---|"]
    for var in variants:
        vs = []
        for n in NS_GRID:
            for Re in RES:
                a = next((r for r in rows if r["n"] == n and r["Re"] == Re
                          and r["N"] == 128 and r["variant"] == var), None)
                b = next((r for r in rows if r["n"] == n and r["Re"] == Re
                          and r["N"] == 256 and r["variant"] == var), None)
                if a and b:
                    vs.append(b["D_over_A"] / a["D_over_A"] - 1.0)
        if vs:
            summary["variation_DoverA_128_256"][var] = {
                "median_pct": 100 * float(np.median(vs)),
                "min_pct": 100 * float(np.min(vs)),
                "max_pct": 100 * float(np.max(vs)), "count": len(vs)}
            lines.append(f"| {var} | {100*np.median(vs):.2f} % | "
                         f"{100*np.min(vs):.2f} % | {100*np.max(vs):.2f} % "
                         f"| {len(vs)} |")

    # --- n*(variante, N, Re) -------------------------------------------
    lines += ["", "## n* (argmin + sommet parabolique) par Re", "",
              "| variante | Re | n* N=128 (argmin/vertex) | n* N=256 "
              "(argmin/vertex) |", "|---|---|---|---|"]
    def curve(var, Re, N):
        return [(n, next((r["D_over_A"] for r in rows if r["n"] == n
                          and r["Re"] == Re and r["N"] == N
                          and r["variant"] == var), None))
                for n in NS_GRID]
    for var in variants:
        for Re in RES:
            cells = []
            for N in (128, 256):
                c = [(n, v) for n, v in curve(var, Re, N) if v is not None]
                ar = min(c, key=lambda t: t[1])[0]
                i = [j for j, (nn, _) in enumerate(c) if nn == ar][0]
                lo = max(0, i - 2); hi = min(len(c), i + 3)
                vtx = vertex([nn for nn, _ in c[lo:hi]],
                             [vv for _, vv in c[lo:hi]])
                cells.append((ar, vtx))
                summary["nstar"].setdefault(var, {}).setdefault(f"Re{Re:g}",
                                                                {})[f"N{N}"] = {
                    "argmin": ar, "vertex": vtx}
            lines.append(f"| {var} | {Re:g} | {cells[0][0]:g} / "
                         f"{cells[0][1]:.2f} | {cells[1][0]:g} / "
                         f"{cells[1][1]:.2f} |")

    # --- Etape 3 : ecart-type de R_D(n) ---------------------------------
    lines += ["", "## Etape 3 — dispersion de R_D(n) sur la grille", "",
              "| variante | Re | N | argmin | std(R_D) % | amplitude "
              "(max-min) % |", "|---|---|---|---|---|---|"]
    for var in variants:
        for Re in RES:
            for N in (128, 256):
                c = [(n, v) for n, v in curve(var, Re, N) if v is not None]
                if not c:
                    continue
                ys = np.array([v for _, v in c], float)
                rd = ys / ys.min()
                summary.setdefault("R_D", {}).setdefault(var, {}).setdefault(
                    f"Re{Re:g}", {})[f"N{N}"] = {
                        "argmin": c[int(np.argmin(ys))][0],
                        "std_pct": 100 * float(np.std(rd)),
                        "amp_pct": 100 * float(np.max(rd) - np.min(rd))}
                lines.append(f"| {var} | {Re:g} | {N} | "
                             f"{c[int(np.argmin(ys))][0]:g} | "
                             f"{100*np.std(rd):.3f} % | "
                             f"{100*(np.max(rd)-np.min(rd)):.3f} % |")

    # --- Point 4 : D/A total sur n=3..8, rapport n=3/n=4 a Re=300 -------
    lines += ["", "## Point 4 — D/A total, creux n=3..8", "",
              "| Re | N | argmin | D/A(3)/D/A(4) |", "|---|---|---|---|"]
    for Re in RES:
        for N in (128, 256):
            def g(n):
                return next((r["D_over_A"] for r in rows if r["n"] == n
                             and r["Re"] == Re and r["N"] == N
                             and r["variant"] == "total"), None)
            sub = [(n, g(n)) for n in (2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0)]
            sub = [(n, v) for n, v in sub if v is not None]
            ar = min(sub, key=lambda t: t[1])[0] if sub else None
            ratio = (g(3.0) / g(4.0)) if (g(3.0) and g(4.0)) else None
            summary["point4_DoverA"][f"Re{Re:g}_N{N}"] = {
                "argmin": ar, "ratio_3_over_4": ratio}
            lines.append(f"| {Re:g} | {N} | {ar:g} | {ratio:.6f} |")

    (OUT_DIR / "etape2_3_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8")
    (OUT_DIR / "etape2_3_resume.md").write_text("\n".join(lines),
                                                encoding="utf-8")
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
