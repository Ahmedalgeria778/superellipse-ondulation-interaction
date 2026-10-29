# -*- coding: utf-8 -*-
"""
06_physics/analyse_coeur_paroi.py
=================================
Comparaison HOMOGENE N=128 vs N=256 avec UN SEUL estimateur pour `n*`, et
decomposition coeur / couche de paroi de D/A.

Contexte (demande reviewer) :
  * l'ancien `n*` a N=128 et a N=256 etait extrait avec deux estimateurs
    differents (fenetre +-2 a N=128, estimateur local a 3 points a N=256) ;
    l'ecart -0.24 a Re=1000 comparait donc deux estimateurs, pas deux
    resolutions.
  * Ici on impose l'estimateur de la Table IV : ajustement parabolique sur
    les points a |n - n_argmin| <= 0.5 (5 points a pas 0.25, 3 points a pas
    0.5), aux DEUX resolutions.
  * Le coeur est la bande d'epaisseur PHYSIQUE fixe d(x,Gamma) >= d0, d0=0.05
    (seule variante convergencee en maillage : mediane 0.52 %).

Donnees : N=128 cherche dans data_all190, data_refine, data_refine_seed ;
          N=256 dans data_conv256.
Sorties : 06_physics/results_coeur_paroi/{coeur_paroi.md,csv,json}
          figs/PH13_core_wall_decomposition.png (via FIG_DIR)
"""
from __future__ import annotations

import csv
import json
import os
import sys
from pathlib import Path

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

_HERE = Path(__file__).resolve().parent
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))
import analyse_etape2_3_coeur as base            # noqa: E402
from etude_convergence_paroi import gradient_valid_mask  # noqa: E402

ROOT = _HERE.parent
DIR128 = (ROOT / "05_dataset" / "data_all190",
          ROOT / "05_dataset" / "data_refine",
          ROOT / "05_dataset" / "data_refine_seed")
DIR256 = (ROOT / "05_dataset" / "data_conv256",)
OUT = _HERE / "results_coeur_paroi"
FIG_DIR = ROOT / "08_operator_results" / "PH13_manuscript_190" / "figs"

GRID128 = (1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75,
           4.0, 4.25, 4.5, 5.0, 6.0, 8.0)
GRID256 = (1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75,
           4.0, 4.25, 4.5, 5.0, 6.0, 8.0)
RE_MAIN = (300.0, 500.0, 1000.0)
RE_LOW = (300.0, 325.0, 350.0)
RE_HIGH = (850.0, 875.0, 900.0)
D0 = 0.05


def case_path(n, Re, N):
    dirs = DIR128 if N == 128 else DIR256
    for d in dirs:
        p = d / f"{base.case_id(n, Re, N)}.npz"
        if p.exists() and "_ckpt" not in p.name:
            return p
    return None


def observables(n, Re, N):
    """(D/A)_total sur A_valid, (D/A)_core d0, part de couche de paroi, meta."""
    p = case_path(n, Re, N)
    if p is None:
        return None
    with np.load(p) as z:
        ux = np.asarray(z["ux"], float)
        uy = np.asarray(z["uy"], float)
        fluid = np.asarray(z["fluid"]).astype(bool)
        xx = np.asarray(z["X"], float)
        yy = np.asarray(z["Y"], float)
    phi, valid, dist, dx = base.phi_and_dist(ux, uy, fluid, xx, Re)
    cell = dx * dx
    D_tot = 2.0 * float(np.sum(phi[valid])) * cell
    A_valid = float(np.count_nonzero(valid)) * cell
    kc = valid & (dist >= D0)
    D_core = 2.0 * float(np.sum(phi[kc])) * cell
    A_core = float(np.count_nonzero(kc)) * cell
    return {"n": n, "Re": Re, "N": N,
            "D_total": D_tot, "A_valid": A_valid, "DA_total": D_tot / A_valid,
            "D_core": D_core, "A_core": A_core, "DA_core": D_core / A_core,
            "wall_share_pct": 100.0 * (D_tot - D_core) / D_tot,
            "path": str(p)}


def series(ns, Re, N, key):
    pts = []
    for n in ns:
        o = observables(n, Re, N)
        if o is not None:
            pts.append((n, o[key]))
    return pts


def vertex_window(pts):
    """Estimateur Table IV : argmin puis polyfit deg 2 sur |n-nmin| <= 0.5."""
    if len(pts) < 3:
        return None
    nmin = min(pts, key=lambda t: t[1])[0]
    win = [(n, v) for n, v in pts if abs(n - nmin) <= 0.5 + 1e-9]
    if len(win) < 3:
        return None
    c = np.polyfit([p[0] for p in win], [p[1] for p in win], 2)
    if c[0] <= 0:
        return {"argmin": nmin, "vertex": None, "window": [p[0] for p in win],
                "margin_pct": None, "n_win": len(win)}
    vtx = float(-c[1] / (2 * c[0]))
    fstar = float(np.polyval(c, vtx))
    best = min(v for n, v in win if abs(n - nmin) > 1e-9)
    return {"argmin": nmin, "vertex": vtx, "window": [p[0] for p in win],
            "margin_pct": 100.0 * (best - fstar) / fstar, "n_win": len(win)}


def vertex_local3(pts):
    i = min(range(len(pts)), key=lambda j: pts[j][1])
    if i == 0 or i == len(pts) - 1:
        return None
    w = pts[i - 1:i + 2]
    c = np.polyfit([p[0] for p in w], [p[1] for p in w], 2)
    if c[0] <= 0:
        return None
    return float(-c[1] / (2 * c[0]))


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    for N, grid in ((128, GRID128), (256, GRID256)):
        for Re in (300.0, 325.0, 350.0, 500.0, 850.0, 875.0, 900.0, 1000.0):
            for n in grid:
                o = observables(n, Re, N)
                if o is not None:
                    rows.append(o)

    with (OUT / "coeur_paroi.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    L = ["# Coeur (d0=0.05) vs couche de paroi -- estimateur Table IV (+-0.5)",
         "",
         "- N=128 : data_all190 + data_refine + data_refine_seed.",
         "- N=256 : data_conv256.",
         "- Estimateur unique : argmin puis parabole sur |n-n_argmin| <= 0.5.",
         ""]

    summary = {"window_vertex": {}, "local3": {}}

    def report(tag, key, res):
        t1 = []
        for Re in res:
            a = vertex_window(series(GRID128, Re, 128, key))
            b = vertex_window(series(GRID256, Re, 256, key))
            va = a["vertex"] if a else None
            vb = b["vertex"] if b else None
            d = (vb - va) if (va is not None and vb is not None) else None
            summary["window_vertex"].setdefault(tag, {})[f"Re{Re:g}"] = {
                "N128": a, "N256": b, "delta": d}
            t1.append((Re, a, b, d))
        L.append(f"## {tag} -- estimateur +-0.5 (meme methode aux 2 N)")
        L.append("")
        L.append("| Re | n* N=128 | n* N=256 | Delta (256-128) | win128 | win256 |")
        L.append("|---|---|---|---|---|---|")
        dmax = 0.0
        for Re, a, b, d in t1:
            wa = ",".join(f"{x:g}" for x in a["window"]) if a else "-"
            wb = ",".join(f"{x:g}" for x in b["window"]) if b else "-"
            av = f"{a['vertex']:.3f}" if a and a["vertex"] is not None else "-"
            bv = f"{b['vertex']:.3f}" if b and b["vertex"] is not None else "-"
            ds = f"{d:+.3f}" if d is not None else "-"
            if d is not None:
                dmax = max(dmax, abs(d))
            L.append(f"| {Re:g} | {av} | {bv} | {ds} | {wa} | {wb} |")
        L.append("")
        L.append(f"Ecart maximal |Delta n*| = **{dmax:.3f}** "
                 f"(critere demande <= 0.15).")
        L.append("")

    report("D/A total (A_valid)", "DA_total", RE_MAIN)
    report("D_core/A_core (d0=0.05)", "DA_core", RE_MAIN)

    # --- 3-point local estimateur, pour continuite avec la redaction --------
    L.append("## n* par estimateur local a 3 points (variante total)")
    L.append("")
    L.append("| Re | n*(N=128) | n*(N=256) |")
    L.append("|---|---|---|")
    for Re in RE_MAIN:
        a = vertex_local3(series(GRID128, Re, 128, "DA_total"))
        b = vertex_local3(series(GRID256, Re, 256, "DA_total"))
        summary["local3"][f"Re{Re:g}"] = {"N128": a, "N256": b}
        L.append(f"| {Re:g} | " + (f"{a:.3f}" if a else "-")
                 + " | " + (f"{b:.3f}" if b else "-") + " |")
    L.append("")

    # --- coeur a N=128 sur la grille raffinee (bandes existantes) ----------
    L.append("## n* du coeur d0=0.05 sur la grille N=128 raffinee")
    L.append("")
    L.append("| Re | argmin | n* (+-0.5) | marge % | n_win |")
    L.append("|---|---|---|---|---|")
    for Re in RE_LOW + RE_HIGH:
        p = series(GRID128, Re, 128, "DA_core")
        v = vertex_window(p)
        summary["window_vertex"].setdefault("core_N128_bandes", {})[
            f"Re{Re:g}"] = v
        if v:
            vv = f"{v['vertex']:.3f}" if v["vertex"] is not None else "n/a"
            mg = f"{v['margin_pct']:.3f}" if v["margin_pct"] is not None else "n/a"
            L.append(f"| {Re:g} | {v['argmin']:g} | {vv} | "
                     f"{mg} | {v['n_win']} |")
        else:
            L.append(f"| {Re:g} | - | - | - | - |")
    L.append("")

    # --- stabilite du coeur 128 -> 256 sur points communs ------------------
    L.append("## Stabilite 128 -> 256 du coeur d0=0.05")
    L.append("")
    L.append("| re | var D_core/A |")
    L.append("|---|---|")
    for Re in RE_MAIN:
        vs = []
        for n in GRID128:
            a = observables(n, Re, 128)
            b = observables(n, Re, 256)
            if a and b:
                vs.append(b["DA_core"] / a["DA_core"] - 1.0)
        if vs:
            L.append(f"| Re={Re:g} | mediane {100*np.median(vs):+.2f} % "
                     f"(min {100*np.min(vs):+.2f} %, max {100*np.max(vs):+.2f} %, "
                     f"{len(vs)} cas) |")
    L.append("")

    (OUT / "coeur_paroi.json").write_text(
        json.dumps(summary, indent=2, default=float), encoding="utf-8")
    (OUT / "coeur_paroi.md").write_text("\n".join(L) + "\n", encoding="utf-8")

    # ---------------------------------------------------------------- figure
    _figure(rows)
    print("\n".join(L))
    print(f"[coeur_paroi] ecrit -> {OUT} et {FIG_DIR}", flush=True)
    return 0


def _figure(rows):
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.4), sharey=False)
    for ax, Re in zip(axes, (300.0, 1000.0)):
        r = [x for x in rows if x["Re"] == Re and x["N"] == 256]
        if not r:
            continue
        r.sort(key=lambda x: x["n"])
        n = [x["n"] for x in r]
        tot = [x["DA_total"] for x in r]
        core = [x["DA_core"] for x in r]
        share = [x["wall_share_pct"] for x in r]
        ax.plot(n, tot, "-o", color="#1f4e79", ms=3.5, lw=1.4,
                label=r"$D/A$ (full cavity, $A_{\rm valid}$)")
        ax.plot(n, core, "-s", color="#c0392b", ms=3.5, lw=1.4,
                label=r"$D_{\rm core}/A_{\rm core}$ ($d_0=0.05$)")
        it = int(np.argmin(tot))
        ic = int(np.argmin(core))
        ax.plot(n[it], tot[it], "*", color="#1f4e79", ms=13,
                markeredgecolor="k", markeredgewidth=0.4, zorder=5)
        ax.plot(n[ic], core[ic], "*", color="#c0392b", ms=13,
                markeredgecolor="k", markeredgewidth=0.4, zorder=5)
        ax.set_xlabel("geometry exponent $n$")
        ax.set_ylabel(r"dissipation density")
        ax.set_title(f"$Re={Re:g}$, $N=256$")
        ax.grid(alpha=0.25)
        ax2 = ax.twinx()
        ax2.plot(n, share, "--", color="#7f8c8d", lw=1.1,
                 label="wall-band share (%)")
        ax2.set_ylabel("wall-band share of $D$ (%)", color="#7f8c8d")
        ax2.tick_params(axis="y", colors="#7f8c8d")
        h1, l1 = ax.get_legend_handles_labels()
        h2, l2 = ax2.get_legend_handles_labels()
        ax.legend(h1 + h2, l1 + l2, fontsize=8, loc="upper center",
                  frameon=False)
    fig.tight_layout()
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG_DIR / "PH13_core_wall_decomposition.png", dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    raise SystemExit(main())
