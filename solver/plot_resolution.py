# -*- coding: utf-8 -*-
"""
solver/plot_resolution.py
=========================
Figures de stabilitte de l'interaction en resolution (N=128 vs N=256).

Sorties (figures/analysis/) :
  - resolution_effect_vs_n.png : dQ(n) aux deux grilles (Q = D, P)
  - resolution_IQ_stability.png: I_D et I_P vs Re + ecart relatif
"""

from __future__ import annotations

import csv
import json
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
for _p in (_HERE, _ROOT):
    if _p not in sys.path:
        sys.path.insert(0, _p)

SCALARS_DIR = os.path.join(_ROOT, "results", "scalars")
FIGDIR = os.path.join(_ROOT, "figures", "analysis")
os.makedirs(FIGDIR, exist_ok=True)

RES = (100, 400, 1000)
NS = (2, 4, 6, 8)
ALPHA = 0.10
PAIR = (2, 4)
COL = {100: "#1f77b4", 400: "#ff7f0e", 1000: "#d62728"}


def load(Re, N):
    """Tableau {n, alpha} -> ligne du CSV pilote."""
    suf = "" if N == 128 else "_tol1e-08"
    path = os.path.join(SCALARS_DIR, "pilot_Re%d_N%d%s.csv" % (Re, N, suf))
    rows = {}
    with open(path, encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            rows[(float(r["n"]), float(r["alpha"]))] = r
    return rows


def effects(rows, key, alpha=ALPHA):
    """dQ(n)/Q(n,0) en % pour chaque n."""
    out = {}
    for n in NS:
        if (n, 0.0) in rows and (n, alpha) in rows:
            q0 = float(rows[(n, 0.0)][key])
            dq = float(rows[(n, alpha)][key]) - q0
            out[n] = 100.0 * dq / q0
    return out


def contrast(rows, key, na=PAIR[0], nb=PAIR[1], alpha=ALPHA):
    """Contraste difference de differences I_Q."""
    d_a = (float(rows[(na, alpha)][key]) - float(rows[(na, 0.0)][key]))
    d_b = float(rows[(nb, alpha)][key]) - float(rows[(nb, 0.0)][key])
    return d_b - d_a


def fig_effects():
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))
    specs = (("D_star", "D*"), ("P_star", "P*"))
    for ax, (key, lab) in zip(axes, specs):
        for Re in RES:
            for N, ls, lw in ((128, "--", 1.6), (256, "-", 2.0)):
                eff = effects(load(Re, N), key)
                xs = sorted(eff)
                ax.plot(xs, [eff[n] for n in xs], ls, color=COL[Re], lw=lw,
                        alpha=0.9 if N == 256 else 0.55,
                        marker="o" if N == 256 else "s", ms=4,
                        label="Re=%d, N=%d" % (Re, N))
        ax.set_xlabel("n (forme globale)")
        ax.set_ylabel("effet relatif de l'ondulation  d%s/%s (%%)"
                      % (lab, lab))
        ax.axhline(0, color="gray", lw=0.8)
        ax.grid(True, alpha=0.25)
        ax.legend(fontsize=7, ncol=2)
    fig.suptitle("Stabilite du contraste d'effet en resolution "
                 "(trait plein N=256, tirets N=128)")
    fig.tight_layout()
    p = os.path.join(FIGDIR, "resolution_effect_vs_n.png")
    fig.savefig(p, dpi=160)
    plt.close(fig)
    return p


def fig_stability():
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.3))
    I = {Re: {"D_star": [], "P_star": []} for Re in RES}
    for Re in RES:
        for N in (128, 256):
            rows = load(Re, N)
            I[Re]["D_star"].append(contrast(rows, "D_star"))
            I[Re]["P_star"].append(contrast(rows, "P_star"))

    ax = axes[0]
    for key, mk in (("D_star", "o"), ("P_star", "s")):
        for i, N in enumerate((128, 256)):
            ys = [abs(I[Re][key][i]) for Re in RES]
            ax.plot(RES, ys, mk, color="#1f77b4",
                    mfc="none" if N == 128 else "#1f77b4",
                    mec="#1f77b4", ms=7, lw=1.6,
                    ls="--" if N == 128 else "-",
                    label="%s, N=%d" % ("$I_D$" if key == "D_star"
                                        else "$I_P$", N))
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Re")
    ax.set_ylabel("$|I_Q|$ (contraste difference de differences)")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8)

    ax = axes[1]
    rels = {}
    for key, mk in (("D_star", "o"), ("P_star", "s")):
        rel = [100.0 * (I[Re][key][1] - I[Re][key][0]) / abs(I[Re][key][0])
               for Re in RES]
        rels[key] = rel
        ax.plot(RES, rel, mk + "-", ms=7, lw=2,
                label="$I_D$" if key == "D_star" else "$I_P$")
    ax.axhline(0, color="gray", lw=0.8)
    ax.axhline(10, color="green", lw=1.2, ls=":", label="seuil 10 %")
    ax.axhline(30, color="red", lw=1.2, ls=":", label="seuil 30 %")
    ax.set_xlabel("Re")
    ax.set_ylabel("changement relatif $N=128\\to256$ (%)")
    ax.set_ylim(-6, 34)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)

    ax = axes[2]
    ratio = [abs(I[Re]["P_star"][1]) / abs(I[Re]["D_star"][1]) for Re in RES]
    ax.plot(RES, ratio, "d-", color="#2ca02c", ms=7, lw=2)
    ax.axhline(1.0, color="gray", lw=0.8)
    ax.set_xlabel("Re")
    ax.set_ylabel("$I_P/I_D$  (coherence puissance / dissipation)")
    ax.set_ylim(0.9, 1.3)
    ax.grid(True, alpha=0.3)

    fig.suptitle("Stabilite et coherence de l'interaction "
                 "(paire 2-4, $\\alpha=0.10$)")
    fig.tight_layout()
    p = os.path.join(FIGDIR, "resolution_IQ_stability.png")
    fig.savefig(p, dpi=160)
    plt.close(fig)

    summary = {"I": I, "rel_pct": rels, "IP_over_ID": dict(zip(RES, ratio))}
    with open(os.path.join(FIGDIR, "resolution_IQ_stability.json"),
              "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2, ensure_ascii=False)
    return p


def main():
    p1 = fig_effects()
    p2 = fig_stability()
    print("[ok] %s" % p1)
    print("[ok] %s" % p2)
    for Re in RES:
        rows256 = load(Re, 256)
        rows128 = load(Re, 128)
        i256 = contrast(rows256, "D_star")
        i128 = contrast(rows128, "D_star")
        rel = 100.0 * (i256 - i128) / abs(i128)
        print("  Re=%-5d I_D(128)=%+.6f  I_D(256)=%+.6f  ecart=%+.2f %%"
              % (Re, i128, i256, rel))
    return 0


if __name__ == "__main__":
    sys.exit(main())
