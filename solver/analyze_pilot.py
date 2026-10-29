# -*- coding: utf-8 -*-
"""
solver/analyze_pilot.py
=======================
Analyse de la campagne (resolution fixe) : figures + table de synthese.

Sorties (figures/analysis/) :
  - effect_vs_n_Re<Re>.png   : effet de l'ondulation dQ(n) pour Q=D,P
  - interaction_vs_pair.png  : I_Q pour les paires (n_a < n_b)
  - fields_<Re>.png          : champs (ux, psi) des cas n*alpha
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
_PROV_A = os.path.join(_ROOT, "provenance", "Lid_driven_ondulation")
for _p in (_HERE, _PROV_A):
    if _p not in sys.path:
        sys.path.insert(0, _p)

SCALARS_DIR = os.path.join(_ROOT, "results", "scalars")
FIELDS_DIR = os.path.join(_ROOT, "results", "fields")
FIGDIR = os.path.join(_ROOT, "figures", "analysis")
os.makedirs(FIGDIR, exist_ok=True)


def load_pilot_table(Re, N, alpha_hi):
    tab = os.path.join(SCALARS_DIR, "pilot_Re%g_N%d.csv" % (Re, N))
    import csv
    rows = {}
    with open(tab, encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            rows[(float(r["n"]), float(r["alpha"]))] = r
    return rows


def fig_effect_vs_n(rows, Re):
    ns = sorted({k[0] for k in rows})
    alphas = sorted({k[1] for k in rows if k[1] > 0})
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    for ax, key, lab in ((axes[0], "D_star", "D*"), (axes[1], "P_star", "P*")):
        for al in alphas:
            xs, dys = [], []
            for n in ns:
                if (n, 0.0) in rows and (n, al) in rows:
                    xs.append(n)
                    dys.append(float(rows[(n, al)][key]) - float(rows[(n, 0.0)][key]))
            # normalize by the alpha=0 value
            xs2, rel = [], []
            for n, d in zip(xs, dys):
                base = float(rows[(n, 0.0)][key])
                xs2.append(n)
                rel.append(100.0 * d / base)
            ax.plot(xs2, rel, "o-", label="alpha=%.2f" % al)
        ax.set_xlabel("n (forme globale)")
        ax.set_ylabel("effet relatif de l'ondulation  (d%s/%s, %%)" % (lab, lab))
        ax.axhline(0, color="gray", lw=0.8)
        ax.grid(True, alpha=0.25)
        ax.legend()
    fig.suptitle("Effet de l'ondulation du fond vs forme globale (Re=%g)" % Re)
    fig.tight_layout()
    out = os.path.join(FIGDIR, "effect_vs_n_Re%g.png" % Re)
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out


def fig_interaction(inter, Re):
    pairs = inter["pairs"]
    names = list(pairs.keys())
    labels, ID, IP = [], [], []
    for k in names:
        labels.append(k.replace("_", " "))
        ID.append(pairs[k]["D_star"]["I"])
        IP.append(pairs[k]["P_star"]["I"])
    x = np.arange(len(names))
    fig, ax = plt.subplots(figsize=(9, 4))
    ax.bar(x - 0.2, ID, 0.4, label="I_D")
    ax.bar(x + 0.2, IP, 0.4, label="I_P")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=30, ha="right", fontsize=8)
    ax.axhline(0, color="gray", lw=0.8)
    ax.set_ylabel("I_Q = dQ(n_b) - dQ(n_a)")
    ax.set_title("Difference de differences (Re=%g, N=%d)" % (Re, inter["N"]))
    ax.grid(True, alpha=0.25, axis="y")
    ax.legend()
    fig.tight_layout()
    out = os.path.join(FIGDIR, "interaction_Re%g.png" % Re)
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out


def fig_fields(Re, N, cases):
    ncol = len(cases)
    fig, axes = plt.subplots(2, ncol, figsize=(3.2 * ncol, 6), squeeze=False)
    for kc, cid in enumerate(cases):
        path = os.path.join(FIELDS_DIR, cid + ".npz")
        if not os.path.exists(path):
            continue
        d = np.load(path)
        ux, uy, fluid = d["ux"], d["uy"], d["fluid"]
        xc, yc = d["xc"], d["yc"]
        psi = np.zeros_like(ux)
        U = np.where(fluid, ux, 0.0)
        Nn = ux.shape[0]
        dx = xc[1] - xc[0]
        for j in range(Nn):
            col = np.nonzero(fluid[:, j])[0]
            if col.size == 0:
                continue
            i0 = col[0]
            acc = 0.0
            for i in range(i0 + 1, Nn):
                if not fluid[i, j]:
                    break
                acc += 0.5 * (U[i, j] + U[i - 1, j]) * dx
                psi[i, j] = acc
        X, Y = np.meshgrid(xc, yc)
        sp = np.where(fluid, np.sqrt(ux ** 2 + uy ** 2), np.nan)
        im0 = axes[0][kc].pcolormesh(X, Y, sp, shading="auto", cmap="viridis")
        axes[0][kc].set_title(cid.replace("se_", "").replace("_Nw1", ""), fontsize=8)
        axes[0][kc].set_aspect("equal")
        fig.colorbar(im0, ax=axes[0][kc], shrink=0.7)
        im1 = axes[1][kc].pcolormesh(X, Y, np.where(fluid, psi, np.nan),
                                     shading="auto", cmap="RdBu_r")
        axes[1][kc].set_aspect("equal")
        fig.colorbar(im1, ax=axes[1][kc], shrink=0.7)
    axes[0][0].set_ylabel("|u|")
    axes[1][0].set_ylabel("psi")
    fig.suptitle("Champs (Re=%g, N=%d)" % (Re, N))
    fig.tight_layout()
    out = os.path.join(FIGDIR, "fields_Re%g_N%d.png" % (Re, N))
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--Re", type=float, default=400)
    ap.add_argument("--N", type=int, default=128)
    args = ap.parse_args()
    rows = load_pilot_table(args.Re, args.N, 0.10)
    f1 = fig_effect_vs_n(rows, args.Re)
    inter_file = os.path.join(SCALARS_DIR, "interaction_Re%g_N%d.json" % (args.Re, args.N))
    f2 = None
    if os.path.exists(inter_file):
        with open(inter_file, encoding="utf-8") as fh:
            inter = json.load(fh)
        f2 = fig_interaction(inter, args.Re)
    ns = sorted({k[0] for k in rows})
    cases = []
    for n in ns:
        for al in (0.0, 0.10):
            cid = "se_n%g_w%.2f_Nw1_Re%d_N%d" % (n, al, int(args.Re), args.N)
            if os.path.exists(os.path.join(FIELDS_DIR, cid + ".npz")):
                cases.append(cid)
    f3 = fig_fields(args.Re, args.N, cases) if cases else None
    print("[ok]", f1)
    if f2:
        print("[ok]", f2)
    if f3:
        print("[ok]", f3)
    return 0


if __name__ == "__main__":
    sys.exit(main())
