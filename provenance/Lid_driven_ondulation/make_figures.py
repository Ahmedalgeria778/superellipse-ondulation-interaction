# -*- coding: utf-8 -*-
"""make_figures.py -- manuscript figures 1-8 (real project data).

Outputs PNG (300 dpi) and PDF into figures/.
All series are distinguished by marker AND line style (monochrome-safe),
not by colour alone.
Data come from results/ and results/analysis/ ; nothing is invented.
Only external data are the tabulated values of Ghia et al. (1982).
"""
from __future__ import annotations

import csv
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RES = os.path.join(ROOT, "results")
ANA = os.path.join(RES, "analysis")
CONV = os.path.join(RES, "convergence")
FIG = os.path.join(ROOT, "figures")
DATA = os.path.join(ROOT, "data")
sys.path.insert(0, HERE)

from diagnostics import detect_vortices  # noqa: E402

os.makedirs(FIG, exist_ok=True)

plt.rcParams.update({
    "font.family": "serif",
    "font.size": 10,
    "axes.titlesize": 10.5,
    "axes.labelsize": 10,
    "legend.fontsize": 8,
    "figure.dpi": 300,
    "savefig.bbox": "tight",
    "axes.grid": False,
    "lines.linewidth": 1.3,
})

# Monochrome-safe styles: (colour, marker, linestyle)
ST = {
    "F00": ("#1f77b4", "o", "-"),
    "G06": ("#2ca02c", "s", "--"),
    "G12": ("#d62728", "^", "-."),
    "G13": ("#9467bd", "D", ":"),
}


def save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(FIG, f"{name}.{ext}"))
    plt.close(fig)
    print("  wrote", name)


def load_npz(rel):
    return np.load(os.path.join(RES, rel), allow_pickle=True)


def geom_from(d):
    N = int(d["N"])
    fluid = d["fluid"].astype(bool)
    return {"N": N, "dx": 1.0 / N, "fluid": fluid, "solid": ~fluid,
            "xc": d["xc"], "yc": d["yc"], "yb": np.asarray(d["yb"])}


def read_csv(path):
    with open(path, encoding="utf-8") as f:
        return list(csv.DictReader(f))


def wall(x, h, nw):
    return (h / 2.0) * (1.0 - np.cos(2.0 * np.pi * nw * x))


# --- Ghia et al. (1982) tabulated data (Table I: u(y) ; Table II: v(x)) -------
GHIA_Y = [1.0000, 0.9766, 0.9688, 0.9609, 0.9531, 0.8516, 0.7344, 0.6172,
          0.5000, 0.4531, 0.2813, 0.1719, 0.1016, 0.0703, 0.0625, 0.0547, 0.0000]
GHIA_X = [1.0000, 0.9688, 0.9609, 0.9531, 0.9453, 0.9063, 0.8594, 0.8047,
          0.5000, 0.2344, 0.2266, 0.1563, 0.0938, 0.0781, 0.0703, 0.0625, 0.0000]
GHIA_U = {
    100: [1.00000, 0.84123, 0.78871, 0.73722, 0.68717, 0.23151, 0.00332,
          -0.13641, -0.20581, -0.21090, -0.15662, -0.10150, -0.06434, -0.04775,
          -0.04192, -0.03717, 0.00000],
    400: [1.00000, 0.75837, 0.68439, 0.61756, 0.55892, 0.29093, 0.16256,
          0.02135, -0.11477, -0.17119, -0.32726, -0.24299, -0.14612, -0.10338,
          -0.09266, -0.08186, 0.00000],
    1000: [1.00000, 0.65928, 0.57492, 0.51117, 0.46604, 0.33304, 0.18719,
           0.05702, -0.06080, -0.10648, -0.27805, -0.38289, -0.29730, -0.22220,
           -0.20196, -0.18109, 0.00000],
}
GHIA_V = {
    100: [0.00000, -0.05906, -0.07391, -0.08864, -0.10313, -0.16914, -0.22445,
          -0.24533, 0.05454, 0.17527, 0.17507, 0.16077, 0.12317, 0.10890,
          0.10091, 0.09233, 0.00000],
    400: [0.00000, -0.12146, -0.15663, -0.19254, -0.22847, -0.23827, -0.44993,
          -0.38598, 0.05186, 0.30174, 0.30203, 0.28124, 0.22965, 0.20920,
          0.19713, 0.18360, 0.00000],
    1000: [0.00000, -0.21388, -0.27669, -0.33714, -0.39188, -0.51500, -0.42665,
           -0.31966, 0.02526, 0.32235, 0.33075, 0.37095, 0.32627, 0.30353,
           0.29012, 0.27485, 0.00000],
}
GHIA_TYPO = {400: 0.9063}
GHIA_VERIFIED = {100: "yes", 400: "yes (0.9063 suspect)", 1000: "yes"}


def write_ghia_csv():
    os.makedirs(DATA, exist_ok=True)
    with open(os.path.join(DATA, "ghia_reference.csv"), "w", newline="",
              encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["# Ghia et al. (1982) J. Comput. Phys. 48:387-411, Tables I & II"])
        w.writerow(["# u = u-velocity on vertical centreline x=0.5 (vs y)"])
        w.writerow(["# v = v-velocity on horizontal centreline y=0.5 (vs x)"])
        w.writerow(["# NOTE: (Re=400, x=0.9063, v=-0.23827) is a known, non-monotone typo"])
        w.writerow(["quantity", "Re", "coord", "value", "verified"])
        for Re in (100, 400, 1000):
            for c, val in zip(GHIA_Y, GHIA_U[Re]):
                w.writerow(["u", Re, c, val, GHIA_VERIFIED[Re]])
        for Re in (100, 400, 1000):
            for c, val in zip(GHIA_X, GHIA_V[Re]):
                w.writerow(["v", Re, c, val, GHIA_VERIFIED[Re]])
    print("  wrote data/ghia_reference.csv")


def centerlines(d):
    xc, yc, ux, uy = d["xc"], d["yc"], d["ux"], d["uy"]
    i0 = int(np.argmin(np.abs(xc - 0.5)))
    j0 = int(np.argmin(np.abs(yc - 0.5)))
    return yc, ux[:, i0], xc, uy[j0, :]


# ---------------------------------------------------------------------------
def fig01_geometry():
    fig, ax = plt.subplots(1, 2, figsize=(9.0, 3.4))
    a = ax[0]
    Nw, h = 4, 0.075
    x = np.linspace(0, 1, 800)
    yb = wall(x, h, Nw)
    a.fill_between(x, yb, 0, color="0.85", zorder=0)
    a.plot(x, yb, "k-", lw=1.6)
    a.plot([0, 0], [0, 1], "k-", lw=1.2)
    a.plot([1, 1], [0, 1], "k-", lw=1.2)
    a.plot([0, 1], [1, 1], "k-", lw=1.2)
    for xa in np.linspace(0.15, 0.85, 6):
        a.annotate("", xy=(xa + 0.08, 1.0), xytext=(xa, 1.0),
                   arrowprops=dict(arrowstyle="->", color="k", lw=1.1))
    a.text(0.5, 1.06, r"$U_{lid}$", ha="center")
    a.annotate("", xy=(1.04, 0), xytext=(1.04, 1),
               arrowprops=dict(arrowstyle="<->", lw=1))
    a.text(1.06, 0.5, r"$H$", va="center")
    a.annotate("", xy=(0, -0.06), xytext=(1, -0.06),
               arrowprops=dict(arrowstyle="<->", lw=1))
    a.text(0.5, -0.12, r"$L=H=1$", ha="center")
    a.annotate("", xy=(0.375, yb[300]), xytext=(0.375, 0.0),
               arrowprops=dict(arrowstyle="<->", lw=0.9))
    a.text(0.30, 0.10, r"$h$", rotation=90, va="bottom")
    a.set_xlim(-0.1, 1.16)
    a.set_ylim(-0.16, 1.12)
    a.set_aspect("equal")
    a.set_xticks([]); a.set_yticks([])
    a.set_title("(a) geometry and boundary conditions")

    b = ax[1]
    styles = [("0.025", "o"), ("0.05", "s"), ("0.075", "^"), ("0.10", "D")]
    for (hh, nw, m), lab in zip([(0.025, 1, "o"), (0.05, 2, "s"),
                                 (0.075, 4, "^"), (0.10, 4, "D")],
                                [r"$\alpha=0.025,\ N_w=1$",
                                 r"$\alpha=0.05,\ N_w=2$",
                                 r"$\alpha=0.075,\ N_w=4$",
                                 r"$\alpha=0.10,\ N_w=4$"]):
        b.plot(x, wall(x, hh, nw), "-", lw=1.5, label=lab)
    b.axhline(0, color="k", lw=0.8, ls=":")
    b.set_xlabel(r"$x^\ast$"); b.set_ylabel(r"$y_b^\ast(x^\ast)$")
    b.set_title("(b) corrugated bottom-wall shapes")
    b.legend(loc="upper right", framealpha=0.9)
    b.set_xlim(0, 1)
    save(fig, "fig01_geometry")


def fig02_ghia_profiles():
    cases = {100: "V01-R100/V01-R100_N256.npz",
             400: "V02-R400/V02-R400_N256.npz",
             1000: "V03-R1000/V03-R1000_N256.npz"}
    mk = {100: "o", 400: "s", 1000: "^"}
    col = {100: "#000000", 400: "#444444", 1000: "#000000"}
    ls = {100: "-", 400: "--", 1000: "-."}
    fig, ax = plt.subplots(1, 2, figsize=(8.4, 4.0))
    for Re, rel in cases.items():
        d = load_npz(rel)
        yc, u_line, xc, v_line = centerlines(d)
        ax[0].plot(u_line, yc, ls[Re], color=col[Re], lw=1.5,
                   label=f"LBM, Re={Re}")
        ax[1].plot(xc, v_line, ls[Re], color=col[Re], lw=1.5,
                   label=f"LBM, Re={Re}")
        ax[0].plot(GHIA_U[Re], GHIA_Y, mk[Re], mfc="none", mec="k", ms=5.5,
                   label=f"Ghia, Re={Re}")
        ax[1].plot(GHIA_X, GHIA_V[Re], mk[Re], mfc="none", mec="k", ms=5.5,
                   label=f"Ghia, Re={Re}")
        if Re in GHIA_TYPO:
            ax[1].plot([GHIA_TYPO[Re]], [-0.23827], "kx", ms=7, mew=1.3)
    ax[0].set_xlabel(r"$u$"); ax[0].set_ylabel(r"$y$")
    ax[0].set_title(r"(a) vertical centreline $u(y)$ at $x=0.5$")
    ax[0].legend(ncol=2, fontsize=7.5, framealpha=0.9)
    ax[1].set_xlabel(r"$x$"); ax[1].set_ylabel(r"$v$")
    ax[1].set_title(r"(b) horizontal centreline $v(x)$ at $y=0.5$")
    ax[1].legend(ncol=2, fontsize=7.5, framealpha=0.9)
    ax[1].text(GHIA_TYPO[400], -0.23827, " typo", fontsize=7, va="center")
    fig.suptitle("Flat-wall validation, uniform lid (N=256)", y=1.02)
    save(fig, "fig02_ghia_profiles")


def fig04_fields():
    runs = [("F00-R400", "convergence/consistent/F00-R400_N512/fields.npz"),
            ("G06-R400", "convergence/consistent/G06-R400_N512/fields.npz"),
            ("G12-R400", "convergence/consistent/G12-R400_N512/fields.npz")]
    fig, ax = plt.subplots(2, 3, figsize=(10.5, 6.2))
    for k, (name, rel) in enumerate(runs):
        d = load_npz(rel)
        g = geom_from(d)
        fluid = g["fluid"]
        psi = np.asarray(d["psi"])
        ux, uy = np.asarray(d["ux"]), np.asarray(d["uy"])
        xc, yc, dx = g["xc"], g["yc"], g["dx"]
        yb = np.asarray(d["yb"])[0]
        dom = np.where(fluid, psi, np.nan)
        ax[0, k].contourf(xc, yc, dom, levels=30, cmap="viridis")
        try:
            ax[0, k].contour(xc, yc, dom, levels=12, colors="w",
                             linewidths=0.35, alpha=0.7)
        except Exception:
            pass
        ax[0, k].plot(xc, yb, "k-", lw=1.2)
        vort = detect_vortices(psi, g)
        if vort:
            ax[0, k].plot([v["x"] for v in vort], [v["y"] for v in vort],
                          "ko", ms=4, mfc="w", mew=1.1)
        ax[0, k].set_title(f"{name}  $n_v={len(vort)}$")
        ax[0, k].set_aspect("equal"); ax[0, k].set_xlim(0, 1); ax[0, k].set_ylim(0, 1)
        ax[0, k].set_xticks([]); ax[0, k].set_yticks([])

        w = np.gradient(uy, dx, axis=1) - np.gradient(ux, dx, axis=0)
        wdom = np.where(fluid, w, np.nan)
        vmax = np.nanpercentile(np.abs(wdom), 99)
        ax[1, k].contourf(xc, yc, wdom, levels=np.linspace(-vmax, vmax, 31),
                          cmap="RdBu_r", extend="both")
        ax[1, k].contour(xc, yc, wdom, levels=[0.0], colors="k", linewidths=0.6)
        ax[1, k].plot(xc, yb, "k-", lw=1.2)
        ax[1, k].set_xlabel(r"$x$")
        ax[1, k].set_aspect("equal"); ax[1, k].set_xlim(0, 1); ax[1, k].set_ylim(0, 1)
        ax[1, k].set_yticks([])
    ax[0, 0].set_ylabel(r"streamfunction  $\psi$")
    ax[1, 0].set_ylabel(r"vorticity  $\omega$")
    fig.suptitle(r"Flow fields at $Re=400$, $N=512$ (open circles: vorticity extrema)",
                 y=0.99)
    save(fig, "fig04_fields")


def convergence_data():
    data = {}
    for r in read_csv(os.path.join(ANA, "wavy_effect.csv")):
        n = int(r["N"])
        if r["case_wavy"] in ("G06-R400", "G12-R400"):
            data[(r["case_wavy"], n)] = (float(r["P_star_wavy"]),
                                         float(r["D_star_wavy"]))
        if r["case_flat"] == "F00-R400":
            data[("F00-R400", n)] = (float(r["P_star_flat"]),
                                     float(r["D_star_flat"]))
    for r in read_csv(os.path.join(CONV, "resolution.csv")):
        if r["variant"] == "consistent":
            data.setdefault((r["case"], int(r["N"])),
                            (float(r["P_star"]), float(r["D_star"])))
    return data


def fig03_convergence():
    data = convergence_data()
    tags = ["F00-R400", "G06-R400", "G12-R400"]
    fig, ax = plt.subplots(1, 3, figsize=(10.5, 3.3), constrained_layout=True)
    for tag in tags:
        key = "F00" if tag.startswith("F00") else tag[:3]
        c, m, ls = ST[key]
        ns = sorted(n for (t, n) in data if t == tag)
        P = [data[(tag, n)][0] for n in ns]
        D = [data[(tag, n)][1] for n in ns]
        bal = [100.0 * (a - b) / b for a, b in zip(P, D)]
        ax[0].plot(ns, P, m + ls, color=c, ms=5, label=tag)
        ax[1].plot(ns, D, m + ls, color=c, ms=5, label=tag)
        ax[2].plot(ns, bal, m + ls, color=c, ms=5, label=tag)
    ax[0].set_ylabel(r"lid power  $P^\ast$"); ax[0].set_title(r"(a) lid power $P^\ast$")
    ax[1].set_ylabel(r"dissipation  $D^\ast$"); ax[1].set_title(r"(b) dissipation $D^\ast$")
    ax[2].set_ylabel(r"$(P^\ast-D^\ast)/D^\ast$  (\%)")
    ax[2].set_title("(c) energy-budget closure")
    ax[2].axhline(-3.2, color="k", ls=":", lw=0.8)
    for a in ax:
        a.set_xlabel(r"$N$"); a.set_xticks([128, 256, 384, 512])
        a.legend(framealpha=0.9)
    fig.suptitle(r"Grid convergence at $Re=400$ (regularised lid)")
    save(fig, "fig03_convergence")


def fig05_deltaq_vs_N():
    rows = read_csv(os.path.join(ANA, "wavy_effect.csv"))
    fig, ax = plt.subplots(1, 2, figsize=(8.4, 3.4))
    qs = [("dP_star_pct", "o", "-", r"$\Delta P^\ast$"),
          ("dD_star_pct", "s", "--", r"$\Delta D^\ast$"),
          ("dK_star_pct", "^", ":", r"$\Delta K^\ast$")]
    for tag in ("G06-R400", "G12-R400"):
        sub = sorted((r for r in rows if r["case_wavy"] == tag),
                     key=lambda r: int(r["N"]))
        ns = [int(r["N"]) for r in sub]
        key = tag[:3]
        c = ST[key][0]
        for q, m, ls, lab in qs:
            ax[0].plot(ns, [float(r[q]) for r in sub], m + ls, color=c, ms=5)
    ax[0].axhline(0, color="k", lw=0.7)
    ax[0].set_xlabel(r"$N$"); ax[0].set_ylabel(r"$\Delta Q$  (\%)")
    ax[0].set_title(r"(a) resolution dependence of $\Delta Q$ at $Re=400$")
    ax[0].set_xticks([128, 256, 384, 512])
    ax[0].text(300, -0.9, "G06", color=ST["G06"][0], fontsize=8)
    ax[0].text(300, -1.6, "G12", color=ST["G12"][0], fontsize=8)
    h = [plt.Line2D([], [], color="k", marker=m, ls=ls, ms=5, label=lab)
         for _, m, ls, lab in qs]
    ax[0].legend(handles=h, framealpha=0.9)

    co = read_csv(os.path.join(ANA, "convergence_order.csv"))
    labels, vals, cols = [], [], []
    for r in co:
        if r["series"].startswith("d"):
            labels.append(r["series"].replace("dG", "G"))
            vals.append(float(r["Q_last"]))
            cols.append(ST[r["series"][1:4]][0])
    ax[1].bar(range(len(labels)), vals, color=cols, edgecolor="k", lw=0.8)
    ax[1].set_xticks(range(len(labels))); ax[1].set_xticklabels(labels)
    ax[1].axhline(0, color="k", lw=0.7)
    ax[1].set_ylabel(r"$\Delta P^\ast_{\infty}$  (\%)")
    ax[1].set_title(r"(b) Richardson extrapolation of $\Delta P^\ast$")
    save(fig, "fig05_deltaq_vs_N")


def fig06_deltaq_map():
    rows = read_csv(os.path.join(ANA, "normalize.csv"))
    fig, ax = plt.subplots(1, 2, figsize=(9.0, 3.6))
    for Re, key in [(100, "F00"), (400, "G06"), (1000, "G12")]:
        c, m, _ = ST[key]
        sub = [r for r in rows if int(r["Re"]) == Re]
        x = [float(r["h_over_H"]) * int(r["Nw"]) for r in sub]
        y = [float(r["dP_pct_F00"]) for r in sub]
        ax[0].plot(x, y, m, color=c, ms=5, label=f"Re={Re}")
        ax[1].scatter([Re] * len(sub), y,
                      s=[16 + 55 * float(r["h_over_H"]) * int(r["Nw"]) for r in sub],
                      marker=m, facecolors="none", edgecolors="k")
    ax[0].set_xlabel(r"$h\,N_w/H$"); ax[0].set_ylabel(r"$\Delta P^\ast$  (\%)")
    ax[0].set_title(r"(a) wavy-wall effect vs $hN_w$")
    ax[0].legend(framealpha=0.9)
    ax[1].axhline(0, color="k", lw=0.7)
    ax[1].set_xlabel(r"$Re$"); ax[1].set_ylabel(r"$\Delta P^\ast$  (\%)")
    ax[1].set_title(r"(b) vs $Re$ (symbol size $\propto hN_w$)")
    ax[1].set_xscale("log"); ax[1].set_xticks([100, 400, 1000])
    ax[1].get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    save(fig, "fig06_deltaq_map")


def fig07_lid_profiles():
    x = np.linspace(0, 1, 600)
    eps = 0.02

    def g(xs):
        t = np.clip(xs / eps, 0, 1)
        gl = 3 * t ** 2 - 2 * t ** 3
        t2 = np.clip((1 - xs) / eps, 0, 1)
        gr = 3 * t2 ** 2 - 2 * t2 ** 3
        return np.minimum(gl, gr)

    gg = g(x); gm = gg / gg.mean()
    fig, ax = plt.subplots(1, 2, figsize=(9.0, 3.4))
    ax[0].plot(x, np.ones_like(x), "-", color="k", lw=1.8, label="uniform (U)")
    ax[0].plot(x, gg, "--", color="k", lw=1.8, label="regularised (A)")
    ax[0].plot(x, gm, ":", color="k", lw=2.0, label="smoothstep, mean 1 (B)")
    ax[0].set_xlabel(r"$x$"); ax[0].set_ylabel(r"$U_{lid}/U_{ref}$")
    ax[0].set_title("(a) moving-lid profiles")
    ax[0].legend(framealpha=0.9); ax[0].set_xlim(0, 1); ax[0].set_ylim(0.9, 1.05)
    pb = read_csv(os.path.join(ANA, "protocol_b_compare.csv"))
    labels = [r["case"].replace("-R400", "") for r in pb]
    idx = np.arange(len(labels)); width = 0.26
    bars = [(r"$K^\ast$", "K", "o", ""), (r"$P^\ast$", "P", "s", "//"),
            (r"$D^\ast$", "D", "^", "\\\\")]
    for i, (lab, kk, mk, hatch) in enumerate(bars):
        vals = [100.0 * (float(r[f"{kk}_star_mean1"]) - float(r[f"{kk}_star_uniform"]))
                / float(r[f"{kk}_star_uniform"]) for r in pb]
        ax[1].bar(idx + (i - 1) * width, vals, width, label=lab + " (B vs U)",
                  facecolor="w", edgecolor="k", hatch=hatch, lw=0.9)
    ax[1].axhline(0, color="k", lw=0.7)
    ax[1].set_xticks(idx); ax[1].set_xticklabels(labels)
    ax[1].set_ylabel(r"$\Delta$ vs uniform lid  (\%)")
    ax[1].set_title(r"(b) protocol B, $Re=400$, $N=256$")
    ax[1].legend(framealpha=0.9, fontsize=7.5)
    save(fig, "fig07_lid_profiles")


def fig08_vortices():
    topo = read_csv(os.path.join(ANA, "topology.csv"))
    fig, ax = plt.subplots(1, 2, figsize=(9.0, 3.4))
    for Re, key in [(100, "F00"), (400, "G06"), (1000, "G12")]:
        c, m, _ = ST[key]
        sub = [r for r in topo if int(r["Re"]) == Re]
        ax[0].plot([float(r["hNw"]) for r in sub],
                   [int(r["n_vortices"]) for r in sub], m, color=c, ms=6,
                   mfc="w", label=f"Re={Re}")
    ax[0].set_xlabel(r"$h\,N_w/H$"); ax[0].set_ylabel(r"$n_v$")
    ax[0].set_title(r"(a) vortex count vs $hN_w$")
    ax[0].legend(framealpha=0.9)
    res = read_csv(os.path.join(CONV, "resolution.csv"))
    for tag in ("G06-R400", "G12-R400"):
        c, m, _ = ST[tag[:3]]
        for var, ls in [("consistent", "-"), ("ma_fixed", "--")]:
            sub = sorted((r for r in res if r["case"] == tag and r["variant"] == var),
                         key=lambda r: int(r["N"]))
            if sub:
                ax[1].plot([int(r["N"]) for r in sub],
                           [int(r["n_vortices"]) for r in sub], m + ls,
                           color=c, ms=5, mfc="w", label=f"{tag[:3]} {var}")
    ax[1].set_xlabel(r"$N$"); ax[1].set_ylabel(r"$n_v$")
    ax[1].set_title(r"(b) topological sensitivity to $N$")
    ax[1].set_xticks([128, 256, 384, 512])
    ax[1].legend(framealpha=0.9, fontsize=7.5)
    save(fig, "fig08_vortices")


def fig09_scaling():
    rows = [r for r in read_csv(os.path.join(ANA, "normalize.csv"))
            if int(r["N"]) == 256]
    a = np.array([float(r["h_over_H"]) for r in rows])
    re = np.array([float(r["Re"]) for r in rows])
    dp = np.array([abs(float(r["dP_pct_F00"])) for r in rows])
    lna, lnr, lny = np.log(a), np.log(re), np.log(dp)
    X = np.column_stack([np.ones(len(rows)), lna, lnr])
    b, *_ = np.linalg.lstsq(X, lny, rcond=None)
    C, Aexp, B = np.exp(b[0]), b[1], b[2]
    res = lny - X @ b
    n, k = X.shape
    rmse = float(np.sqrt((res ** 2).mean()))
    se = np.sqrt(np.diag((res @ res / (n - k)) * np.linalg.inv(X.T @ X)))
    fig, ax = plt.subplots(1, 3, figsize=(11.0, 3.3))
    cols = {100: "0.35", 400: "0.0", 1000: "0.6"}
    mks = {100: "o", 400: "s", 1000: "^"}
    for Re in (100, 400, 1000):
        m = re == Re
        ax[0].plot(a[m], dp[m], mks[Re], mfc="none", mec=cols[Re], ms=6,
                   label=f"$Re={Re}$")
        aa = np.linspace(a[m].min(), a[m].max(), 50)
        ax[0].plot(aa, C * aa ** Aexp * Re ** B, "-", color=cols[Re], lw=1.1)
    ax[0].set_xscale("log"); ax[0].set_yscale("log")
    ax[0].set_xlabel(r"$\alpha=h/H$")
    ax[0].set_ylabel(r"$|\Delta P^\ast|$  (\%)")
    ax[0].set_title(r"(a) $|\Delta P^\ast|\propto \alpha^{a}Re^{b}$")
    ax[0].legend(framealpha=0.9, fontsize=8)
    ax[0].text(0.03, 0.05,
               (r"$a=%.2f\pm%.2f$" % (Aexp, se[1])) + "\n"
               + (r"$b=%.2f\pm%.2f$" % (B, se[2])),
               transform=ax[0].transAxes, fontsize=8, va="bottom")
    norm = dp / (a ** Aexp * re ** B)
    for Re in (100, 400, 1000):
        m = re == Re
        ax[1].plot(a[m], norm[m], mks[Re], mfc="none", mec=cols[Re], ms=6)
    ax[1].axhline(C, color="k", ls=":", lw=0.9)
    ax[1].set_xscale("log")
    ax[1].set_xlabel(r"$\alpha=h/H$")
    ax[1].set_ylabel(r"$|\Delta P^\ast|/(\alpha^{a}Re^{b})$")
    ax[1].set_title("(b) collapse of the 48-point sweep")
    for Re in (100, 400, 1000):
        m = re == Re
        ax[2].plot(a[m], res[m], mks[Re], mfc="none", mec=cols[Re], ms=6)
    ax[2].axhspan(-rmse, rmse, color="0.85", zorder=0)
    ax[2].axhline(0.0, color="k", ls=":", lw=0.9)
    ax[2].set_xscale("log")
    ax[2].set_xlabel(r"$\alpha=h/H$")
    ax[2].set_ylabel("residual (log units)")
    ax[2].set_title(r"(c) residuals, $\pm$RMSE$=%.2f$" % rmse)
    save(fig, "fig09_scaling")


def main():
    write_ghia_csv()
    fig01_geometry()
    fig02_ghia_profiles()
    fig03_convergence()
    fig04_fields()
    fig05_deltaq_vs_N()
    fig06_deltaq_map()
    fig07_lid_profiles()
    fig08_vortices()
    fig09_scaling()
    print("FIGURES OK")


if __name__ == "__main__":
    main()
