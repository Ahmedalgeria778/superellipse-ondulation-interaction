# -*- coding: utf-8 -*-
"""
solver/make_paper_figs.py
=========================
Regenerates the paper figures that depend on the primary contrast
I_D(2,4) at alpha = 0, taken directly from the JSON scalar records written
by the fixed-tau campaign (results/scalars/se_*_tol1e-08.json).

Figures produced (written to results/paper_figures/, then copied into the
manuscript tree):

  fig6_resolution.pdf   I_D(2,4) vs N = 128/256/384 for Re = 100/400/1000
  fig9_resweep.pdf      (a) log-log |I_D(2,4)| vs Re (N=128 sweep)
                        (b) local exponent gamma_loc = dln|I_D|/dln Re
  figS_contrasts_res.pdf (a) I_D(2,4) vs Re at each N
                         (b) all successive transitions at Re=400
                         (c) estimator discrepancy delta_IP/ID(2,4)

All quantities are the pure alpha = 0 contrasts (each geometry compared
with its own flat-... no: n=2 vs n=4, both at alpha = 0), matching the
definitions and numbers quoted in Sections IV.C, IV.D and the SM.
"""

from __future__ import annotations

import json
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
SCALARS = os.path.join(_ROOT, "results", "scalars")
OUT = os.path.join(_ROOT, "results", "paper_figures")
os.makedirs(OUT, exist_ok=True)

RE_FIXED = (100, 400, 1000)
NS_FIXED = (128, 256, 384)
RE_SWEEP = list(range(100, 1001, 50))  # 100..1000 ; gamma_loc ends at 950 (central diff)


def rec(n, a, Re, N):
    f = os.path.join(SCALARS, "se_n%g_w%.2f_Nw1_Re%d_N%d_tol1e-08.json" % (n, a, Re, N))
    if not os.path.exists(f):
        return None
    with open(f, encoding="utf-8") as fh:
        return json.load(fh)


def contrast(n_a, n_b, a, Re, N, key="D_star"):
    ra, rb = rec(n_a, a, Re, N), rec(n_b, a, Re, N)
    if ra is None or rb is None:
        return None
    return 100.0 * (rb[key] - ra[key]) / ra[key]


# --------------------------------------------------------------- fig 6
def fig6():
    fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.9))
    for ax, Re in zip(axes, RE_FIXED):
        ys, xs = [], []
        for N in NS_FIXED:
            v = contrast(2, 4, 0.0, Re, N)
            if v is not None:
                xs.append(N)
                ys.append(v)
        ax.plot(xs, ys, "o-", color="#1f77b4", ms=7, lw=1.8)
        for x, y in zip(xs, ys):
            ax.annotate("%.4f" % y, (x, y), textcoords="offset points",
                        xytext=(0, 8), ha="center", fontsize=7.5)
        ax.set_xscale("log", base=2)
        ax.set_xticks(list(NS_FIXED))
        ax.set_xticklabels([str(N) for N in NS_FIXED])
        ax.set_xlabel("grid $N$")
        ax.set_ylabel(r"$I_D(2,4)$  $[\%]$")
        ax.set_title(r"$Re=%d$" % Re)
        ax.grid(alpha=0.3)
    fig.suptitle(r"Resolution sensitivity of $I_D(2,4)$ at $\alpha=0$")
    fig.tight_layout()
    p = os.path.join(OUT, "fig6_resolution.pdf")
    fig.savefig(p)
    plt.close(fig)
    return p


# --------------------------------------------------------------- fig 9
def sweep_series():
    Res, Is = [], []
    for Re in RE_SWEEP:
        v = contrast(2, 4, 0.0, Re, 128)
        if v is not None:
            Res.append(Re)
            Is.append(v)
    return np.array(Res, float), np.array(Is, float)


def gamma_loc(Re, I):
    lnRe, lnI = np.log(Re), np.log(I)
    g, Rg = [], []
    for k in range(1, len(Re) - 1):
        g.append((lnI[k + 1] - lnI[k - 1]) / (lnRe[k + 1] - lnRe[k - 1]))
        Rg.append(Re[k])
    return np.array(Rg), np.array(g)


def fig9():
    Re, I = sweep_series()
    Rg, g = gamma_loc(Re, I)
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.2))
    ax = axes[0]
    ax.loglog(Re, I, "o-", color="#d62728", ms=6, lw=1.6)
    ax.set_xlabel("$Re$")
    ax.set_ylabel(r"$I_D(2,4)$  $[\%]$")
    ax.set_title(r"(a) growth of the primary contrast ($N=128$)")
    ax.grid(alpha=0.3, which="both")

    ax = axes[1]
    ax.plot(Rg, g, "s-", color="#2ca02c", ms=6, lw=1.8)
    ax.axhline(1.0, color="gray", lw=1.0, ls="--")
    # crossing of unity by central difference
    xc = None
    for k in range(len(g) - 1):
        if (g[k] - 1.0) * (g[k + 1] - 1.0) <= 0:
            xc = Rg[k] + (1.0 - g[k]) * (Rg[k + 1] - Rg[k]) / (g[k + 1] - g[k])
            break
    if xc is not None:
        ax.axvline(xc, color="0.4", lw=1.0, ls=":")
        ax.annotate(r"$\gamma_{\rm loc}=1$ near $Re\approx%.0f$" % xc,
                    (xc, 1.6), fontsize=8, color="0.25")
    ax.set_xscale("log")
    ax.set_xlabel("$Re$")
    ax.set_ylabel(r"$\gamma_{\rm loc}=d\ln|I_D|/d\ln Re$")
    ax.set_title(r"(b) local exponent (central differences)")
    ax.grid(alpha=0.3, which="both")
    fig.tight_layout()
    p = os.path.join(OUT, "fig9_resweep.pdf")
    fig.savefig(p)
    plt.close(fig)
    return p, Re, I, Rg, g, xc


# --------------------------------------------------------------- fig S
def figS():
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.0))

    ax = axes[0]
    mk = {128: ("s", "none"), 256: ("o", "none"), 384: ("D", "#1f77b4")}
    for N in NS_FIXED:
        xs, ys = [], []
        for Re in RE_FIXED:
            v = contrast(2, 4, 0.0, Re, N)
            if v is not None:
                xs.append(Re)
                ys.append(v)
        if xs:
            m, fc = mk[N]
            ax.plot(xs, ys, m + ("-" if N == 384 else "--"), color="#1f77b4",
                    ms=7, lw=1.6, mfc=fc, mec="#1f77b4",
                    label="$N=%d$" % N)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("$Re$")
    ax.set_ylabel(r"$I_D(2,4)$  $[\%]$")
    ax.set_title(r"(a) primary contrast")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, which="both")

    ax = axes[1]
    trans = [(2, 4), (4, 6), (6, 8)]
    for (na, nb), c in zip(trans, ("#1f77b4", "#ff7f0e", "#d62728")):
        xs, ys = [], []
        for N in NS_FIXED:
            vd = contrast(na, nb, 0.0, 400, N, "D_star")
            if vd is not None:
                xs.append(N)
                ys.append(abs(vd))
        if xs:
            ax.plot(xs, ys, "o-", color=c, ms=7, lw=1.8,
                    label=r"$(%d,%d)$" % (na, nb))
    ax.set_xscale("log", base=2)
    ax.set_yscale("log")
    ax.set_xticks(list(NS_FIXED))
    ax.set_xticklabels([str(N) for N in NS_FIXED])
    ax.set_xlabel("grid $N$")
    ax.set_ylabel(r"$|I_D|$  $[\%]$")
    ax.set_title(r"(b) transitions at $Re=400$")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, which="both")

    ax = axes[2]
    for N in NS_FIXED:
        xs, ys = [], []
        for Re in RE_FIXED:
            d = contrast(2, 4, 0.0, Re, N, "D_star")
            p = contrast(2, 4, 0.0, Re, N, "P_star")
            if d is not None and p is not None and abs(d) > 0:
                xs.append(Re)
                ys.append(abs(100.0 * (p / d - 1.0)))
        if xs:
            m, fc = mk[N]
            ax.plot(xs, ys, m + ("-" if N == 384 else "--"), color="#2ca02c",
                    ms=7, lw=1.6, mfc=fc, mec="#2ca02c", label="$N=%d$" % N)
    ax.set_xscale("log")
    ax.set_xlabel("$Re$")
    ax.set_ylabel(r"$\delta_{IP/ID}(2,4)$  $[\%]$")
    ax.set_title(r"(c) estimator discrepancy")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, which="both")

    fig.tight_layout()
    p = os.path.join(OUT, "figS_contrasts_res.pdf")
    fig.savefig(p)
    plt.close(fig)
    return p


def main():
    p6 = fig6()
    p9, Re, I, Rg, g, xc = fig9()
    pS = figS()
    print("[ok] %s" % os.path.relpath(p6, _ROOT))
    print("[ok] %s" % os.path.relpath(p9, _ROOT))
    print("[ok] %s" % os.path.relpath(pS, _ROOT))
    print("[resweep] I_D(2,4) alpha=0: %.4f%% (Re=%d) -> %.4f%% (Re=%d)"
          % (I[0], Re[0], I[-1], Re[-1]))
    print("[resweep] gamma_loc: %.3f near Re=%d -> %.3f near Re=%d ; crossing unity %s"
          % (g[0], Rg[0], g[-1], Rg[-1],
             ("Re~%.0f" % xc) if xc is not None else "n/a"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
