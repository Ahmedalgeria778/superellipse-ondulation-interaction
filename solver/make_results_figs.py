# -*- coding: utf-8 -*-
"""
solver/make_results_figs.py
===========================
Regenerates the four results-section figures that depend on the primary
campaign scalars (results/scalars/se_*_tol1e-08.json) and on the saved
lattice fields (results/fields/se_*_Re400_N384_tol1e-08.npz):

  fig4_energetic.pdf        P*, D*, K* vs n at N=384 (alpha=0 solid, 0.1 dashed),
                            one panel per Re (100/400/1000)
  fig5_contrasts.pdf        (a) Re=400 three successive transitions (2,4)(4,6)(6,8):
                            I_D, I_P, I_area on a symlog ordinate
                            (b) primary (2,4) at Re=100/400/1000
  fig7_dissipation_map.pdf  local dissipation density 2 S:S at Re=400, alpha=0, N=384:
                            (left) n=2, (middle) n=4, (right) difference on the n=2 region
  fig8_reference.pdf        primary (2,4), alpha=0, N=384: I_D, I_D/A (intrinsic),
                            and the extensive area-only reference I_area, symlog

All quantities use the exact solver convention reproduced from
diagnostics.compute_diagnostics_fields, verified to machine precision against
the scalar JSONs. Outputs go straight into paper_pof_new/figures/.

Usage: python -m solver.make_results_figs
"""

from __future__ import annotations

import json
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.ticker import NullFormatter

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
for _p in (_HERE, _ROOT, os.path.join(_ROOT, "provenance", "Lid_driven_ondulation")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import combined_geometry as CG                        # noqa: E402
from geometry import lid_profile_regularized          # noqa: E402

SCAL = os.path.join(_ROOT, "results", "scalars")
FIELDS = os.path.join(_ROOT, "results", "fields")
FIGS = os.path.join(_ROOT, "paper", "figures")

RE_ALL = (100, 400, 1000)
NS_ALL = (2, 4, 6, 8)
ALPHAS = (0.0, 0.10)
EPS_LID = 0.02


# ------------------------------------------------------------------ loaders
def rec(n, a, Re, N):
    fn = "se_n%g_w%.2f_Nw1_Re%d_N%d_tol1e-08.json" % (n, a, Re, N)
    p = os.path.join(SCAL, fn)
    if not os.path.exists(p):
        return None
    with open(p, encoding="utf-8") as fh:
        return json.load(fh)


def contrast(n_a, n_b, a, Re, N, key):
    ra, rb = rec(n_a, a, Re, N), rec(n_b, a, Re, N)
    if ra is None or rb is None:
        return None
    return 100.0 * (rb[key] - ra[key]) / ra[key]


def area_ref(n_a, n_b, a, Re, N):
    """Extensive area-only reference I_{D,area} = (A_b - A_a)/A_a (%)."""
    ra, rb = rec(n_a, a, Re, N), rec(n_b, a, Re, N)
    if ra is None or rb is None:
        return None
    return 100.0 * (rb["fluid_area"] - ra["fluid_area"]) / ra["fluid_area"]


def intrinsic(n_a, n_b, a, Re, N):
    """Per-unit-area dissipation density contrast I_{D/A} (%)."""
    ra, rb = rec(n_a, a, Re, N), rec(n_b, a, Re, N)
    if ra is None or rb is None:
        return None
    da = ra["D_star"] / ra["fluid_area"]
    db = rb["D_star"] / rb["fluid_area"]
    return 100.0 * (db - da) / da


# ------------------------------------------------------------- fig 4
def fig4():
    fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.8))
    mk = {"P_star": ("o", "P*"), "D_star": ("s", "D*"), "K_star": ("^", "K*")}
    for ax, Re in zip(axes, RE_ALL):
        present = set()
        for a, ls in zip(ALPHAS, ("-", "--")):
            for key, (m, lab) in mk.items():
                xs, ys = [], []
                for n in NS_ALL:
                    r = rec(n, a, Re, 384)
                    if r is not None:
                        xs.append(n)
                        ys.append(r[key])
                        present.add(n)
                if not xs:
                    continue
                ax.plot(xs, ys, m + ls, ms=5, lw=1.4,
                        label=("%s" % lab) if a == 0.0 else None)
        ax.set_yscale("log")
        ax.set_xticks(sorted(present))
        ax.set_xlabel("exponent $n$")
        ax.set_ylabel(r"$P^*,D^*,K^*$")
        ax.set_title(r"$Re=%d$" % Re)
        ax.grid(alpha=0.3, which="both")
        if len(present) < len(NS_ALL):
            ax.annotate("only $n=2,4$ available",
                        xy=(0.5, 0.03), xycoords="axes fraction",
                        ha="center", va="bottom", fontsize=7, style="italic",
                        color="0.35")
    # combined legend: markers + line styles
    h = [plt.Line2D([], [], marker=m, ls="none", color="k", ms=5, label=lab)
         for m, lab in (("o", r"$P^*$"), ("s", r"$D^*$"), ("^", r"$K^*$"))]
    h += [plt.Line2D([], [], ls="-", color="k", lw=1.4, label=r"$\alpha=0$"),
          plt.Line2D([], [], ls="--", color="k", lw=1.4, label=r"$\alpha=0.1$")]
    axes[0].legend(handles=h, fontsize=8, ncol=2, frameon=False)
    fig.tight_layout()
    p = os.path.join(FIGS, "fig4_energetic.pdf")
    fig.savefig(p)
    fig.savefig(p.replace(".pdf", ".png"), dpi=600)
    plt.close(fig)
    return p


# ------------------------------------------------------------- fig 5
def _symlog(ax):
    ax.set_yscale("symlog", linthresh=1e-4)
    ax.yaxis.set_minor_formatter(NullFormatter())


def fig5():
    fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.6))

    # panel (a): three transitions at Re=400
    ax = axes[0]
    trans = [(2, 4), (4, 6), (6, 8)]
    xpos = np.arange(len(trans))
    idr = [contrast(a, b, 0.0, 400, 384, "D_star") for a, b in trans]
    ipr = [contrast(a, b, 0.0, 400, 384, "P_star") for a, b in trans]
    iar = [area_ref(a, b, 0.0, 400, 384) for a, b in trans]
    ax.plot(xpos, idr, "o-", color="#1f77b4", ms=6, lw=1.5, label=r"$I_D$")
    ax.plot(xpos, ipr, "s--", color="#ff7f0e", ms=6, lw=1.5, label=r"$I_P$")
    ax.plot(xpos, iar, "^:", color="#2ca02c", ms=6, lw=1.5,
            label=r"$I_{D,\mathrm{area}}$")
    ax.set_xticks(xpos)
    ax.set_xticklabels([r"$(%d,%d)$" % t for t in trans])
    _symlog(ax)
    ax.set_xlabel("transition")
    ax.set_ylabel(r"contrast  $[\%]$")
    ax.set_title(r"(a) $Re=400$")
    ax.legend(fontsize=8, frameon=False)
    ax.grid(alpha=0.3, which="both")

    # panel (b): primary (2,4) vs Re
    ax = axes[1]
    xs = list(RE_ALL)
    idv = [contrast(2, 4, 0.0, Re, 384, "D_star") for Re in xs]
    ipv = [contrast(2, 4, 0.0, Re, 384, "P_star") for Re in xs]
    iav = [area_ref(2, 4, 0.0, Re, 384) for Re in xs]
    ax.plot(xs, idv, "o-", color="#1f77b4", ms=6, lw=1.5, label=r"$I_D$")
    ax.plot(xs, ipv, "s--", color="#ff7f0e", ms=6, lw=1.5, label=r"$I_P$")
    ax.plot(xs, iav, "^:", color="#2ca02c", ms=6, lw=1.5,
            label=r"$I_{D,\mathrm{area}}$")
    ax.set_xscale("log")
    ax.set_xticks(list(RE_ALL))
    ax.set_xticklabels([str(r) for r in RE_ALL])
    _symlog(ax)
    ax.set_xlabel("$Re$")
    ax.set_ylabel(r"contrast  $[\%]$")
    ax.set_title(r"(b) primary $(2,4)$")
    ax.legend(fontsize=8, frameon=False)
    ax.grid(alpha=0.3, which="both")

    fig.tight_layout()
    p = os.path.join(FIGS, "fig5_contrasts.pdf")
    fig.savefig(p)
    fig.savefig(p.replace(".pdf", ".png"), dpi=600)
    plt.close(fig)
    return p


# ------------------------------------------------------------- fig 8
def fig8():
    fig, ax = plt.subplots(figsize=(4.2, 3.6))
    xs = list(RE_ALL)
    idv = [contrast(2, 4, 0.0, Re, 384, "D_star") for Re in xs]
    ida = [intrinsic(2, 4, 0.0, Re, 384) for Re in xs]
    iav = [area_ref(2, 4, 0.0, Re, 384) for Re in xs]
    ax.plot(xs, idv, "o-", color="#1f77b4", ms=7, lw=1.6, label=r"$I_D$")
    ax.plot(xs, ida, "s--", color="#d62728", ms=7, lw=1.6,
            label=r"$I_{D/A}$")
    ax.plot(xs, iav, "^:", color="#2ca02c", ms=7, lw=1.6,
            label=r"$I_{D,\mathrm{area}}$")
    for x, y in zip(xs, ida):
        ax.annotate("%.1f" % y, (x, y), textcoords="offset points",
                    xytext=(0, -12), ha="center", fontsize=7)
    ax.set_xscale("log")
    ax.set_xticks(list(RE_ALL))
    ax.set_xticklabels([str(r) for r in RE_ALL])
    _symlog(ax)
    ax.set_xlabel("$Re$")
    ax.set_ylabel(r"contrast  $[\%]$")
    ax.legend(fontsize=8, frameon=False)
    ax.grid(alpha=0.3, which="both")
    fig.tight_layout()
    p = os.path.join(FIGS, "fig8_reference.pdf")
    fig.savefig(p)
    fig.savefig(p.replace(".pdf", ".png"), dpi=600)
    plt.close(fig)
    return p


# ------------------------------------------------------------- fig 7
def _gradients(U, V, geom, Uln):
    N, dx = geom["N"], geom["dx"]
    fl = geom["fluid"].astype(bool)
    dudx = np.gradient(U, dx, axis=1)
    dudy = np.gradient(U, dx, axis=0)
    dvdx = np.gradient(V, dx, axis=1)
    dvdy = np.gradient(V, dx, axis=0)
    for j in range(N):
        if fl[N - 1, j]:
            dudy[N - 1, j] = (Uln[j] - U[N - 1, j]) / (0.5 * dx)
            dvdy[N - 1, j] = (0.0 - V[N - 1, j]) / (0.5 * dx)
    for i in range(N):
        if fl[i, 0]:
            dudx[i, 0] = (0.0 - U[i, 0]) / (0.5 * dx)
            dvdx[i, 0] = (0.0 - V[i, 0]) / (0.5 * dx)
        if fl[i, N - 1]:
            dudx[i, N - 1] = (0.0 - U[i, N - 1]) / (0.5 * dx)
            dvdx[i, N - 1] = (0.0 - V[i, N - 1]) / (0.5 * dx)
    yb = geom["yb"]
    for j in range(N):
        col = np.nonzero(fl[:, j])[0]
        if col.size == 0:
            continue
        i0 = col[0]
        dist = (i0 + 0.5) * dx - yb[0, j]
        if dist > 1e-12:
            dudy[i0, j] = (0.0 - U[i0, j]) / dist
            dvdy[i0, j] = (0.0 - V[i0, j]) / dist
    return dudx, dudy, dvdx, dvdy, fl


def _dissipation_density(n, alpha, Re):
    fn = "se_n%g_w%.2f_Nw1_Re%d_N384_tol1e-08.npz" % (n, alpha, Re)
    d = np.load(os.path.join(FIELDS, fn), allow_pickle=True)
    N = int(d["N"])
    geom = CG.build_combined_geometry(N, float(d["n"]), float(d["alpha"]),
                                      int(d["Nw"]))
    fl = geom["fluid"].astype(bool)
    U = np.where(fl, np.array(d["ux"]), 0.0)
    V = np.where(fl, np.array(d["uy"]), 0.0)
    Uln = lid_profile_regularized(geom["xc"], EPS_LID)
    dudx, dudy, dvdx, dvdy, fl = _gradients(U, V, geom, Uln)
    S11, S22, S12 = dudx, dvdy, 0.5 * (dudy + dvdx)
    dens = 2.0 * (S11 ** 2 + S22 ** 2 + 2.0 * S12 ** 2)
    return dens, fl, geom


def fig7():
    Re = 400
    d2, f2, g2 = _dissipation_density(2, 0.0, Re)
    d4, f4, g4 = _dissipation_density(4, 0.0, Re)
    # common region = n=2 fluid mask; both grids are N=384 aligned
    common = f2 & f4
    d2m = np.where(common, d2, np.nan)
    d4m = np.where(common, d4, np.nan)
    diff = np.where(common, d4 - d2, np.nan)

    xc, yc = g2["xc"], g2["yc"]
    ext = [0, 1, 0, 1]
    fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.7))
    vmax = np.nanpercentile(d2m, 99.5)
    for ax, field, ttl, cmap, kw in (
            (axes[0], d2m, r"$n=2$", "inferno", {}),
            (axes[1], d4m, r"$n=4$", "inferno", {}),
            (axes[2], diff, r"$(n{=}4)-(n{=}2)$ (on $n{=}2$ region)", "RdBu_r",
             {"diverging": True})):
        if kw.get("diverging"):
            lim = np.nanpercentile(np.abs(diff), 99.0)
            norm = TwoSlopeNorm(vcenter=0.0, vmin=-lim, vmax=+lim)
            im = ax.imshow(field, origin="lower", extent=ext, cmap=cmap,
                           norm=norm, interpolation="nearest", aspect="equal")
        else:
            im = ax.imshow(field, origin="lower", extent=ext, cmap=cmap,
                           vmin=0.0, vmax=vmax, interpolation="nearest",
                           aspect="equal")
        # wall outline
        ax.plot(xc, g2["yb"][0, :], "k-", lw=1.2)
        ax.set_title(ttl, fontsize=9)
        ax.set_xlabel("$x$")
        cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, shrink=0.9)
        cb.ax.tick_params(labelsize=7)
        if kw.get("diverging"):
            cb.set_label(r"$\Delta\,2S{:}S$", fontsize=8)
        else:
            cb.set_label(r"$2S{:}S$", fontsize=8)
    axes[0].set_ylabel("$y$")
    fig.tight_layout()
    p = os.path.join(FIGS, "fig7_dissipation_map.pdf")
    fig.savefig(p)
    fig.savefig(p.replace(".pdf", ".png"), dpi=600)
    plt.close(fig)
    # report the added/common partition for the record
    added = f4 & ~f2
    add_d = float(np.sum(d4[added]) * g4["dx"] ** 2)
    tot4 = float(np.sum(d4[f4]) * g4["dx"] ** 2)
    print("[fig7] added region (n=4 only): area share %.2f%% of A_f(n=4), "
          "D share %.4f%% of D*(n=4)"
          % (100.0 * added.sum() / f4.sum(), 100.0 * add_d / tot4))
    return p


# ------------------------------------------------------------- fig S (SM)
def figS_deltaD():
    """SM version: local dissipation contrast Delta(n=4 - n=2) on the n=2
    region, one panel per Re (100/400/1000), vector output."""
    fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.7))
    for ax, Re in zip(axes, RE_ALL):
        d2, f2, g2 = _dissipation_density(2, 0.0, Re)
        d4, f4, _ = _dissipation_density(4, 0.0, Re)
        common = f2 & f4
        diff = np.where(common, d4 - d2, np.nan)
        lim = np.nanpercentile(np.abs(diff), 99.0)
        norm = TwoSlopeNorm(vcenter=0.0, vmin=-lim, vmax=+lim)
        im = ax.imshow(diff, origin="lower", extent=[0, 1, 0, 1], cmap="RdBu_r",
                       norm=norm, interpolation="nearest", aspect="equal")
        ax.plot(g2["xc"], g2["yb"][0, :], "k-", lw=1.2)
        ax.set_title(r"$Re=%d$" % Re, fontsize=9)
        ax.set_xlabel("$x$")
        cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, shrink=0.9)
        cb.ax.tick_params(labelsize=7)
        cb.set_label(r"$\Delta\,2S{:}S$", fontsize=8)
    axes[0].set_ylabel("$y$")
    fig.tight_layout()
    out = os.path.join(_ROOT, "paper", "supplementary", "figures",
                       "figS_deltaD_maps.pdf")
    fig.savefig(out)
    fig.savefig(out.replace(".pdf", ".png"), dpi=600)
    plt.close(fig)
    return out


def main():
    print("[ok]", os.path.relpath(fig4(), _ROOT))
    print("[ok]", os.path.relpath(fig5(), _ROOT))
    print("[ok]", os.path.relpath(fig7(), _ROOT))
    print("[ok]", os.path.relpath(fig8(), _ROOT))
    print("[ok]", os.path.relpath(figS_deltaD(), _ROOT))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
