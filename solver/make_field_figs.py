# -*- coding: utf-8 -*-
"""
solver/make_field_figs.py
=========================
Regenerates the streamline maps coloured by velocity magnitude, in the same
visual convention as the main-text Fig. 3 (2 rows x 4 columns at Re=400):
streamplot of the converged field, line colour = |u| on a common scale,
panel title carries (n, alpha, fluid area).

  fig3_fields.pdf              Re=400, alpha=0 and 0.1, n=2..8        (main)
  figS_fields_alpha0.pdf       rows Re=100/400/1000, n=2..8, alpha=0   (SM)
  figS_fields_alpha01.pdf      same for alpha=0.1                      (SM)

Usage: python -m solver.make_field_figs [--skip-main]
"""

from __future__ import annotations

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
FLD = os.path.join(_ROOT, "results", "fields")
MAINF = os.path.join(_ROOT, "paper", "figures")
SUPPF = os.path.join(_ROOT, "paper", "supplementary", "figures")

RES = (100, 400, 1000)
EXPO = (2, 4, 6, 8)


def path(n, alpha, Re):
    return os.path.join(FLD, "se_n%d_w%.2f_Nw1_Re%d_N384_tol1e-08.npz" % (n, alpha, Re))


def load(n, alpha, Re):
    p = path(n, alpha, Re)
    if not os.path.exists(p):
        return None
    d = np.load(p)
    speed = np.sqrt(d["ux"] ** 2 + d["uy"] ** 2)
    fluid = d["fluid"].astype(bool)
    return dict(xc=d["xc"], yc=d["yc"], ux=np.where(fluid, d["ux"], 0.0),
                uy=np.where(fluid, d["uy"], 0.0), speed=np.where(fluid, speed, np.nan),
                af=float(fluid.mean()))


def panel(ax, d, vmax, title):
    ax.set_aspect("equal")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xticks([0, 0.5, 1.0])
    ax.set_yticks([0, 0.5, 1.0])
    ax.tick_params(labelsize=6)
    ax.set_title(title, fontsize=7)
    if d is None:
        ax.text(0.5, 0.5, "not simulated", ha="center", va="center", fontsize=7)
        return
    ax.streamplot(d["xc"], d["yc"], d["ux"], d["uy"],
                  color=np.nan_to_num(np.clip(d["speed"] / vmax, 0.0, 1.0)),
                  cmap="viridis", density=1.5, linewidth=0.8)


def add_colorbar(fig, sm, axes):
    """Vertical bar on a dedicated axes, spanning exactly the panel block."""
    box = axes[0, 0].get_position().union(
        [a.get_position() for a in axes.ravel()])
    cax = fig.add_axes([box.x1 + 0.012, box.y0, 0.013, box.height])
    return fig.colorbar(sm, cax=cax)


def grid(alpha, out, title):
    data = {(Re, n): load(n, alpha, Re) for Re in RES for n in EXPO}
    finite = [np.nanmax(v["speed"]) for v in data.values() if v is not None]
    vmax = float(np.nanmax(finite))
    fig, axes = plt.subplots(len(RES), len(EXPO), figsize=(10, 7.4))
    for i, Re in enumerate(RES):
        for j, n in enumerate(EXPO):
            d = data[(Re, n)]
            lab = "$n=%d$, $Re=%d$, $A_f=%.3f$" % (
                n, Re, d["af"]) if d is not None else "$n=%d$, $Re=%d$" % (n, Re)
            panel(axes[i, j], d, vmax, lab)
            axes[i, j].set_ylabel("$Re=%d$" % Re if j == 0 else "")
    sm = plt.cm.ScalarMappable(cmap="viridis", norm=plt.Normalize(vmin=0, vmax=vmax))
    fig.suptitle(title, fontsize=9)
    fig.tight_layout(rect=(0, 0, 0.955, 0.96))
    cbar = add_colorbar(fig, sm, axes)
    cbar.set_label(r"$\|\mathbf{u}\|$", fontsize=8)
    cbar.ax.tick_params(labelsize=6)
    fig.savefig(out)
    plt.close(fig)
    print("[ok]", os.path.relpath(out, _ROOT), " vmax=%.3f" % vmax)


def main_grid(out):
    data = {(a, n): load(n, a, 400) for a in (0.0, 0.1) for n in EXPO}
    vmax = float(np.nanmax([np.nanmax(v["speed"]) for v in data.values()]))
    fig, axes = plt.subplots(2, len(EXPO), figsize=(10, 5.0))
    for i, a in enumerate((0.0, 0.1)):
        for j, n in enumerate(EXPO):
            panel(axes[i, j], data[(a, n)], vmax,
                  "$n=%d$, $\\alpha$=%.1f" % (n, a))
    sm = plt.cm.ScalarMappable(cmap="viridis", norm=plt.Normalize(vmin=0, vmax=vmax))
    fig.suptitle("$Re=400$, $N=384$: streamlines colored by $|\\mathbf{u}|$", fontsize=9)
    fig.tight_layout(rect=(0, 0, 0.955, 0.95))
    cbar = add_colorbar(fig, sm, axes)
    cbar.set_label(r"$\|\mathbf{u}\|$", fontsize=8)
    cbar.ax.tick_params(labelsize=6)
    fig.savefig(out)
    plt.close(fig)
    print("[ok]", os.path.relpath(out, _ROOT))


if __name__ == "__main__":
    import sys
    os.makedirs(SUPPF, exist_ok=True)
    grid(0.0, os.path.join(SUPPF, "figS_fields_alpha0.pdf"),
         "$N=384$, $\\alpha=0$: streamlines colored by $\\|\\mathbf{u}\\|$")
    grid(0.1, os.path.join(SUPPF, "figS_fields_alpha01.pdf"),
         "$N=384$, $\\alpha=0.1$: streamlines colored by $\\|\\mathbf{u}\\|$")
    if "--with-main" in sys.argv:
        main_grid(os.path.join(MAINF, "fig3_fields.pdf"))
