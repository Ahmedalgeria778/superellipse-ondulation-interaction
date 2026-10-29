# -*- coding: utf-8 -*-
"""
solver/plot_fem_xcheck.py
=========================
Independent-solver cross-validation figure (referee request: second solver).

Compares, for Re=400, alpha=0, the geometric contrast I(2,4) of the kinetic
energy (K*), lid power (P*), and viscous dissipation (D*) obtained by the
MRT-LBM and by a conforming Taylor-Hood FEM (P2/P1) solver.

Both solvers are taken at their FINEST grid:
    MRT-LBM  : N = 384
    FEM      : nx = ny = 256

All FEM entries are NATIVE (element-wise strain-rate / lid-power /
kinetic-energy integrals, no projection).  The native FEM data are read
from results/fem/convergence_native.json (assembled from the per-run JSON
records of the mesh-convergence campaign).

Output:
    paper_pof_new/figures/fig10_fem_xcheck.pdf (+ .png 600 dpi)

Usage: python -m solver.plot_fem_xcheck
"""

from __future__ import annotations

import json
import os

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)

import numpy as np                                     # noqa: E402
import matplotlib                                      # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt                        # noqa: E402

SCAL = os.path.join(_ROOT, "results", "scalars")
FEMD = os.path.join(_ROOT, "results", "fem")
OUT = os.path.join(_ROOT, "paper", "figures", "fig10_fem_xcheck.pdf")

# --- FEM native contrast reference (from the mesh-convergence campaign) ---
FEM_NATIVE_NX = 256
LBM_N = 384


def load(path):
    with open(path, "r", encoding="utf-8") as fh:
        return json.load(fh)


def contrast(a2, a4):
    return 100.0 * (a4 - a2) / a2


def main():
    fem_all = load(os.path.join(FEMD, "convergence_native.json"))

    def fem(n, nx, key):
        return fem_all["%d_%d_%d" % (400, n, nx)][key]

    L = {n: load(os.path.join(
        SCAL, "se_n%d_w0.00_Nw1_Re400_N%d_tol1e-08.json" % (n, LBM_N)))
        for n in (2, 4)}

    labels = [r"$K^*$", r"$P^*$", r"$D^*$"]
    lbm = [
        contrast(L[2]["K_star"], L[4]["K_star"]),
        contrast(L[2]["P_star"], L[4]["P_star"]),
        contrast(L[2]["D_star"], L[4]["D_star"]),
    ]
    fem = [
        contrast(fem(2, FEM_NATIVE_NX, "K"), fem(4, FEM_NATIVE_NX, "K")),
        contrast(fem(2, FEM_NATIVE_NX, "P"), fem(4, FEM_NATIVE_NX, "P")),
        contrast(fem(2, FEM_NATIVE_NX, "D"), fem(4, FEM_NATIVE_NX, "D")),
    ]

    x = np.arange(len(labels))
    w = 0.36

    fig, ax = plt.subplots(figsize=(3.4, 2.8))
    c_lbm = "#1f4e79"
    c_fem = "#c0504d"
    b1 = ax.bar(x - w / 2, lbm, w, color=c_lbm, edgecolor="k",
                linewidth=0.5, label=r"MRT-LBM ($N=384$)")
    b2 = ax.bar(x + w / 2, fem, w, color=c_fem, edgecolor="k",
                linewidth=0.5, label=r"FEM (native, $nx=256$)")

    for bars in (b1, b2):
        for r in bars:
            h = r.get_height()
            ax.annotate("%.3f" % h, (r.get_x() + r.get_width() / 2, h),
                        ha="center", va="bottom", fontsize=6.0,
                        xytext=(0, 1), textcoords="offset points")

    ax.set_ylabel(r"contrast $I(2,4)$  [\%]", fontsize=8)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylim(0, max(max(lbm), max(fem)) * 1.28)
    ax.tick_params(axis="y", labelsize=7)
    ax.legend(fontsize=6.0, frameon=False, loc="upper right",
              handlelength=1.2, borderaxespad=0.2)
    ax.grid(axis="y", ls=":", lw=0.4, alpha=0.5)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

    fig.tight_layout(pad=0.3)
    fig.savefig(OUT, dpi=600)
    fig.savefig(OUT.replace(".pdf", ".png"), dpi=600)
    print("wrote", OUT)
    print("LBM (N=%d):" % LBM_N, ["%.4f" % v for v in lbm])
    print("FEM (nx=%d, native):" % FEM_NATIVE_NX, ["%.4f" % v for v in fem])


if __name__ == "__main__":
    main()
