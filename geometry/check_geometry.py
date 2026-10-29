# -*- coding: utf-8 -*-
"""
geometry/check_geometry.py
==========================
Validation geometrique (continuite, aire, limites, pente) + figures.
Repere boite : x in [-a,a], y in [0,H], lid en y=H, epaules du fond en y=b.
"""

from __future__ import annotations

import os
import json

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import sys
_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

from combined_geometry import y_b, dy_b_dx, area_fluid, min_yb  # noqa: E402
from configs import PILOT_4CASES                                  # noqa: E402

FIGDIR = os.path.join(os.path.dirname(_HERE), "figures", "geometry")
os.makedirs(FIGDIR, exist_ok=True)


def check_case(cfg, Nx=8192):
    x = np.linspace(-cfg.a, cfg.a, Nx + 1)
    yb = y_b(x, cfg.n, cfg.alpha, cfg.Nw, cfg.a, cfg.b)
    Af = area_fluid(cfg.n, cfg.alpha, cfg.Nw, cfg.a, cfg.b, cfg.H, Nx=4096)
    dy = dy_b_dx(x, cfg.n, cfg.alpha, cfg.Nw, cfg.a, cfg.b)
    dy = dy[np.isfinite(dy)]
    return {
        "case_id": cfg.case_id,
        "n": cfg.n, "alpha": cfg.alpha, "Nw": cfg.Nw,
        "a": cfg.a, "b": cfg.b, "H": cfg.H,
        "A_f": Af,
        "y_b_min": float(np.min(yb)),
        "y_b(x=-a)": float(yb[0]),
        "y_b(x=+a)": float(yb[-1]),
        "slope_max_abs": float(np.max(np.abs(dy))) if dy.size else 0.0,
        "depth_center": float(cfg.H - np.min(yb)),
        "depth_wall": float(cfg.H - cfg.b),
    }


def plot_case(cfg, Nx=512):
    x = np.linspace(-cfg.a, cfg.a, Nx + 1)
    yb = y_b(x, cfg.n, cfg.alpha, cfg.Nw, cfg.a, cfg.b)
    fig, ax = plt.subplots(figsize=(6, 3.5))
    ax.fill_between(x, yb, cfg.H, facecolor="#f0f7ff", alpha=0.6, label="fluide")
    ax.plot(x, yb, lw=1.6, color="C0", label="y_b(x)")
    ax.axhline(cfg.H, ls="-", lw=1.2, color="k", label="lid (mobile)")
    ax.axvline(-cfg.a, ls="--", lw=0.7, color="gray")
    ax.axvline(cfg.a, ls="--", lw=0.7, color="gray")
    ax.set_xlim(-cfg.a - 0.01, cfg.a + 0.01)
    ax.set_ylim(-0.02, cfg.H + 0.05)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(cfg.case_id)
    ax.grid(True, alpha=0.25)
    ax.legend(loc="upper center", ncol=3, fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, f"{cfg.case_id}_profile.png"), dpi=150)
    plt.close(fig)


def main():
    res = [check_case(c) for c in PILOT_4CASES]
    for c in PILOT_4CASES:
        plot_case(c)
    for r in res:
        print(r)
    out = os.path.join(FIGDIR, "geometry_checks.json")
    with open(out, "w", encoding="utf-8") as fh:
        json.dump(res, fh, indent=2, ensure_ascii=False)
    print(f"[ok] ecrit {out}")


if __name__ == "__main__":
    main()
