# -*- coding: utf-8 -*-
"""
solver/referee_mech_analysis.py
===============================
Analyses demandees par le rapport du referee (points 1 et 2) :

  (A) part de la couche-couvercle (bande superieure, 10 % de la hauteur)
      dans D* (absolu) pour n=2, N=384, aux 3 Reynolds  -> corrige
      l'affirmation "74 % qui monte avec Re".
  (B) partition exacte dA / dB du contraste (2,4) a alpha=0.1,
      aux 3 Reynolds (auparavant seulement alpha=0).
  (C) densite de dissipation dans la zone AJOUTEE (region dB) rapportee
      a la densite moyenne de D*, i.e. le "degre de mort" de la zone ajoutee.

Sortie : results/referee/mech_analysis.json + affichage.

Usage : python -m solver.referee_mech_analysis
"""

from __future__ import annotations

import json
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
for _p in (_HERE, _ROOT):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import numpy as np                                     # noqa: E402

from solver.combined_geometry import build_combined_geometry  # noqa: E402

FIELDS = os.path.join(_ROOT, "results", "fields")
OUTDIR = os.path.join(_ROOT, "results", "referee")
os.makedirs(OUTDIR, exist_ok=True)

N = 384
A = 0.5
B = 0.5


def load(n, alpha, Re):
    fn = os.path.join(FIELDS,
                      "se_n%d_w%.2f_Nw1_Re%d_N%d_tol1e-08.npz" % (n, alpha, Re, N))
    d = np.load(fn)
    return d


def gradients(U, V, geom, U_lid_norm):
    Nn = geom["N"]
    d = geom["dx"]
    fluid = geom["fluid"].astype(bool)
    dudx = np.gradient(U, d, axis=1)
    dudy = np.gradient(U, d, axis=0)
    dvdx = np.gradient(V, d, axis=1)
    dvdy = np.gradient(V, d, axis=0)
    for j in range(Nn):
        if fluid[Nn - 1, j]:
            dudy[Nn - 1, j] = (U_lid_norm[j] - U[Nn - 1, j]) / (0.5 * d)
            dvdy[Nn - 1, j] = (0.0 - V[Nn - 1, j]) / (0.5 * d)
    for i in range(Nn):
        if fluid[i, 0]:
            dudx[i, 0] = (0.0 - U[i, 0]) / (0.5 * d)
            dvdx[i, 0] = (0.0 - V[i, 0]) / (0.5 * d)
        if fluid[i, Nn - 1]:
            dudx[i, Nn - 1] = (0.0 - U[i, Nn - 1]) / (0.5 * d)
            dvdx[i, Nn - 1] = (0.0 - V[i, Nn - 1]) / (0.5 * d)
    yb = geom["yb"]
    for j in range(Nn):
        col = np.nonzero(fluid[:, j])[0]
        if col.size == 0:
            continue
        i0 = col[0]
        dist = (i0 + 0.5) * d - yb[0, j]
        if dist > 1e-12:
            dudy[i0, j] = (0.0 - U[i0, j]) / dist
            dvdy[i0, j] = (0.0 - V[i0, j]) / dist
    return dudx, dudy, dvdx, dvdy, fluid


def lid_profile_regularized(x, eps=0.02):
    x = np.asarray(x, dtype=np.float64)
    g = np.ones_like(x)
    lo = x < eps
    hi = x > 1.0 - eps
    tl = x[lo] / eps
    g[lo] = 3.0 * tl ** 2 - 2.0 * tl ** 3
    th = (1.0 - x[hi]) / eps
    g[hi] = 3.0 * th ** 2 - 2.0 * th ** 3
    return g


def diss_density(d, geom):
    """Return (phi, fluid) where phi = 2(S11^2+S22^2+2S12^2), D* = sum(phi)*dx*dy."""
    ux = d["ux"]
    uy = d["uy"]
    U_lid = lid_profile_regularized(geom["xc"], 0.02)
    U = np.where(geom["fluid"].astype(bool), ux, 0.0)
    V = np.where(geom["fluid"].astype(bool), uy, 0.0)
    dudx, dudy, dvdx, dvdy, fluid = gradients(U, V, geom, U_lid)
    S11, S22, S12 = dudx, dvdy, 0.5 * (dudy + dvdx)
    phi = 2.0 * (S11 ** 2 + S22 ** 2 + 2.0 * S12 ** 2)
    return phi, fluid


def lid_band_share(Re, frac=0.10):
    d = np.load(os.path.join(
        FIELDS, "se_n2_w0.00_Nw1_Re%d_N%d_tol1e-08.npz" % (Re, N)))
    geom = build_combined_geometry(N, 2, 0.0, 1, a=A, b=B)
    phi, fluid = diss_density(d, geom)
    dx = geom["dx"]
    Y = geom["Y"]
    top = (Y >= 1.0 - frac) & fluid
    Dtot = float(np.sum(phi[fluid]) * dx * dx)
    Dtop = float(np.sum(phi[top]) * dx * dx)
    return Dtot, Dtop, Dtop / Dtot


def partition(n2, n4, alpha, Re):
    """Exact dA/dB partition: dB = region added by flattening (n4 fluid and not n2),
    dA = common region (both fluid).  Returns sums of the local difference
    phi(n4)-phi(n2) over each region."""
    d2 = np.load(os.path.join(
        FIELDS, "se_n2_w%.2f_Nw1_Re%d_N%d_tol1e-08.npz" % (alpha, Re, N)))
    d4 = np.load(os.path.join(
        FIELDS, "se_n4_w%.2f_Nw1_Re%d_N%d_tol1e-08.npz" % (alpha, Re, N)))
    g2 = build_combined_geometry(N, n2, alpha, 1, a=A, b=B)
    g4 = build_combined_geometry(N, n4, alpha, 1, a=A, b=B)
    phi2, f2 = diss_density(d2, g2)
    phi4, f4 = diss_density(d4, g4)
    dx = g2["dx"]
    common = f2 & f4
    added = f4 & (~f2)
    dA = float(np.sum((phi4 - phi2)[common]) * dx * dx)
    dB = float(np.sum(phi4[added]) * dx * dx)
    A_common = float(common.sum()) * dx * dx
    A_added = float(added.sum()) * dx * dx
    D4 = float(np.sum(phi4[f4]) * dx * dx)
    return dict(dA=dA, dB=dB, net=dA + dB, D4=D4,
                A_common=A_common, A_added=A_added,
                dens_added=dB / A_added, dens_mean=D4 / (g4["fluid_area"]))


def main():
    out = {"lid_band": {}, "partition_alpha": {}}
    print("=== (A) lid top-10%% band share of D* (n=2, N=%d) ===" % N)
    for Re in (100, 400, 1000):
        Dtot, Dtop, share = lid_band_share(Re)
        out["lid_band"][Re] = dict(D=Dtot, Dtop=Dtop, share_pct=100 * share)
        print("  Re=%4d  D*=%.4f  Dtop=%.4f  share=%.2f%%" % (Re, Dtot, Dtop, 100 * share))

    for alpha in (0.0, 0.1):
        out["partition_alpha"][alpha] = {}
        print("\n=== (B/C) dA/dB partition (2,4), alpha=%.1f, N=%d ===" % (alpha, N))
        for Re in (100, 400, 1000):
            r = partition(2, 4, alpha, Re)
            out["partition_alpha"][alpha][Re] = r
            dens_ratio = r["dens_added"] / r["dens_mean"]
            print("  Re=%4d  dA=%+.6f dB=%+.6f net=%+.6f  dB/D4=%.3f%%  "
                  "A_added=%.4f (%.2f%%)  dens_added/dens_mean=%.4f"
                  % (Re, r["dA"], r["dB"], r["net"],
                     100 * r["dB"] / r["D4"], r["A_added"],
                     100 * r["A_added"] / (r["A_common"] + r["A_added"]),
                     dens_ratio))

    with open(os.path.join(OUTDIR, "mech_analysis.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print("\nwrote", os.path.join(OUTDIR, "mech_analysis.json"))


if __name__ == "__main__":
    main()
