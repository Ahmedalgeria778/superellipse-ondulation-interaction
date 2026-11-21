# -*- coding: utf-8 -*-
"""Regenerate the SM dA/dB-vs-resolution figure for the 2026-11-21 V1 revision.

Adapted from make_dadb_fig.py; reads the extended partition recomputed by
extract_revision_dAdB.py (which now includes N=768 at Re=400 and Re=1000).
Only this project's own partition data is plotted.
"""
import json, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
data = json.load(open(os.path.join(HERE, "revision_dAdB.json"), encoding="utf-8"))
OUT = os.path.join(HERE, "supplementary", "figures", "figS_dadb_res.pdf")

by_re = {}
for r in data:
    by_re.setdefault(r["Re"], []).append(r)

fig, axes = plt.subplots(1, 3, figsize=(10.0, 3.0), sharey=False)
S = 1e4
for ax, Re in zip(axes, (100, 400, 1000)):
    rs = sorted(by_re[Re], key=lambda d: d["N"])
    Ns = [r["N"] for r in rs]
    dA = [r["dA_frac"] * S for r in rs]
    dB = [r["dB_frac"] * S for r in rs]
    DD = [r["Delta_frac"] * S for r in rs]
    ax.axhline(0, color="grey", lw=0.7)
    ax.plot(Ns, dB, "o-", color="#1b6ca8", label=r"$\mathrm{d}B$")
    ax.plot(Ns, dA, "s-", color="#c0392b", label=r"$\mathrm{d}A$")
    ax.plot(Ns, DD, "^-", color="#2c2c2c", lw=1.6, label=r"$\Delta D=\mathrm{d}A+\mathrm{d}B$")
    ax.set_xscale("log", base=2)
    ax.set_xticks(Ns)
    ax.set_xticklabels([str(n) for n in Ns])
    ax.set_xlabel("$N$")
    ax.set_title("$Re=%d$" % Re)
    ax.grid(alpha=0.25, ls=":")
    lo = min(dA + [0]) * 1.15
    hi = max(dB + DD) * 1.25
    ax.set_ylim(lo, hi)
axes[0].set_ylabel(r"$10^{-4}\,D^*(2,0)$")
axes[0].legend(frameon=False, fontsize=8, loc="center right")
fig.tight_layout()
fig.savefig(OUT)
print("wrote", OUT)
