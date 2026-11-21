# -*- coding: utf-8 -*-
"""Final style/consistency pass on the 2026-11-21 V1 revision.

Removes meta-discourse and corrects one numerically inconsistent FEM--LBM gap
(0.008 -> 0.002 percentage points at Re=1000), which the native FEM records
and Section III.C of the manuscript both give as 0.002.

Strictly idempotent: a replacement is applied only if the old text occurs
exactly once; if the new text is already present and the old absent, it is
skipped.  Every new string is chosen so that it does not contain the old one.
"""
from __future__ import annotations

from pathlib import Path

HERE = Path(__file__).resolve().parent


def sub(path, old, new):
    p = HERE / path
    text = p.read_text(encoding="utf-8")
    found = text.count(old)
    if found == 0:
        if new in text:
            return
        raise SystemExit("POLISH %s: old text missing and new text absent: %r"
                         % (path, old[:80]))
    if found != 1:
        raise SystemExit("POLISH %s: expected 1 occurrence, found %d: %r"
                         % (path, found, old[:80]))
    p.write_text(text.replace(old, new), encoding="utf-8")


# --- introduction -----------------------------------------------------------
sub("sections/01_intro.tex",
    r"The central original contribution of this work is the exact $\mathrm{d}A/\mathrm{d}B$ spatial partition of the dissipation contrast. While spatial decompositions are common in mathematical analysis, their application to the energetic accounting of cavity flows---separating the added-region dissipation from the common-region redistribution---provides direct physical insight into why the dissipation response is sub-extensive relative to the area change. This partition demonstrates that the added shoulder region, despite occupying $7.4\%$ of the fluid area at $Re=400$, contributes only $0.2\%$ of the total dissipation, while the common-region redistribution is itself a near-cancellation of opposing lid-band and core-band contributions.",
    r"An exact cell-by-cell partition separates the dissipation contrast into an added-region term $\mathrm{d}B$ and a common-region redistribution $\mathrm{d}A$. At $Re=400$ the added shoulder region occupies $7.4\%$ of the fluid area but contributes $0.2\%$ of the total dissipation, while the common-region redistribution is a near-cancellation of opposing lid-band and core-band contributions. The sub-extensive response follows from the low mean dissipation density in the added region.")

sub("sections/01_intro.tex",
    r"The contribution of the present work is to combine a continuous superelliptic interpolation of the bottom wall with",
    r"The study combines a continuous superelliptic interpolation of the bottom wall with")

sub("sections/01_intro.tex",
    r"This combination provides a complete energetic accounting that was not previously available for corrugated or curved-wall cavity flows.",
    r"This combination gives an energetic accounting for corrugated and curved-wall cavity flows.")

# --- formulation ------------------------------------------------------------
sub("sections/02_formulation.tex",
    r"This is why the analysis is built on $I_D$ and uses $I_P$ only as a cross-check.",
    r"The analysis is therefore built on $I_D$, with $I_P$ retained as a cross-check.")

sub("sections/02_formulation.tex",
    r"Note that the estimator discrepancy $\delta_{IP/ID}$ of Section~\ref{subsec:consistency}",
    r"The estimator discrepancy $\delta_{IP/ID}$ of Section~\ref{subsec:consistency}")

# --- results ----------------------------------------------------------------
sub("sections/04_results.tex",
    r"The global integrals show \emph{that} the dissipation responds weakly to the bottom geometry; we now localize \emph{where} the change occurs.",
    r"The local dissipation density localizes the geometric response.")

sub("sections/04_results.tex",
    r"This algebraic exactness must be distinguished from two statements that it does not imply.",
    r"Two qualifications apply to this identity.")

sub("sections/04_results.tex",
    r"; no curve contains extrapolated, clipped, or missing points.",
    r".")

sub("sections/04_results.tex",
    r"its quantitative details must be distinguished from the primary contrast $I_D(2,4)$ established at $N=384$",
    r"its quantitative details differ from the primary contrast $I_D(2,4)$ established at $N=384$")

sub("sections/04_results.tex",
    r"This is the physical content of the crossover: the super-linear growth at low $Re$ is the signature of the common-region redistribution $\mathrm{d}A$ rising fastest, the flattening near $Re\approx346$ is where $\mathrm{d}A$ begins to saturate, and the sub-linear growth at high $Re$ is where the linear added-region term $\mathrm{d}B$ dominates and $\mathrm{d}A$ declines. The crossover is therefore a direct consequence of the two-term structure of the partition.",
    r"The super-linear growth at low $Re$ follows the common-region redistribution $\mathrm{d}A$ rising fastest; the flattening near $Re\approx346$ marks its saturation; the sub-linear growth at high $Re$ reflects the linear added-region term $\mathrm{d}B$ dominating as $\mathrm{d}A$ declines.")

sub("sections/04_results.tex",
    r"Taken together, the three levels of the analysis are coherent.",
    r"The three levels of the analysis are mutually consistent.")

sub("sections/04_results.tex",
    r"finest-grid FEM--LBM gaps are $0.002$ percentage point at $Re=400$ and $0.008$ at $Re=1000$",
    r"finest-grid FEM--LBM gaps are $0.002$ percentage point at both $Re=400$ and $Re=1000$")

# --- conclusions ------------------------------------------------------------
sub("sections/05_conclusions.tex",
    r"reproduce the primary contrast to $0.002$ percentage points at $Re=400$ and $0.008$ at $Re=1000$ on the finest meshes of each solver",
    r"reproduce the primary contrast to $0.002$ percentage points at $Re=400$ and $Re=1000$ on the finest meshes of each solver")

# --- supplementary ----------------------------------------------------------
sub("supplementary/sections/s3_convergence.tex",
    r"both are grid effects that cancel in the paired contrast, exactly as the main text argues.",
    r"both are grid effects that cancel in the paired contrast.")

print("style polish applied")
