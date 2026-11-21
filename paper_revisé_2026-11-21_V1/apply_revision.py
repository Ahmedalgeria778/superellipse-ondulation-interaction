# -*- coding: utf-8 -*-
"""Apply the 2026-11-21 V1 revision to the manuscript copy.

Every edit is checked for uniqueness, so the script fails rather than silently
changing the wrong passage.  Numbers come from revision_data.py (all scalar
records) and extract_revision_dAdB.py (six-band partition at N=768).
"""
from __future__ import annotations

from pathlib import Path

HERE = Path(__file__).resolve().parent
EDITED = []


def load(path):
    p = HERE / path
    return p, p.read_text(encoding="utf-8")


def save(path, text):
    (HERE / path).write_text(text, encoding="utf-8")
    EDITED.append(path)


def sub(path, old, new, count=1):
    p, text = load(path)
    found = text.count(old)
    if found == 0 and text.count(new) >= count:
        return  # already applied by an earlier run
    if found != count:
        raise SystemExit("SUB %s: expected %d occurrence(s), found %d: %r"
                         % (path, count, found, old[:90]))
    save(path, text.replace(old, new))


def line(path, prefix, new):
    p, text = load(path)
    lines = text.split("\n")
    hits = [i for i, l in enumerate(lines) if l.strip().startswith(prefix)]
    if not hits and new in lines:
        return  # already applied by an earlier run
    if len(hits) != 1:
        raise SystemExit("LINE %s: prefix %r matched %d lines" % (path, prefix, len(hits)))
    lines[hits[0]] = new
    save(path, "\n".join(lines))


def insert_after(path, prefix, new_lines):
    p, text = load(path)
    lines = text.split("\n")
    hits = [i for i, l in enumerate(lines) if l.strip().startswith(prefix)]
    if not hits and all(l in lines for l in new_lines):
        return  # already applied by an earlier run
    if len(hits) != 1:
        raise SystemExit("INSERT %s: prefix %r matched %d lines" % (path, prefix, len(hits)))
    i = hits[0] + 1
    lines[i:i] = new_lines
    save(path, "\n".join(lines))


# --------------------------------------------------------------- main.tex
sub("main.tex",
    r"Simulations use grids $N=128$--$384$, with cases at $N=512$ and a fourth level at $N=768$ for the primary transition at $Re=100$.",
    r"Simulations use grids $N=128$--$384$, extended to $N=512$ at $Re=400$ and $1000$ and to $N=768$ for the primary transition ($n=2$ and $4$, $\alpha=0$) at all three Reynolds numbers.")

sub("main.tex",
    r"while total integrated viscous dissipation rises by $0.0218\%$, $0.3249\%$, and $0.5909\%$, respectively. An area-only extensive model",
    r"while total integrated viscous dissipation rises by $0.0218\%$, $0.3249\%$, and $0.5909\%$, respectively, at the $N=384$ production grid. At $N=768$ the same transition gives $0.0223\%$, $0.3244\%$, and $0.5920\%$. An area-only extensive model")

sub("main.tex",
    r"Paired-contrast uncertainties over the available resolutions ($N=128$--$768$ at $Re=100$, $N=128$--$512$ at $Re=400$ and $1000$) are bounded by $7.2\times10^{-4}$, $2.4\times10^{-3}$, and $1.4\times10^{-2}$ percentage points, respectively; the $Re=100$ figure includes the $N=768$ level, the corresponding $N=128$--$384$ span being $3.4\times10^{-4}$, and the $Re=1000$ figure is the $N=128$--$512$ span, the $N=128$--$384$ span being $1.3\times10^{-2}$.",
    r"Paired-contrast uncertainties, taken as the span of $I_D(2,4)$ over all available resolutions ($N=128$--$768$ at each Reynolds number for this transition), are $7.2\times10^{-4}$, $2.4\times10^{-3}$, and $1.4\times10^{-2}$ percentage points at $Re=100$, $400$, and $1000$; the corresponding $N=128$--$384$ spans are $3.4\times10^{-4}$, $2.4\times10^{-3}$, and $1.3\times10^{-2}$ percentage points. The four-band partition recomputed at $N=768$ gives $\mathrm{d}A=+12.62$ and $\mathrm{d}B=+19.82$ at $Re=400$ and $+9.99$ and $+49.21$ at $Re=1000$, in units of $10^{-4}D^*(2,0)$; the net contrast differs from its $N=384$ value by $0.14\%$ and $0.19\%$, respectively.")

sub("main.tex",
    r"($N=128$--$768$ at $Re=100$, $N=128$--$512$ at $Re=400$ and $1000$), transitions, and amplitudes (Table~SM.X);",
    r"($N=128$--$768$ at all three Reynolds numbers for the primary transition), transitions, and amplitudes (Table~SM.X);")

sub("main.tex",
    r"Energy-balance audit at the $N=512$ extension ($Re=400$ and $1000$, $\alpha=0$) (Table~SM.XII);",
    r"Energy-balance audit at the $N=512$ extension ($Re=400$ and $1000$, $\alpha=0$) and at the $N=768$ extension ($n=2$ and $4$, $\alpha=0$, all three Reynolds numbers) (Table~SM.XII);")

sub("main.tex",
    r"with all four bands recomputed at every available resolution (up to $N=512$ where available) and the cancellation index,",
    r"with all four bands recomputed at every available resolution (up to $N=768$ at all three Reynolds numbers) and the cancellation index,")

# --------------------------------------------------------------- 04_results.tex
R = "sections/04_results.tex"

line(R, "Two robust facts follow.",
     r"Two features of the partition are robust. The added region contributes $dB>0$ at every Reynolds number, rising from $+6.99\times10^{-4}$ at $Re=100$ to $+4.87\times10^{-3}$ at $Re=1000$ on the production grid, and to $+6.96\times10^{-4}$ and $+4.92\times10^{-3}$ at $N=768$, in units of $D^*(2,0)$. It covers about $7.4\%$ of the total fluid area at $Re=400$ yet contributes only $0.20\%$ of $D^*$, so its mean dissipation density is about $2.7\%$ of the cavity average (about $1\%$ at $Re=100$ and $7\%$ at $Re=1000$). The redistribution term $dA$ is a small difference of much larger, opposing contributions from the four disjoint bands of the common region (Table~SM.XIII): at $Re=400$ the lid band ($y>0.9$) contributes $+2.44\times10^{-3}$ and the interior core $-2.83\times10^{-3}$, whose near-cancellation plus the bottom ($+6.21\times10^{-4}$) and side/corner ($+1.04\times10^{-3}$) bands leaves $dA=+1.28\times10^{-3}$; at $Re=1000$ the lid band is $+5.55\times10^{-3}$ and the core $-6.60\times10^{-3}$. The fluid added by flattening therefore sits in a slow, weakly sheared shoulder and dissipates little. The redistribution term is negative only at $Re=100$ ($-4.81\times10^{-4}$) and changes sign near $Re\approx145$; the added-zone term never changes sign. At $Re\ge400$ both terms are positive, so the small net contrast arises from cancellation inside $dA$, between the lid and core bands, rather than from cancellation between $dA$ and $dB$.")

line(R, "The grid robustness of this mechanism is assessed",
     r"The grid robustness of this mechanism is assessed in Table~SM.XV and Fig.~SM.5, where all four bands are recomputed cell by cell with the identical exact-partition code at every available resolution: $N=128$, $256$, $384$, $512$ (the last at $Re=400$ and $1000$) and $N=768$ for $n=2$ and $4$ at all three Reynolds numbers; the $N=384$ rows reproduce the headline partition to machine precision. The cancellation index $C=\sum|\text{bands}|/|dA|$ grows with Reynolds number and refinement, from $4.0$ at $Re=100$, to $5.0\to5.5$ at $Re=400$ and $11.3\to14.7$ at $Re=1000$ between $N=128$ and $N=768$, reflecting the increasing internal cancellation within $dA$. The net contrast $\Delta D=dA+dB$, which equals $I_D(2,4)$, remains stable to within $3.2\%$ ($Re=100$), $0.73\%$ ($Re=400$) and $2.5\%$ ($Re=1000$) of its value over all available resolutions, because the grid-induced changes in $dA$ and $dB$ partially cancel. At $N=768$ the partition gives $dA=-4.73$ and $dB=+6.96$ at $Re=100$, $+12.62$ and $+19.82$ at $Re=400$, and $+9.99$ and $+49.21$ at $Re=1000$, in units of $10^{-4}D^*(2,0)$; the net contrasts are $0.02232\%$, $0.32442\%$ and $0.59201\%$, within $2.5\%$, $0.14\%$ and $0.19\%$ of the $N=384$ values. The sign and order of magnitude of $dA$ and $dB$ are robust on every grid; their individual magnitudes carry a grid sensitivity rising with $Re$ (the $dA$ span is $1.9\%$, $9.0\%$ and $30\%$ of its value at $Re=100$, $400$ and $1000$), which reflects the internal cancellation.")

sub(R, r"This negative sign is resolved. The contrast is negative on all three grids",
       r"The contrast is negative on every computed grid")

sub(R, r"($-0.00371\%$, $-0.00360\%$, $-0.00371\%$ at $N=128$, $256$, $384$), with a span of $1.1\times10^{-4}$",
       r"($-0.00371\%$, $-0.00360\%$, $-0.00371\%$ and $-0.00364\%$ at $N=128$, $256$, $384$ and $512$), with a span of $1.1\times10^{-4}$")

sub(R, r"The $Re=1000$ value is therefore the first resolved negative trailing step.",
       r"The $Re=1000$ value is a resolved negative trailing step.")

sub(R,
    r"Uncertainty is the span of the contrast itself over the available resolutions at fixed ($Re$, transition): $N=128$--$768$ for $(2,4)$ at $Re=100$, $N=128$--$384$ for every other entry, extended to $N=128$--$512$ at $Re=400$ and $Re=1000$, $\alpha=0$. The $1.3\times10^{-2}$~pp entry for $(2,4)$ at $Re=1000$ is therefore the $N=128$--$384$ span; including the $N=512$ run raises it to $1.4\times10^{-2}$~pp (Table~SM.X), a difference that changes no conclusion.",
    r"Uncertainty is the span of the contrast itself over the available resolutions at fixed ($Re$, transition): $N=128$--$768$ for $(2,4)$ at all three Reynolds numbers, $N=128$--$384$ for the trailing transitions, extended to $N=128$--$512$ at $Re=400$ and $Re=1000$, $\alpha=0$. The $1.4\times10^{-2}$~pp entry for $(2,4)$ at $Re=1000$ is the $N=128$--$768$ span; over $N=128$--$384$ it is $1.3\times10^{-2}$~pp (Table~SM.X), a difference that changes no conclusion.")

sub(R,
    r"$1000$ & $(2,4)$ & $+0.5909$ & $+0.6022$ & $1.9$  & $+7.93$ & $1.3\times 10^{-2}$ \\",
    r"$1000$ & $(2,4)$ & $+0.5909$ & $+0.6022$ & $1.9$  & $+7.93$ & $1.4\times 10^{-2}$ \\")

line(R, "The primary contrast $I_D(2,4)$ at $\\alpha=0$ is evaluated",
     r"The primary contrast $I_D(2,4)$ at $\alpha=0$ is evaluated with the relaxation time held fixed per Reynolds number (diffusive scaling, $U\propto1/N$), so that the comparison isolates the grid. The production values at $N=384$ are $+0.02178\%$, $+0.32488\%$ and $+0.59087\%$ at $Re=100$, $400$ and $1000$; the additional levels are $N=128$, $256$, $512$ (the last at $Re=400$ and $1000$) and $N=768$ for $n=2$ and $4$ at all three Reynolds numbers. The $N=768$ values are $+0.02232\%$, $+0.32442\%$ and $+0.59201\%$. Between $N=256$ and $N=384$ the relative change is $-0.7\%$, $+0.08\%$ and $+0.54\%$; at $N=768$ the contrast lies $+2.5\%$, $-0.14\%$ and $+0.19\%$ from the production value. The sequence rises monotonically up to $N=384$ at $Re=400$ and $1000$ and continues to rise at $Re=1000$; at $Re=400$ the $N=768$ value falls slightly below the $N=384$ maximum, and at $Re=100$ the fourth significant figure is non-monotone. The spans over $N=128$--$768$ are $7.2\times10^{-4}$~pp ($3.2\%$ of the finest-grid value) at $Re=100$, $2.4\times10^{-3}$~pp ($0.73\%$) at $Re=400$ and $1.4\times10^{-2}$~pp ($2.5\%$) at $Re=1000$; over the three production grids alone they are $3.4\times10^{-4}$, $2.4\times10^{-3}$ and $1.3\times10^{-2}$~pp. These spans are the conservative numerical uncertainties quoted throughout. The absolute integrals are less stable: between $N=384$ and $N=768$ the reference dissipation rises by $0.86\%$, $0.84\%$ and $0.92\%$ at the three Reynolds numbers, so the two-finest-grid $0.2\%$ criterion is not met and $N=768$ remains a resolution extension rather than a converged mesh (Section~\ref{subsec:scope}).")

sub(R,
    r"between $N=256$ and $N=384$. These changes are larger than the changes in the contrast itself",
    r"between $N=256$ and $N=384$, and by $0.86\%$, $0.84\%$ and $0.92\%$ between $N=384$ and $N=768$ at $Re=100$, $400$ and $1000$. These changes are larger than the changes in the contrast itself")

line(R, "\\caption{(Color online) Resolution sensitivity of the primary contrast",
     r"\caption{(Color online) Resolution sensitivity of the primary contrast $I_D(2,4)$ at $\alpha=0$. Panels: (a) $Re=100$; (b) $Re=400$; (c) $Re=1000$. Each panel includes all computed levels, $N=128$, $256$, $384$ and $768$, together with $N=512$ at $Re=400$ and $1000$. The $N=768$ value lies $+2.5\%$, $-0.14\%$ and $+0.19\%$ from the $N=384$ value at the three Reynolds numbers. The spans over all available resolutions quoted in the text are the conservative numerical uncertainties on $I_D(2,4)$ used throughout.}")

line(R, "\\emph{(iii) Resolution.}",
     r"\emph{(iii) Resolution.} The span of $I_D(2,4)$ over all available resolutions is retained as the conservative numerical uncertainty: $7.2\times10^{-4}$~pp at $Re=100$ ($N=128$--$768$), $2.4\times10^{-3}$~pp at $Re=400$ and $1.4\times10^{-2}$~pp at $Re=1000$ (both $N=128$--$768$). At $N=768$ the absolute reference dissipation still changes by $0.84\%$ ($Re=400$) and $0.92\%$ ($Re=1000$) relative to $N=384$, so the two-finest-grid $0.2\%$ criterion is not met and no grid-converged limit is claimed. The refinement ratios are not uniform and the fitted effective order is unstable (Section~VII of the Supplementary Material); no Richardson extrapolation or grid-convergence index is quoted as a result. Over the production grids alone the spans are $3.4\times10^{-4}$, $2.4\times10^{-3}$ and $1.3\times10^{-2}$~pp, and are given wherever that distinction matters.")

sub(R,
    r"Because the finite-element results agree with the lattice-Boltzmann contrasts to within $0.002$ percentage point at the finest meshes, this independent solution confirms the primary contrast across solvers. A broader cross-solver sweep over the full parameter matrix remains a follow-up.",
    r"The finest-grid FEM--LBM gaps are $0.002$ percentage point at $Re=400$ and $0.008$ at $Re=1000$, smaller than the spread of the finite-element results across their own meshes. The comparison establishes cross-solver consistency; a sweep over the full parameter matrix would require additional finite-element calculations.")

sub(R,
    r"from $-6.2\%$ at $N=128$ to $-2.0\%$ at $N=384$ at $Re=400$ (Section~VII.A of the Supplemental Material).",
    r"from $-6.2\%$ at $N=128$ to $-2.0\%$ at $N=384$ and $-0.9\%$ at $N=768$ at $Re=400$ (Section~VII.A of the Supplemental Material); at $N=768$ the imbalance is $-0.89\%$ at $Re=400$ and $-0.85\%$ at $Re=1000$.")

sub(R, r"Crucially, this estimator agreement does not replace the resolution analysis:",
       r"This estimator agreement does not constrain the grid error:")

line(R, "The results above establish the robustness of the paired contrasts",
     r"Scope and limitations. (i) The span of $I_D(2,4)$ over all available resolutions ($N=128$--$768$ for $n=2$ and $4$, $\alpha=0$, at each Reynolds number) is retained as the conservative numerical uncertainty: $7.2\times10^{-4}$, $2.4\times10^{-3}$ and $1.4\times10^{-2}$~pp at $Re=100$, $400$ and $1000$. The $N=768$ runs are converged to residuals below $10^{-8}$ and are treated as resolution extensions, not converged meshes: the absolute integrals still change by $0.86\%$, $0.84\%$ and $0.92\%$ between $N=384$ and $N=768$, so the $0.2\%$ two-finest-grid criterion is not met. The contrast itself changes by $+2.5\%$, $-0.14\%$ and $+0.19\%$ over the same refinement. An indicative extrapolation is reported in Section~VII of the Supplementary Material: the fitted effective order varies from $0.70$ to $1.07$ for $D^*$ and is not stable under removal of the coarsest grid, so no asymptotic regime is established. (ii) The absolute integrals carry a residual grid dependence of order $1\%$ between $N=256$ and $N=384$, and about $0.9\%$ between $N=384$ and $N=768$; grid-independent absolute energetics are not claimed. (iii) The Taylor--Hood comparison is a cross-solver consistency check, not a convergence proof for either solver, because the finite-element variation across its own meshes exceeds the FEM--LBM gap. (iv) The scope is limited to the tested superellipse family, single undulation count $N_w=1$ and three Reynolds numbers; the unmatched Mach numbers of the coarse $n\ge6$ runs at $Re=100$ preclude a grid-converged assignment of the trailing-transition signs there. The $N=768$ series covers only $n=2$ and $4$, so the trailing transitions remain scoped to $N\le512$. (v) The $N=512$ extension is available at $Re=400$ and $Re=1000$, $\alpha=0$; it confirms the stability of the primary contrast to within $3\times10^{-5}$~pp between $N=384$ and $N=512$ at $Re=400$ and strengthens the negative sign of the $(6,8)$ transition at $Re=1000$. (vi) The $Re=1000$ results concern the steady two-dimensional solution of the MRT-LBM equations. The square lid-driven cavity undergoes a three-dimensional transition only above $Re$ of order $2000$; weak three-dimensional effects cannot be excluded near the square limit ($n=8$), so the conclusions at this Reynolds number describe the two-dimensional steady state.")

sub(R,
    r"The purpose is to make explicit which conclusions are resolution-backed and which remain qualitative or partial.",
    r"Each row gives the data supporting the claim and the resolution or modelling limitation attached to it.")

sub(R,
    r"extended to $N=128$--$512$ where that extension is available and to $N=128$--$768$ for the primary contrast at $Re=100$).",
    r"extended to $N=128$--$512$ for the trailing transitions at $Re=400$ and $1000$ and to $N=128$--$768$ for the primary contrast at all three Reynolds numbers).")

sub(R, r"$I_D(2,4)>0$ at $Re=100,400,1000$ & $N=384$ values; spans &",
       r"$I_D(2,4)>0$ at $Re=100,400,1000$ & $N=384$ and $N=768$ values; spans &")

sub(R, r"four bands at $N=128$--$512$ (Table~SM.XV)",
       r"four bands at $N=128$--$768$ (Table~SM.XV)", count=2)

# --------------------------------------------------------------- 05_conclusions.tex
C = "sections/05_conclusions.tex"

sub(C,
    r"but integrated dissipation rises by only $0.0218\%$, $0.3249\%$, and $0.5909\%$ at $Re=100$, $400$, and $1000$.",
    r"but integrated dissipation rises by only $0.0218\%$, $0.3249\%$, and $0.5909\%$ at $Re=100$, $400$, and $1000$ on the $N=384$ production grid. At $N=768$ the corresponding increases are $0.0223\%$, $0.3244\%$, and $0.5920\%$.",)

line(C, "An exact spatial partition separates",
     r"An exact spatial partition separates dissipation in the newly opened shoulder region from redistribution over the region common to both geometries. At $Re=400$, the added region occupies $7.4\%$ of the fluid area but contributes $0.2\%$ of the reference dissipation, while the common-region redistribution is a near-cancellation of opposing lid-band and core-band contributions. The common-region contribution is negative at $Re=100$ and positive at $Re=400$ and $1000$; the net contrast remains positive throughout. The four bands were recomputed with the same exact-partition code at $N=128$, $256$, $384$, $512$ where available and $N=768$ at all three Reynolds numbers (Table~SM.XV, Fig.~SM.5). At $N=768$ the partition gives $\mathrm{d}A=-4.73$ and $\mathrm{d}B=+6.96$ at $Re=100$, $+12.62$ and $+19.82$ at $Re=400$, and $+9.99$ and $+49.21$ at $Re=1000$, in units of $10^{-4}D^*(2,0)$; the corresponding contrasts are $0.02232\%$, $0.32442\%$ and $0.59201\%$, within $2.5\%$, $0.14\%$ and $0.19\%$ of the $N=384$ values. The cancellation index rises from $4.0$ at $Re=100$ to $14.7$ at $Re=1000$, while the net contrast remains stable. The cell-by-cell identity $\Delta D=\mathrm{d}A+\mathrm{d}B$ is algebraically exact and closes to machine precision at every $N$. The numerical accuracy of the individual terms is weaker, because $\mathrm{d}A$ is the residual of larger opposing band contributions and its relative uncertainty grows with the cancellation index; the discrete partition as a whole is stable to a few percent, but the separate continuum limits of $\mathrm{d}A$ and $\mathrm{d}B$ are not established. The grid sensitivity of $\mathrm{d}A$ reaches $1.9\%$, $9.0\%$ and $30\%$ of its value at $Re=100$, $400$ and $1000$. At $Re=100$, the energy-balance defect varies by $5.0$ points between $N=128$ and $N=768$ while $I_D(2,4)$ remains within $7.2\times10^{-4}$~pp (Table~SM.XIX), showing that the defect cancels in the paired contrast.")

line(C, "The response is concentrated in the $n=2\\to4$ transition",
     r"The response is concentrated in the $n=2\to4$ transition at all Reynolds numbers studied, which carries about $96.2\%$ of the summed absolute contrasts at $Re=400$. The full exponent sequence is available at $N\le384$ at all three Reynolds numbers: at $\alpha=0$ the $n=4\to6$ contrast is $+0.00296\%$, $+0.0129\%$, and $+0.0103\%$, and the $n=6\to8$ contrast is $+0.00042\%$, $+6\times10^{-5}\%$, and $-0.0037\%$ at $Re=100$, $400$, and $1000$. Resolution is judged transition by transition against each contrast's own span. At $Re=1000$ both trailing steps are resolved, the $(6,8)$ step by a factor of $34$ with the same negative sign on all computed grids; at $Re=400$ that step lies inside its own span, and at $Re=100$ the $(4,6)$ step mixes Mach regimes on the coarse grids and is quoted as an upper bound. At $\alpha=0.1$ both trailing transitions are negative and individually resolved at $Re=400$ and $1000$. Increasing undulation amplitude weakens the primary contrast. The primary $\alpha=0$ contrast varies by at most $3.2\%$ of its own value over $N=128$--$768$ at each Reynolds number, while the absolute integrals retain roughly $1\%$ grid sensitivity at the production level and about $0.9\%$ between $N=384$ and $N=768$. Independent Taylor--Hood calculations reproduce the primary contrast to $0.002$ percentage points at $Re=400$ and $0.008$ at $Re=1000$ on the finest meshes of each solver; the finite-element spread across its own meshes exceeds these gaps, so the comparison establishes consistency rather than convergence. Scope is limited to the tested superellipse family and a single undulation count; the unmatched Mach numbers of the coarse $n\ge6$ runs at $Re=100$ and the absence of $n=6,8$ fields at $N=768$ preclude universal scaling or grid-independent claims for the trailing transitions. Numerical uncertainties on $I_D(2,4)$ are $7.2\times10^{-4}$, $2.4\times10^{-3}$ and $1.4\times10^{-2}$ percentage points at $Re=100$, $400$ and $1000$, each measured over $N=128$--$768$; the spans over $N=128$--$384$ are $3.4\times10^{-4}$, $2.4\times10^{-3}$ and $1.3\times10^{-2}$~pp. The trailing contrasts have their own, generally smaller spans, collected in Table~SM.X.")

line(C, "Several targeted extensions would quantify",
     r"Several targeted extensions would reduce the conservative spans reported here. A converged Richardson fit, with an effective order stable under removal of the coarsest grid, would replace the span by an extrapolated value with an explicit grid-convergence index; the present sequences do not provide such an order. Extending the $N=768$ series to $n=6$ and $8$ would give the trailing transitions the same resolution treatment as the primary one. A $Re=100$ Taylor--Hood calculation would test cross-solver agreement at the low-Reynolds end. Finally, a finer-grid Reynolds sweep would establish whether the crossover $Re\approx346$, currently a descriptive feature of the $N=128$ sweep, survives refinement.")

# --------------------------------------------------------------- supplement
S3 = "supplementary/sections/s3_convergence.tex"

sub(S3, r"The convergence history for representative cases is shown in Fig.~\ref{suppfig:residual}.",
     r"The convergence history for representative cases is shown in Fig.~\ref{suppfig:residual}. Six further runs at $N=768$ ($n=2$ and $4$, $\alpha=0$, $Re=100$, $400$ and $1000$) reached residuals between $9.7\times10^{-9}$ and $9.96\times10^{-9}$ after $1.17$--$3.38\times10^{6}$ iterations.")

sub(S3, r"against $+0.5909\%$ at $N=384$ in the lattice-Boltzmann production set:",
       r"against $+0.5909\%$ at $N=384$ in the lattice-Boltzmann production set and $+0.5920\%$ at $N=768$:")

S6 = "supplementary/sections/s6_contrasts.tex"
sub(S6,
    r"The $N=512$ data at $Re=1000$ confirm the monotonicity of the $(2,4)$ contrast ($+0.592\%$) and strengthen the negative sign of the $(6,8)$ transition ($-0.0036\%$).",
    r"The $N=512$ data at $Re=1000$ confirm the monotonicity of the $(2,4)$ contrast ($+0.592\%$) and strengthen the negative sign of the $(6,8)$ transition ($-0.0036\%$). The $N=768$ level extends the primary transition at all three Reynolds numbers, giving $+0.02232\%$, $+0.32442\%$ and $+0.59201\%$ at $Re=100$, $400$ and $1000$; only $n=2$ and $4$ were run at that resolution, so the trailing transitions remain scoped to $N\le512$.")

S7 = "supplementary/sections/s7_consistency.tex"
sub(S7, r"The $N=768$ extension at $Re=100$ continues the decay, reaching $-0.86\%$ for $n=2$.",
       r"The $N=768$ extension continues the decay, reaching $-0.86\%$ at $Re=100$, $-0.89\%$ at $Re=400$ and $-0.85\%$ at $Re=1000$ for $n=2$.")

sub(S7, r"; $0.32251\%$, $0.32462\%$, $0.32488\%$, $0.32487\%$ at $Re=400$;",
       r"; $0.32251\%$, $0.32462\%$, $0.32488\%$, $0.32487\%$, $0.32442\%$ at $Re=400$;")

sub(S7, r"and $0.57758\%$, $0.58769\%$, $0.59087\%$, $0.59176\%$ at $Re=1000$",
       r"and $0.57758\%$, $0.58769\%$, $0.59087\%$, $0.59176\%$, $0.59201\%$ at $Re=1000$")

sub(S7, r"(values at $N=128$, $256$, $384$, $512$; $N=512$ available only at $Re=400$ and $Re=1000$)",
       r"(values at $N=128$, $256$, $384$, $512$, $768$; $N=512$ available only at $Re=400$ and $1000$)")

sub(S7, r"A fourth level, $N=768$, was computed at $Re=100$ (converged, residual $9.9\times10^{-9}$) and gives $+0.02232\%$.",
       r"The fifth level, $N=768$, is converged at all three Reynolds numbers (residuals below $10^{-8}$) and gives $+0.02232\%$, $+0.32442\%$ and $+0.59201\%$.")

sub(S7, r"($N=128$--$768$ at $Re=100$, $N=128$--$512$ at $Re=400$ and $1000$)",
       r"($N=128$--$768$ at all three Reynolds numbers)", count=2)

sub(S7, r"($N = 128$--$768$ at $Re=100$, $N = 128$--$512$ at $Re=400$ and $1000$)",
       r"($N = 128$--$768$ at all three Reynolds numbers)")

sub(S7, r"and by a further $+0.86\%$ between $N=384$ and $N=768$ at $Re=100$.",
       r"and by a further $+0.86\%$, $+0.84\%$ and $+0.92\%$ between $N=384$ and $N=768$ at $Re=100$, $400$ and $1000$.")

sub(S7,
    r"The $N=512$ data at $Re=1000$ confirm the monotonicity of the $(2,4)$ contrast ($+0.592\%$ at $N=512$) and strengthen the negative sign of the $(6,8)$ transition ($-0.0036\%$ at $N=512$).",
    r"The $N=512$ data at $Re=1000$ confirm the monotonicity of the $(2,4)$ contrast ($+0.592\%$ at $N=512$) and strengthen the negative sign of the $(6,8)$ transition ($-0.0036\%$ at $N=512$); at $N=768$ the primary contrast reaches $+0.5920\%$ while the trailing transitions remain unavailable at that resolution.")

sub(S7, r"$N=512$ was not run at $Re=100$; the fourth level $N=768$ (converged, $Re=100$ only) is listed as a resolution extension.",
       r"$N=512$ was not run at $Re=100$; the $N=768$ level (converged for $n=2$ and $4$ at all three Reynolds numbers) is listed as a resolution extension.")

sub(S7, r"up to $\sim 27\%$ at $Re=1000$", r"up to $\sim 30\%$ at $Re=1000$")

sub(S7, r"up to $\sim27\%$ at $Re=1000$", r"up to $\sim30\%$ at $Re=1000$")

sub(S7, r"the $N=768$ point there is the converged resolution extension",
       r"the $N=768$ points are converged resolution extensions")

sub(S7, r"$11.3\to14.3$ ($Re=1000$) from $N=128$ to $N=512$",
       r"$11.3\to14.7$ ($Re=1000$) from $N=128$ to $N=768$")

sub(S7, r"At $Re=100$ the fourth level is the most informative:",
       r"At $Re=100$ the finest level is the most informative:")

sub(S7,
    r"At $Re=400$ the $N=512$ extension leaves the span unchanged ($2.4\times10^{-3}$~pp); at $Re=1000$ it widens the span to $1.4\times10^{-2}$~pp over $N=128$--$512$ ($2.4\%$ of the finest-grid value). At $Re=100$ the converged $N=768$ level widens the span to $7.2\times10^{-4}$~pp ($3.2\%$): the contrast itself remains inside a few percent on four grids, whereas the absolute integrals still change by $0.87\%$ between $N=384$ and $N=768$, so the $0.2\%$ two-finest-grid criterion is not met and the fourth level is reported as a resolution extension rather than as evidence of grid convergence.",
    r"With $N=768$ included, the spans become $2.4\times10^{-3}$~pp at $Re=400$ and $1.4\times10^{-2}$~pp ($2.5\%$ of the finest-grid value) at $Re=1000$, both over $N=128$--$768$; at $Re=100$ the span widens from $3.4\times10^{-4}$ to $7.2\times10^{-4}$~pp ($3.2\%$). The absolute integrals still change by $0.84\%$ ($Re=400$) and $0.92\%$ ($Re=1000$) between $N=384$ and $N=768$, so the $0.2\%$ two-finest-grid criterion is not met and the finest level is reported as a resolution extension.")

sub(S7,
    r"using the available sequences ($N=128,256,384,768$ at $Re=100$; $N=128,256,384,512$ at $Re=400$ and $1000$).",
    r"using the available sequences ($N=128,256,384,768$ at $Re=100$; $N=128,256,384,512$ at $Re=400$ and $1000$). The $N=768$ point is not used in these fits, so that the reported orders remain comparable across Reynolds numbers.")

sub(S7,
    r"A fourth LBM level at $N=768$ for $n=2$ and $n=4$ at $Re=100$ ($\alpha=0$) has been completed, and both runs are converged:",
    r"The $N=768$ level for $n=2$ and $n=4$ ($\alpha=0$) is converged at all three Reynolds numbers. At $Re=100$ the runs give:")

sub(S7, r"The four-level sequence also settles the question left open by the intermediate checkpoint:",
       r"The four-level sequence resolves the intermediate-checkpoint excursion:")

sub(S7,
    r"The Taylor--Hood FEM runs at $Re=1000$ for $n=2$ and $n=4$ up to $n_x=256$ are complete (Section~\ref{suppsec:crosscheck}).",
    r"At $Re=400$ and $Re=1000$ the $N=768$ runs give $D^*=23.62949$ and $31.83541$ for $n=2$, and $23.70615$ and $32.02387$ for $n=4$, with residuals below $10^{-8}$; the primary contrasts are $+0.32442\%$ and $+0.59201\%$, against $+0.32488\%$ and $+0.59087\%$ at $N=384$. The Taylor--Hood FEM runs at $Re=1000$ for $n=2$ and $n=4$ up to $n_x=256$ are complete (Section~\ref{suppsec:crosscheck}).")

S9 = "supplementary/sections/s9_repro.tex"
sub(S9,
    r"raw converged LBM fields (244 \texttt{.npz}, including the fourth level at $N=768$ for $n=2$ and $n=4$ at $Re=100$)",
    r"raw converged LBM fields (244 \texttt{.npz}, including the $N=768$ level for $n=2$ and $n=4$ at $Re=100$)")

sub(S9,
    r"A DOI-bearing Zenodo snapshot of the \texttt{v1.1.0} tree will be added upon acceptance;",
    r"The four $N=768$ scalar records at $Re=400$ and $1000$ produced for the present revision are stored in \texttt{results/scalars/}; their raw fields postdate the v1.1.0 archive and are not included in it. A DOI-bearing Zenodo snapshot of the \texttt{v1.1.0} tree will be added upon acceptance;")

S10 = "supplementary/sections/s10_defect_analysis.tex"
sub(S10,
    r"and by a further factor of $2.3$ to $0.87\%$ at $N=768$ and $Re=100$,",
    r"and by a further factor of about $2.3$ to $0.87\%$ at $N=768$ and $Re=100$ ($0.89\%$ at $Re=400$ and $0.85\%$ at $Re=1000$),")

sub(S10,
    r"$2.4\times10^{-3}$~pp ($0.7\%$) at $Re=400$ over $N=128$--$384$.",
    r"$2.4\times10^{-3}$~pp ($0.7\%$) at $Re=400$ over $N=128$--$768$.")

sub(S10,
    r"$\approx 4.2$ points at $Re=400$ (over $N=128$--$384$)",
    r"$\approx 5.3$ points at $Re=400$ (over $N=128$--$768$)")

# --------------------------------------------------------------- SM contrast table
T = "supplementary/contrasts_table.tex"

sub(T,
    r"the range of $I_D$ over the available grids, $N=128$, $256$, $384$, and $512$ (the last resolution available only at $Re=400$ and $Re=1000$, $\alpha=0$), plus the fourth level $N=768$ for the $(2,4)$ contrast at $Re=100$ and $\alpha=0$.",
    r"the range of $I_D$ over the available grids, $N=128$, $256$, $384$, $512$ (the last available only at $Re=400$ and $Re=1000$, $\alpha=0$) and $768$, the last being available for the $(2,4)$ contrast at $\alpha=0$ at all three Reynolds numbers.")

insert_after(T, r"400  & 512 & $(2,4)$ & 0.00",
             [r"400  & 768 & $(2,4)$ & 0.00 & $+0.32442$ & $+0.32563$ &  0.37 & $+0.00122$ & $+0.00237$ \\"])

sub(T,
    r"1000 & 512 & $(2,4)$ & 0.00 & $+0.59176$ & $+0.59987$ &  1.37 & $+0.00811$ & $+0.01418$ \\",
    r"1000 & 512 & $(2,4)$ & 0.00 & $+0.59176$ & $+0.59987$ &  1.37 & $+0.00811$ & $+0.01418$ \\" + "\n" +
    r"1000 & 768 & $(2,4)$ & 0.00 & $+0.59201$ & $+0.59813$ &  1.03 & $+0.00612$ & $+0.01443$ \\")

# --------------------------------------------------------------- SM defect table
D = "supplementary/sections/s10_defect_analysis.tex"
# The final table rows are handled by replacing the three-row block.
sub(D,
    r"$400$ & $128$ & $6.204$ & $6.193$ & $0.011$ & $0.32251$ \\" + "\n" +
    r"      & $256$ & $3.199$ & $3.194$ & $0.005$ & $0.32462$ \\" + "\n" +
    r"      & $384$ & $2.041$ & $2.038$ & $0.003$ & $0.32488$ \\",
    r"$400$ & $128$ & $6.204$ & $6.193$ & $0.011$ & $0.32251$ \\" + "\n" +
    r"      & $256$ & $3.199$ & $3.194$ & $0.005$ & $0.32462$ \\" + "\n" +
    r"      & $384$ & $2.041$ & $2.038$ & $0.003$ & $0.32488$ \\" + "\n" +
    r"      & $768$ & $0.892$ & $0.891$ & $0.001$ & $0.32442$ \\")

print("edited:", ", ".join(sorted(set(EDITED))))
