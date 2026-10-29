# Lid-Driven Cavity with a Corrugated Floor — Mechanical Energy Budget

Numerical study of the steady incompressible flow in a two-dimensional lid-driven
cavity whose **bottom wall is sinusoidally corrugated**, with a controlled
**moving-lid velocity profile**. The work quantifies two effects on the global
mechanical-energy budget:

- **(A)** the geometrical effect of the wall corrugation, and
- **(B)** the kinematic effect of the prescribed lid profile.

The reference solver is a D2Q9 **multiple-relaxation-time (MRT) lattice
Boltzmann** scheme with an **interpolated (Bouzidi) bounce-back** treatment of the
curved wall. An independent **Taylor–Hood P2/P1 finite-element** solver (steady
Navier–Stokes, Picard linearisation) is used for cross-validation.

Author: **Ahmed Beniaiche**, Laboratory of Fluid Mechanics, Ecole Militaire
Polytechnique, Bordj El Kiffan, Algiers, Algeria.
Contact: `beniaiche_ahmed25@yahoo.fr` — ORCID [0000-0002-9613-7049](https://orcid.org/0000-0002-9613-7049).

---

## Overview

The cavity is the unit square (`L = H = 1`). The bottom wall is

```
y_b(x) = (h/2) * [ 1 - cos(2*pi*Nw*x/L) ],   Nw in {1, 2, 3, 4},   alpha = h/H.
```

The amplitude ratio is `alpha in {0.025, 0.05, 0.075, 0.10}`; the sixteen
combinations `(alpha, Nw)` are labelled `G01..G16`, with the flat-wall reference
`F00`. Three lid profiles are compared at fixed Reynolds number: a **uniform**
profile, a **corner-regularised smoothstep** profile (`eps = 0.02`), and the same
smoothstep **renormalised to unit mean speed** (protocol B).

The global dimensionless diagnostics are

```
K* = 0.5 * int |u*|^2 dA*              (kinetic energy)
D* = 2 * int S*:S* dA*                 (viscous dissipation; S from finite differences)
P* = int U*_lid (du*/dy*)|lid dx*      (mechanical power injected by the lid)
balance = (P* - D*) / D*               (closure residual, -> 0 at steady state)
```

## Key results

1. **The corrugation effect is small but stable.** At `Re = 400` the corrugation
   reduces `P*` by `-0.17%` (G06) to `-0.45%` (G12), and `K*` by `-0.32%` to
   `-1.14%`. The differential quantity varies by `<= 0.011` percentage points
   between `N = 128` and `512`, and is insensitive to the numerical variant
   (`consistent` vs `Mach`-fixed, `<= 0.053%`). The effect grows with `Re` and
   with the amplitude `h`.
2. **The lid-profile effect dominates.** At equal mean lid speed, the spatial
   distribution of the imposed velocity changes `P*` by about `-16.3%` and `D*` by
   about `-24.4%`; renormalising the mean speed adds about `-3.75` percentage
   points. The uniform lid also degrades the global power balance
   (`-12.6%` vs `-3.2%` for smooth profiles), which is interpreted as a difference
   in numerical closure quality, not in physical efficiency.
3. **The absolute budget is not fully grid-converged** over `N = 128..512`
   (`P*` still moves by a few percent). Physical conclusions therefore rest on
   matched-resolution differences and Richardson diagnostics rather than on
   extrapolated absolute values.

A full, self-contained analysis of the data is in
[`REPORT.md`](REPORT.md).

## Repository structure

```
src/                     Solver and analysis code
  geometry.py            Undulated wall, fluid mask, Bouzidi boundary links
  lbm_wavy.py            D2Q9 MRT lattice-Boltzmann solver (Numba kernel) + WavyMRT
  diagnostics.py         psi/omega, K*, P*, D*, vortex detection
  cases.py               48-case matrix (G01..G16 x Re) and control groups
  run_case.py            Run one case, save fields and table rows to results/
  reanalyze.py           Re-analyse saved fields into the summary table
  fem_taylor_hood.py     Taylor-Hood P2/P1 finite-element solver (validation)
  wavy_effect.py         Differential effect dQ = (Q_wavy - Q_flat)/Q_flat
  convergence_order.py   Observed orders and Richardson extrapolation
  scaling_law.py         Empirical correlation of the corrugation effect
  make_tables.py         LaTeX/Markdown tables from the analysis CSVs
  make_figures.py        Manuscript figures (PNG + PDF)
  ...
configs/                 (reserved)
data/                    Reference data (Ghia et al. 1982)
results/                 Steady-state fields (*.npz) and analysis CSVs
results_fneq/            Non-equilibrium dissipation cross-check
figures/                 Generated figures
manuscrit/               LaTeX manuscript (AIP / REVTeX 4.2) + compiled PDFs
REPORT.md                Detailed internal analysis report
```

## Requirements and installation

Python 3.10+ with the packages in [`requirements.txt`](requirements.txt)
(`numpy`, `numba`, `scipy`, `matplotlib`).

```bash
python -m venv .venv
# Windows
.venv\Scripts\python -m pip install -r requirements.txt
# Linux / macOS
.venv/bin/python -m pip install -r requirements.txt
```

The LaTeX manuscript requires a TeX distribution with REVTeX 4.2 and the AIP
substyle (`revtex4-2`, `aip4-2`, `aipnum4-2`); the manuscript directory ships the
needed style files.

## Usage

Run a single LBM case and save its fields:

```bash
.venv\Scripts\python src\run_case.py G06-R400 --N 256
```

Run the finite-element solver on one case:

```bash
.venv\Scripts\python src\fem_taylor_hood.py FEM-G06-R400 --nx 64 --ny 64 --N 160
```

Regenerate the analysis tables and figures:

```bash
.venv\Scripts\python src\make_tables.py     # writes manuscrit/tables/*.tex and *.md
.venv\Scripts\python src\make_figures.py    # writes figures/*.png and *.pdf
```

Build the manuscript:

```bash
cd manuscrit
pdflatex manuscrit && bibtex manuscrit && pdflatex manuscrit && pdflatex manuscrit
```

## Case identifiers

```
G01..G16    increasing h/H and Nw; main cases x Re = 48
F00-R{Re}   flat-wall reference (regularised lid), baseline for P*
EQ05/EQ10   equal-area flat wall (volume control)
L1-...      uniform vs regularised lid sensitivity
V01..V03    flat-wall validation (uniform lid)
FEM-...     finite-element cases
```

## Extracted quantities

See the formulas in [Overview](#overview). Vortices are detected as local extrema
of the stream function `psi` (primary plus weak secondary eddies); their count is
reported as a qualitative topological diagnostic only, since the weakest
secondary eddies lie close to the detection threshold.

## Validation

- Flat-wall LBM reproduces the reference MRT solver to `max|du| ~ 1e-6`.
- Primary-vortex centre agrees with Ghia et al. (1982):
  `V01-R100`: `(0.620, 0.734)` FEM / `(0.613, 0.738)` LBM / Ghia `(0.617, 0.734)`;
  `V03-R1000`: `(0.526, 0.568)` FEM / `(0.535, 0.566)` LBM / Ghia `(0.531, 0.563)`.
- LBM vs FEM on `G06-R400`: `K*` and `P*` agree to better than `2%`.
- The non-equilibrium dissipation estimate `D*_fneq` reproduces `D*_grad u` to
  `0.15%`, confirming that the measured dissipation is not an artefact of the
  finite-difference strain-rate estimate.

## Data

Steady-state fields are stored as compressed NumPy archives (`*.npz`) under
`results/`; summary quantities are in `results/scalar_table.json` (78 validated
rows, all `converged` with `checks_ok = True`) and in the CSVs under
`results/analysis/` and `results/convergence/`.

## Manuscript

The manuscript (`manuscrit/manuscrit.tex`, AIP/REVTeX 4.2, single-column preprint)
and its Supplementary Material (`manuscrit/supplementary_material.tex`) are
included with their compiled PDFs. Tables are generated from the analysis CSVs by
`src/make_tables.py`; no number is transcribed by hand.

## Citation

If you use this code or data, please cite the manuscript and this repository as
described in [`CITATION.cff`](CITATION.cff).

## License and copyright

- **Source code** (`src/`, scripts): MIT License — see [`LICENSE`](LICENSE).
- **Data** (`results/`, `results_fneq/`, `data/`, `figures/`): Creative Commons
  Attribution 4.0 International (CC BY 4.0) — see
  [`LICENSE-DATA`](LICENSE-DATA).
- **Manuscript and its figures/tables** (`manuscrit/`): © 2026 Ahmed Beniaiche.
  All rights reserved. Not covered by the code or data licences above.

See [`NOTICE`](NOTICE) for the full breakdown.
