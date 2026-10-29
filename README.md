# Superelliptic bottom-shaping vs. lid-driven-cavity dissipation

**Bottom-wall shaping redistributes rather than adds viscous dissipation in a lid-driven cavity: a superelliptic energy-budget analysis.**

Reproducible source code, manuscript, figures and derived results for the study of how a
superelliptic bottom family interacts with a local bed undulation in a lid-driven cavity,
quantified through the mechanical energy budget (`P*`, `D*`, `K*`) and an exact cell-by-cell
`dA/dB` dissipation partition.

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

---

## Scientific question

> Does the energetic effect of a local bed undulation depend on the global shape of a
> lid-driven cavity, and can this interaction be explained by a redistribution of viscous
> dissipation between near-wall regions and the flow core?

Testable hypotheses:

- **H0 — additivity**: the undulation effect is independent of the global shape `n`.
- **H1 — geometric interaction**: the undulation effect depends on `n`.
- **H2 — parietal mechanism**: any interaction is associated with a redistribution of
  dissipation near the bed / in the core (not assumed true).

Interaction test (difference of differences), for two shapes `n_a`, `n_b` and an undulation
`(alpha, Nw)`:

```
I_Q = [Q(n_b, alpha) - Q(n_b, 0)] - [Q(n_a, alpha) - Q(n_a, 0)]
```

applied separately to `Q = P` (power input) and `Q = D` (dissipation), with identical `Re`,
resolution, boundary conditions and normalisation. See `docs/02_protocole_interaction.md`.

## Headline result

Flattening the bottom (`n = 2 -> 4`) grows the fluid area by **7.93 %**, yet total viscous
dissipation rises by only **+0.0218 %** (`Re = 100`), **+0.3249 %** (`Re = 400`) and
**+0.5909 %** (`Re = 1000`). An area-only extensive model overpredicts these contrasts by
factors of **13 to 364**. The exact `dA/dB` partition explains the sub-extensive response:
the newly opened shoulder region adds fluid that barely dissipates, while the top near-lid
layer concentrates **74–85 %** of all dissipation, and the common-region redistribution `dA`
is itself a Reynolds-dependent near-cancellation (negative at `Re = 100`, zero-crossing near
`Re ≈ 145`, saturating above `Re ≈ 600`).

## Method at a glance

- **Geometry**: superelliptic bottom family `n ∈ {2, 4, 6, 8}` with optional sinusoidal bed
  undulation `alpha ∈ {0, 0.1}`; flat moving lid, fixed vertical side walls.
- **Solver**: D2Q9 multiple-relaxation-time lattice Boltzmann with Bouzidi bounce-back;
  resolutions `N ∈ {128, 256, 384}` (selected cases to `N = 512`).
- **Cross-validation**: independent conforming Taylor–Hood P2/P1 finite-element solver with
  native element-wise dissipation diagnostics.
- **Diagnostics**: energy-budget indicators `P*`, `D*`, `K*`; local dissipation density
  `phi = 2(S11² + S22² + 2S12²)`; four-band `dA/dB` partition (lid / bottom / side / core).

## Repository layout

```
.
├── README.md, LICENSE, CITATION.cff, requirements.txt
├── .gitignore, .gitattributes
├── solver/          # LBM-MRT core, geometry, diagnostics, case runner, analysis scripts
├── geometry/        # combined superellipse + undulation geometry and validation
├── fem/             # Taylor–Hood P2/P1 finite-element cross-validation solver
├── audit/           # unified case inventory + generator
├── docs/            # audit, geometry spec, interaction protocol, campaign plan, per-stage notes
├── results/         # derived, version-controlled outputs (scalars, JSON/CSV, per-run summaries)
├── paper/           # manuscript source + compiled PDF (main + supplementary) + figures
└── provenance/      # timestamped, hashed copies of reused files from the two base projects
```

Root-level post-processing scripts: `extract_dAdB_multires_v2.py` (multi-resolution `dA/dB`
partition), `make_dadb_fig.py` (resolution-dependence figure), `postprocess_dAdB.py`.

## Environment

Python 3.14+, `numpy>=1.26`, `numba>=0.61`, `scipy>=1.13`, `matplotlib>=3.9`
(see `requirements.txt`).

```bash
python -m venv .venv
source .venv/bin/activate      # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

## Reproducing the analysis

1. **Run a campaign** (or reuse `results/`): the case matrix and executor live in `solver/`
   (`run_case.py`, `cases.py`, `run_campaign.py`). Field output is written to `results/fields/`.
2. **Compute contrasts and the `dA/dB` partition**:
   ```bash
   python -m solver.audit_contrasts
   python -m solver.audit_balance
   python -m solver.referee_mech_analysis
   python extract_dAdB_multires_v2.py       # full four-band partition, N = 128..512
   ```
3. **Regenerate figures** used in the paper: scripts in `solver/` (`make_*_figs.py`,
   `plot_*.py`) plus `make_dadb_fig.py`.

> The scripts expect the `.npz` field arrays produced in step 1 under `results/fields/`
> and `results/fem/`. Those raw arrays are **not committed** (see Data availability);
> re-running a campaign regenerates them.

## Building the paper

LaTeX with `revtex4-2` (MiKTeX or TeX Live):

```bash
cd paper
pdflatex main.tex && bibtex main && pdflatex main.tex && pdflatex main.tex
cd supplementary
pdflatex supp_main.tex && bibtex supp_main && pdflatex supp_main.tex && pdflatex supp_main.tex
```

Both compiled PDFs (`paper/main.pdf`, `paper/supplementary/supp_main.pdf`) are committed for
reference. Tables are numbered `SM.I–SM.XIX` and figures `SM.1–SM.6` in the supplement.

## Data availability

Raw converged velocity-field arrays (242 LBM `.npz`) and native FEM fields (28 `.npz`),
together ≈ 278 MB, are distributed as a **GitHub Release asset** rather than inside the git
tree, so that `git clone` stays lightweight:

- Release: <https://github.com/Ahmedalgeria778/superellipse-ondulation-interaction/releases/tag/v1.0.0>
- Direct download: <https://github.com/Ahmedalgeria778/superellipse-ondulation-interaction/releases/download/v1.0.0/superellipse_lbm_fem_fields_v1.0.zip>

The archive (`superellipse_lbm_fem_fields_v1.0.zip`) contains `fields/`, `fem/`, a
`DATA_README.md` (file naming convention and npz keys) and a `MANIFEST.txt` (size + truncated
SHA-256 per file). All version-controlled `results/` content in this repository (scalar
budgets, contrast tables, band partitions, audit JSON) is derived directly from the paper's
reported numbers. To regenerate the raw fields instead, re-run the campaign with the scripts
above.

## Provenance

This work reuses solver components from two earlier private projects
(`Lid_driven_ondulation`, `superellipse_LBM_PI_operator`). Those originals are never modified;
the reused files are copied under `provenance/` with timestamps and hashes recorded in
`provenance/ORIGINS.md`.

## Citation

If you use this code or data, please cite the repository and the associated paper — see
[`CITATION.cff`](CITATION.cff).

## License

Released under the **MIT License**. See [LICENSE](LICENSE).
