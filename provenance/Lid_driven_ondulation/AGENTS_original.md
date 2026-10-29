# AGENTS.md

## Always recompile after any `.tex` edit

After editing any `.tex` source, you MUST rebuild the PDF **and** replace the
deliverable copies. Do not stop at editing the source.

```powershell
# Main manuscript (run in manuscrit/)
pdflatex -interaction=nonstopmode manuscrit.tex
bibtex manuscrit
pdflatex -interaction=nonstopmode manuscrit.tex
pdflatex -interaction=nonstopmode manuscrit.tex

# Supplementary material (run in manuscrit/)
pdflatex -interaction=nonstopmode supplementary_material.tex
bibtex supplementary_material
pdflatex -interaction=nonstopmode supplementary_material.tex
pdflatex -interaction=nonstopmode supplementary_material.tex
```

Deliverable PDFs in `manuscrit/`: `manuscrit.pdf`, `manuscrit_aip.pdf`
(identical copies), `supplementary_material.pdf`.

## PDF lock workaround (very common here)

`manuscrit.pdf` is frequently locked by **`explorer.exe`** (Explorer preview
pane), not only by Adobe Acrobat. Symptoms: `! I can't write on file
'manuscrit.pdf'` with `no output PDF file produced`, even though the build
appears to "succeed" (the old PDF still exists).

Workflow:
1. Build to a temporary jobname to sidestep the lock:
   `pdflatex -jobname=build_m ...` → `build_m.pdf`.
2. Verify `build_m.pdf`, then replace the deliverable.
3. If the target is locked, find the locker with the Restart Manager API
   (`rstrtmgr.dll`); if it is `explorer.exe`, run
   `Stop-Process -Name explorer -Force`, copy the file, then
   `Start-Process explorer.exe`.
4. Delete temporary `build_*` artefacts afterwards (`build_mNotes.bib` etc.).

Notes:
- `pdflatex` writes a `.log` with the **job name**; never redirect stdout to the
  same name as the log.
- The `MiKTeX: major issue: So far, you have not checked for MiKTeX updates`
  line is a harmless warning, not an error. Check `Output written on` in the
  `.log` to confirm success.

## Repo layout

- `manuscrit/manuscrit.tex` — main text (revtex4-2, `[aip,pof,preprint]`,
  `\bibliographystyle{aipnum4-2t}`).
- `manuscrit/supplementary_material.tex` — supplement (article, `[11pt]`,
  `\bibliographystyle{aipnum4-2t}`, S-numbered tables/figures).
- `manuscrit/aipnum4-2t.bst` — **patched copy of the AIP bst** (see below).
- `manuscrit/tables/*.tex` — table fragments, **auto-generated** by
  `src/make_tables.py` (never edit generator and output inconsistently).
- `figures/*.pdf|png` — figures, **auto-generated** by `src/make_figures.py`.
- `src/*.py` — solver post-processing and figure/table generation.

## Editorial conventions enforced

- **Figure files are numbered in appearance order.** `fig03_convergence` is
  Fig. 3, `fig04_fields` is Fig. 4, etc. Keep file names, labels
  (`fig:figNN_*`), and the call order in `make_figures.py` aligned.
- **Nomenclature** (main text) lists Roman, Greek, and Abbreviation entries.
  Greek section ends page 2; `\newpage` then Abbreviations on page 3
  (`longtable` is incompatible with revtex4-2 — use two `tabular`s).
- **Every abbreviation must be defined at its first occurrence**, including in
  the abstract and conclusion: `D2Q9`, `MRT`, `LBM`, `RMS`, `RMSE`, `LOGO`,
  `LAAO`.
- **`K_{\max}` is renamed `\kappa_{\max}`** (wall curvature) everywhere to avoid
  confusion with `K^{\ast}` (kinetic energy).
- `\delta_c` (vortex-centre distance) is used in Table III — keep it in the
  nomenclature.
- Abstract must stay **< 250 words** (journal-style count; math counts as one).
- Reference `choaire2026` is verified real (DOI 10.55592/cilamce2025.v5i.14082,
  published 2026, CILAMCE vol. 5). Do not remove it.

## Bibliography (PoF style)

- Reference format must match *Physics of Fluids*: authors, **article title in
  quotes**, abbreviated journal, **bold volume**, **inclusive pages**, year —
  e.g. `U. Ghia, ..., "High-Re solutions ...," J. Comput. Phys. 48, 387–411
  (1982).`
- `aipnum4-2t.bst` is `aipnum4-2.bst` with the `control.init` defaults changed
  to `control.title := 0` and `control.pages := 1` (article titles + page
  ranges). This is needed because the bundled AIP substyle disables titles by
  default and only re-enables them via a `@CONTROL` entry that in this MiKTeX
  install is written only when `\bibliographystyle` is **not** called — which is
  fragile and impossible in the `article`-class supplement. The patched bst makes
  main text and supplement identical with no warnings.
- Journal names in `references.bib` are stored **abbreviated** (AIP style); the
  bst passes them through verbatim. Full names such as "Journal of Computational
  Physics" would render unabbreviated, so keep the abbreviations.
- Do not add `\nocite{aip41Control}` to the supplement: the control entry then
  fires but natbib emits spurious "Citation undefined" warnings.

## Author / submission facts

- Author: **Ahmed Beniaiche**, ORCID 0000-0002-9613-7049.
- Affiliation (identical in main text and supplement): *Laboratory of Fluid
  Mechanics, Ecole Militaire Polytechnique, Algiers, Algeria*.
- Corresponding-author marker: `\author{Ahmed Beniaiche\textsuperscript{a)}}`
  with `\email[Corresponding author: ]{beniaiche\_ahmed25@yahoo.fr}`.
- Data availability: "available from the corresponding author upon reasonable
  request" (deliberate user choice — no public deposit).
- Licences: code MIT (`LICENSE`), data CC BY 4.0 (`LICENSE-DATA`), manuscript
  © 2026 Ahmed Beniaiche.
