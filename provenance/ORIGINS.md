# ORIGINS — traçabilité des fichiers réutilisés

Ce projet **ne modifie aucun fichier** des deux bases d'origine. Les fichiers
ci-dessous ont été **copiés** dans `provenance/` comme référence en lecture seule.
Les chemins d'origine sont donnés à titre d'audit.

- Projet A : `C:\Users\user\Documents\local_projects\Default Project\Lid_driven_ondulation`
- Projet B : `C:\Users\user\Documents\local_projects\Default Project\superellipse_LBM_PI_operator`

Empreinte SHA-256 tronquée à 16 caractères hexadécimaux ; taille en octets.

## Projet A — `Lid_driven_ondulation` (fond ondulé)

Chemin d'origine des sources : `<A>\src\<fichier>` sauf mention contraire.

| Fichier (`provenance/Lid_driven_ondulation/`) | Taille | SHA-256 (16) |
|---|---:|---|
| geometry.py | 8178 | EDD0E7046FAD1E46 |
| lbm_wavy.py | 12271 | E14C441500417ADF |
| diagnostics.py | 11147 | 949FBA207DB2373F |
| diag_fneq.py | 5148 | 728E18592CD25AD4 |
| run_case.py | 12828 | 313EB160017DF40D |
| cases.py | 3543 | 3B4B132B43B72FAC |
| checks.py | 5718 | 9B6D08369950BEE0 |
| fem_taylor_hood.py | 18587 | 6F02AB43BD5514FF |
| reanalyze.py | 6340 | E984B78BD1035605 |
| analyze.py | 6625 | D364B5EEE10F03E5 |
| make_figures.py | 21416 | E5F75D74C515B630 |
| make_tables.py | 14975 | 4B71F7943F3760D3 |
| convergence_study.py | 7760 | 365FE89979754A3C |
| convergence_order.py | 4685 | C83FE03B18E50007 |
| wavy_effect.py | 3673 | 93582C578E6E2B75 |
| README_original.md | 8053 | 276B7BD8E01E7E5F |
| REPORT_original.md | 13761 | 69C28E59E426AC81 |
| AGENTS_original.md | 5364 | 2F794E90A983FBAA |
| CITATION_original.cff | 1607 | 9F4C460035F1C951 |

## Projet B — `superellipse_LBM_PI_operator` (cavité superelliptique)

Les sous-dossiers de `provenance/superellipse_LBM_PI_operator/` reproduisent
l'arborescence d'origine (`01_geometry`, `02_LBM`, `06_physics`,
`06_neural_operator`). Chemin d'origine : `<B>\<sous-dossier>\<fichier>`.

| Fichier (`provenance/superellipse_LBM_PI_operator/`) | Taille | SHA-256 (16) |
|---|---:|---|
| 01_geometry/superellipse.py | 20747 | 5EF812AEBBE799F8 |
| 02_LBM/solver_curved.py | 8177 | E1802002F9CA1CE0 |
| 02_LBM/mrt_collision.py | 3100 | AD306C0CE56928BF |
| 02_LBM/geometry_links.py | 2930 | 0DB0CAEFF2E3CEAA |
| 02_LBM/streaming.py | 2033 | BE3C7FDC8D0FD6BA |
| 02_LBM/bouzidi_linear.py | 2037 | B9FADDF54A71DB03 |
| 02_LBM/moving_wall.py | 2493 | 027CF19A35A5CA53 |
| 02_LBM/lbm_geometry_grid.py | 13523 | BFB8A33FA0C916CD |
| 02_LBM/lbm_mrt_reference_square.py | 51887 | CFAA3A55D03E28EA |
| 06_physics/analyse_coeur_paroi.py | 11471 | 9EA7EF47315CCB7C |
| 06_physics/analyse_etape2_3_coeur.py | 9886 | 66265425C856E795 |
| 06_physics/donnees_papier.py | 32783 | E1708651434B45E0 |
| 06_physics/diag_dissipation_alg.py | 9517 | 68E8134C26D8BDF1 |
| 06_neural_operator/physics_observables.py | 30275 | 9C0D56138B39838F |

## Fichiers recopiés dans l'arborescence de travail

| Fichier de travail | Origine |
|---|---|
| `geometry/superellipse.py` | Projet B `01_geometry/superellipse.py` (classe d'évaluation de la frontière superelliptique) |
| `fem/fem_taylor_hood.py` | Projet A `src/fem_taylor_hood.py` (solveur FEM P2/P1) |

## Interpréteur utilisé (lecture seule)

`C:\Users\user\Documents\local_projects\Default Project\Lid_driven_ondulation\.venv\Scripts\python.exe`
(Python 3.14.5, numpy 2.5.4, numba 0.68.0, scipy 1.18.1, matplotlib 3.11.2).
Utilisé uniquement pour exécuter les scripts de **ce** projet.
