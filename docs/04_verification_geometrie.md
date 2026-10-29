# 04 — Vérification géométrique des 4 cas pilotes (Option A)

Généré par `geometry/check_geometry.py`.

## Résultats

| Cas | ID | n | alpha | y_b(±a) | y_b_min | A_f (num, Nx=2048) | A_f_analytic_base (α=0) | Écart | slope_max_abs (dy/dx) |
|---|---|---|---|---|---|---|---|---|---|
| A | se_n4_w0p00_Nw1_Re400_N256 | 4.0 | 0.00 | 0.0 | -0.5 | 1.46348 | 1.46352 | ~4e-5 | ~108 |
| B | se_n4_w0p10_Nw1_Re400_N256 | 4.0 | 0.10 | 0.0 | -0.44888 | 1.41348 | 1.46352 | -0.05004 (–3.42%) | ~108 |
| C | se_n8_w0p00_Nw1_Re400_N256 | 8.0 | 0.00 | 0.0 | -0.5 | 1.48912 | 1.48923 | ~1e-4 | ~128 |
| D | se_n8_w0p10_Nw1_Re400_N256 | 8.0 | 0.10 | 0.0 | -0.47940 | 1.43912 | 1.48923 | -0.05011 (–3.36%) | ~128 |

## Observations

- `y_b(±a)=0` pour tous les cas : condition de fermeture verticale (extrémités plates) satisfaite.
- `y_b_min >= -b = -0.5$` et augmente légèrement pour n=8 + α=0.10 (effet de l'ondulation).
- Pente `|dy/dx|$` élevée près des "coins" arrondis (superellipse) — attendu. À prendre en compte pour résolution N (Bouzidi efficace ; contrôle par N=256/384/512).
- Variation de surface fluide `A_f$` entre α=0 et α=0.10 : **~3.4%** pour chaque n (diminution car ondulation remonte le fond localement en moyenne). C'est modéré. Pour l'interaction `I_Q = [Q_b-Q_{b,0}]-[Q_a-Q_{a,0}]$, la différence de différences devrait être peu biaisée par cet effet de volume ; on documentera néanmoins dans l'analyse.

## Figures

- `figures/geometry/se_n4_w0p00_Nw1_Re400_N256_profile.png`
- `figures/geometry/se_n4_w0p10_Nw1_Re400_N256_profile.png`
- `figures/geometry/se_n8_w0p00_Nw1_Re400_N256_profile.png`
- `figures/geometry/se_n8_w0p10_Nw1_Re400_N256_profile.png`

## Conclusion

Géométrie combinée **validée** (continuité aux bords, formes cohérentes). Passage à la phase suivante : adaptation du solveur LBM-MRT (héritage Projet A + éléments Projet B) pour calculer **P, D, K** avec diagnostics complets.
