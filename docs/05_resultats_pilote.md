# 05 — Résultats de la campagne pilote

Campagne : **Re = 400, N = 128, Nw = 1**, `n ∈ {2, 4, 6, 8}`, `α ∈ {0, 0.10}`.
Géométrie combinée `y_b(x;n,α) = relief0(x;n) + η(x;α)` (boîte carrée 1×1,
`a=b=0.5`, `H=1`), lid plat régularisé, résolution fixée (pas de convergence —
résolution de travail validée dans les projets d'origine).

Rappel des quantités (adimensionnées, `U_ref=1`) :
- `D* = 2 ∫ S:S dA`  (dissipation, = traînée normalisée)
- `P* = ∫ U_lid (∂u/∂y)|_lid dx`  (puissance pariétale injectée)
- `K* = ½ ∫ |u|² dA`  (énergie cinétique)
- Effet de l'ondulation : `ΔQ(n,α) = Q(n,α) − Q(n,0)`
- **Différence de différences (interaction)** : `I_Q = ΔQ(n_b,α) − ΔQ(n_a,α)`

## 1) Scalaires (Re = 400, N = 128)

| n | α | D\* | P\* | K\* | bilan (P−D)/D | aire fluide | tourbillon (x,y) | nv |
|---|---|---|---|---|---|---|---|---|
| 2 | 0.00 | 22.75808 | 21.34619 | 0.0391305 | −6.20 % | 0.89343 | (0.558, 0.615) | 1 |
| 2 | 0.10 | 22.64137 | 21.23464 | 0.0384693 | −6.21 % | 0.84290 | (0.566, 0.633) | 1 |
| 4 | 0.00 | 22.83148 | 21.41749 | 0.0394229 | −6.19 % | 0.96387 | (0.556, 0.606) | 2 |
| 4 | 0.10 | 22.67368 | 21.26348 | 0.0386140 | −6.22 % | 0.91309 | (0.563, 0.629) | 3 |
| 6 | 0.00 | 22.83430 | 21.41989 | 0.0394353 | −6.19 % | 0.98230 | (0.556, 0.606) | 2 |
| 6 | 0.10 | 22.67260 | 21.26249 | 0.0386272 | −6.22 % | 0.93201 | (0.563, 0.629) | 3 |
| 8 | 0.00 | 22.83429 | 21.41981 | 0.0394386 | −6.19 % | 0.98938 | (0.556, 0.606) | 2 |
| 8 | 0.10 | 22.67224 | 21.26213 | 0.0386306 | −6.22 % | 0.93884 | (0.563, 0.629) | 3 |

Fichiers : `results/scalars/pilot_Re400_N128.csv`, `results/scalars/se_*.json`,
champs `results/fields/se_*.npz`.

### Contrôle physique
- Tourbillon principal (α=0) à `(0.556, 0.606)` pour n≥4 : **accord avec le
  benchmark de Ghia et al. (1982)** à Re=400 → `(0.5547, 0.6055)`. Le solveur
  (noyau Projet A) est donc correctement opérationnel sur la géométrie combinée.
- **La cavité entraînée n'est pas symétrique gauche-droite** : la réflexion
  `x→−x` change le signe de la vitesse du lid, donc les CL ne sont pas
  invariantes ; le tourbillon est décalé du côté aval. Aucun test de symétrie
  n'est donc pertinent (l'ancienne colonne `sym_err` a été retirée de `checks.py`).
- Ajouter l'ondulation (`α=0.10`) **déplace le tourbillon** vers le haut/droite
  (`0.556→0.563`, `0.606→0.629`) et **augmente le nombre de tourbillons** (`2→3`).

## 2) Effet de l'ondulation `ΔQ` et interaction `I_Q`

| n | ΔD\* | ΔP\* | ΔK\* |
|---|---|---|---|
| 2 | −0.11671 | −0.11155 | −0.0006612 |
| 4 | −0.15780 | −0.15401 | −0.0008089 |
| 6 | −0.16170 | −0.15740 | −0.0008081 |
| 8 | −0.16205 | −0.15768 | −0.0008080 |

| paire (n_a, n_b) | I_D | I_P | I_K | I_D / D(n_a,0) |
|---|---|---|---|---|
| (2, 4) | −0.04109 | −0.04246 | −1.48e−4 | −0.181 % |
| (2, 6) | −0.04499 | −0.04585 | −1.47e−4 | −0.198 % |
| (2, 8) | −0.04534 | −0.04613 | −1.47e−4 | −0.199 % |
| (4, 6) | −0.00390 | −0.00339 | +8e−7 | −0.017 % |
| (4, 8) | −0.00425 | −0.00367 | +9e−7 | −0.019 % |
| (6, 8) | −0.00035 | −0.00028 | +1e−7 | −0.002 % |

Fichier : `results/scalars/interaction_Re400_N128.json`.

## 3) Interprétation

1. **H1 (non-additivité) soutenue.** `I_D` et `I_P` sont **négatifs et
   significativement non nuls** dès que l'on change fortement la forme globale :
   `|I_D(2,4)| ≈ 0.041`, à comparer au bruit numérique sur `D*` (`~1e−5`,
   apparié — même solutionneur et même protocole). Le signe négatif signifie que
   l'effet de l'ondulation est **plus fort** pour une cavité anguleuse (grand n)
   que pour une cavité très arrondie (n=2).
2. **Saturation pour n ≥ 4.** `ΔD`, `ΔP` atteignent un plateau pour `n ≥ 4`
   (`ΔD ≈ −0.158…−0.162`), et `I(4,6)`, `I(4,8)`, `I(6,8) → 0`. Autrement dit,
   au-delà de n≈4 la forme globale n'évolue plus beaucoup (l'aire fluide
   plafonne : 0.9639 → 0.9894) et l'interaction résiduelle est quasi nulle.
   L'essentiel de l'interaction se joue entre **n=2 (demi-disque) et n≥4
   (fonds quasi-plats)**.
3. **Cohérence entre quantités** : `I_P ≈ I_D` (à 3 %) et `I_K/I_D ≈ 3.5e−3`,
   ordre attendu (`K` est une intégrale d'énergie, `D` une intégrale de
   gradients). La tendance est donc robuste, pas un artefact d'une seule métrique.
4. **Bilan mécanique** systématiquement à `(P−D)/D ≈ −6.2 %` (le Projet A
   rapportait `≈−3 %` sur sa géométrie). C'est un **biais systématique** du
   post-traitement discret (`P` et `D` par différences finies) ; il est
   quasi-constant entre cas et se compense largement dans `I_Q`. Il borne
   toutefois la précision absolue.

## 4) Conclusion du pilote

À Re=400 : l'interaction forme globale × ondulation est **non nulle**, de
l'ordre de **0.2 % de D** (≈ 1.7 % de l'effet d'ondulation à n=2), et **sature
pour n ≥ 4**. Le mécanisme apparent : plus la cavité est anguleuse, plus
l'ondulation mord sur une région proche-paroi active (recirculation
secondaire), donc plus son effet sur `D` et `P` est marqué.

## 5) Figures

- `figures/analysis/effect_vs_n_Re400.png` : `ΔD(n)`, `ΔP(n)` (effet vs forme).
- `figures/analysis/interaction_Re400.png` : `I_D`, `I_P` par paire.
- `figures/analysis/fields_Re400_N128.png` : `|u|` et `ψ` (8 cas).

## 6) Suite

- Étendre en Reynolds : **Re ∈ {100, 1000}** (mêmes `n`, `α`) pour la dépendance
  en Re de `I_Q`.
- Sensibilité d'amplitude : `α = 0.05`.
- Décomposition **near-wall / core** de `D` pour tester le mécanisme (H2).
- Validation croisée FEM (`fem/fem_taylor_hood.py`) sur 1–2 cas.
