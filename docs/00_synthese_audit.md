# 00 — Synthèse d'audit des deux bases (Étape 1)

Objectif : établir un **inventaire unifié** et les **écarts méthodologiques**
entre `Lid_driven_ondulation` (Projet A) et `superellipse_LBM_PI_operator`
(Projet B) avant de concevoir la géométrie combinée.

## 1) Matériel et statut

- Projet A — cavité entraînée à fond ondulé. Solveur **D2Q9-MRT** + **Bouzidi**
  (frontières courbes), lid régularisé. Diagnostics complets : `K, P, D`,
  **bilan mécanique** `r = (P-D)/D`. Décomposition **near-wall/core** disponible
  (`diag_fneq.py`). FEM Taylor–Hood P2/P1 présente. Champs `.npz` complets.
- Projet B — cavité **superelliptique** `|x/a|^n + |y/b|^n < 1`. Solveur courbe
  (Bouzidi). **190 cas** (10 n × 19 Re) à N=128 + raffinements N=256/384/512.
  Diagnostics : `D, D_core, K, enstrophie, vorticité, aire de masque vs exacte`.
  **Absence totale de P (pression/traction pariétale/puissance exacte)**
  (`pressure_available = False` dans tous les cas de production).
- **Stratégie** : fusionner **l'apport A** (P + décomposition énergétique) et
  **l'apport B** (générateur de fond superelliptique, traitement des coins
  arrondis/nettes selon n, infrastructure LBM courbe) dans un **nouveau solveur**
  dédié, **sans modifier** les projets d'origine.

## 2) Inventaire unifié

Généré par `audit/build_inventory.py` — sorties :
`audit/inventaire_cas.csv`, `audit/synthese_inventaire.json`.

| Total | Projet A | Projet B |
|---|---:|---:|
| 382 | 99 | 283 |

Répartition :
- A : `wavy_square` (81) + `wavy_square_conv` (18) ; Re∈{100,400,1000} ; α∈{0,0.025,0.05,0.075,0.10} ; Nw∈{1,2,3,4} ; N∈{128,256,384,512}.
- B : `superellipse` (283 uniques) ; n∈{1.25,...,8.0} ; Re∈{100..1000} ; N majoritairement {128,256,384,512} selon corpus. Tous avec champs `.npz`.

## 3) Données disponibles par projet

| Indicateur | Projet A | Projet B | Remarque |
|---|---|---|---|
| `K` (énergie cinétique) | oui | oui | normalisations à vérifier si on compare entre projets |
| `D` (dissipation totale) | oui | oui | définitions spatiales proches ; vérifier conventions |
| `P` (puissance lid/tractions) | **oui** | **non** | écart bloquant pour un bilan `P-D=0` complet côté B |
| Bilan `(P-D)/D` | oui (balances < 5%) | inconnu | à reproduire dans le nouveau solveur |
| `D_core` / décomposition | oui (`diag_fneq`) | oui (`analyse_coeur_paroi`) | définitions différentes → **ne pas réutiliser directement** |
| Champs `ux,uy,rho,psi` | oui | oui | formats `.npz` ; layouts à inspecter si on veut réanalyser B |
| Aire fluide | oui (fluid_area) | oui (area_mask, area_exact) | erreur masque croît avec `n` (corners) : 0.19% (n=2) → ~0.95% (n=8) |

> **Conclusion tables B** : `physics_observables.csv` est riche pour `D,K` et
> statistiques, mais **n'inclut ni P, ni drag, ni traction pariétale**. Toute
> campagne combinée devra donc calculer `P = ∫_{lid} u_w · (σ n) ds`
> explicitement (approche Projet A à adapter).

## 4) Conventions géométriques observées (à unifier)

| Élément | Projet A (carré + fond ondulé) | Projet B (superellipse seule) |
|---|---|---|
| Domain | `[0,1]×[y_b(x),1]` typiquement | `Ω_n = {|x/a|^n + |y/b|^n < 1}`, `a=b` souvent |
| Fond | `y_b(x; α,Nw)= (α/2)(1-cos(2π Nw x))` s'annule à x=0,1 si centré | arc superelliptique `y_b0(x;n)= -b[1-(|x|/a)^n]^{1/n}` |
| Lid | horizontal mobile (u_top) ; forme profilée selon variantes | arc supérieur "tangent" (mode `upper`) dans la construction B |
| Parois latérales | verticales (x=const) + fond ondulé | frontière curviligne continue (superellipse) |
| Point de référence | `H=1`, `U_ref` normalisation | `L_char`, `U_lid` en unités LBM |
| Re | `Re = U_ref H / ν` | `Re = U_lid L_char / ν_lbm` (avec conventions spécifiques) |

## 5) Scénarios pilotes (cible)

Alignés sur la demande :
- `n_a = 4` vs `n_b = 8` (formes contrastées : quasi-rectangulaire vs plus "piquée")
- Ondulation `(α=0, 0.10)`, `Nw=1`
- `Re=400`, résolution `N=256` (contrôle minimal à poursuivre)
- Lid plat horizontal pour tous les cas (choix A validé).

Liste pilotée : `(n=4, α=0), (n=4, α=0.10, Nw=1), (n=8, α=0), (n=8, α=0.10, Nw=1)`

## 6) Risques identifiés

1. **Dénominateurs** : incohérences de nommage cas B (`n4p0` vs `n4`) → déjà corrigées dans l'inventaire.
2. **Normalisation** : `U_ref/H` vs `U_lid/L_char`. À imposer **systématiquement** dans le nouveau solveur (un seul jeu de conventions).
3. **P manquant** : le nouveau solveur **doit calculer P** sur lid plat (intégrale de contrainte de cisaillement/pression sur frontière supérieure) ; réutiliser formules diagnostics Projet A.
4. **Décomposition near-wall/core** : définitions distinctes. Préférer **implémenter une seule** (distance à toute frontière, seuil unifié) pour comparer proprement.
5. **Aire/masque** : erreur croissante avec `n` (coins très pointus, n grand). Nécessite inspection maillage + convergence N.
6. **Reproduction** : ne **pas** extrapoler B→combiné. **Valider** sur sous-ensemble plat (`α=0`) vs B pour n=4,8 Re=400 avant d'activer ondulation.

## 7) Prochaines étapes (cf plan)

**Étape 2** : géométrie combinée (`geometry/`) — `superellipse.py` (existant, copié), `combined_geometry.py` (fond `y_b0 + eta`), `wall_functions.py`, `configs.py`, `check_geometry.py`. Démarrage immédiat.
