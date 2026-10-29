# 03 — Plan de campagne numérique

## 1) Objectifs

1. Implémenter géométrie combinée + solveur LBM-MRT (héritage A+B).
2. Valider sur cas de référence (`α=0`, `n in {4,8}`, `Re=400`) vs Projet B (D,K).
3. Réaliser campagne pilote 4 cas (A–D), calculer `I_P, I_D`, estimer incertitude.
4. Étendre selon résultat (Re, N, n supplémentaires).

## 2) Paramètres fixes

| Paramètre | Valeur (adim H=1) | Justification |
|---|---|---|
| `H` | 1.0 | référence Projet A |
| `L = 2a` | 1.0 → `a = 0.5` | cavité globale 1×1 centrée |
| `x` | `[-0.5, 0.5]` | symétrie naturelle superellipse |
| `b/H` | 0.5 (provisoire) | profondeur modérée, surface fluide comparable entre n ; à confirmer par check_geometry |
| `U_ref` | 1.0 | lid plat |
| `ν` | `U_ref H / Re` | cohérent |
| `lid` | plat horizontal `u=(1,0), y=1` | choix A (isole interaction) |
| `latérales x=±0.5` | no-slip, verticales | choix A |
| `fond y_b` | no-slip, courbe | combiné |

## 3) Campagne pilote (prioritaire)

| Cas | ID | n | alpha | Nw | Re | N | État |
|---|---|---|---|---|---|---|---|
| A | `se_n4_w0p00_Nw1_Re400_N256` | 4.0 | 0.00 | 1 | 400 | 256 | à exécuter |
| B | `se_n4_w0p10_Nw1_Re400_N256` | 4.0 | 0.10 | 1 | 400 | 256 | à exécuter |
| C | `se_n8_w0p00_Nw1_Re400_N256` | 8.0 | 0.00 | 1 | 400 | 256 | à exécuter |
| D | `se_n8_w0p10_Nw1_Re400_N256` | 8.0 | 0.10 | 1 | 400 | 256 | à exécuter |

Cible résidu < 1e-8, `|P-D|/D < 0.02`. Générer : scalaires (`results/scalars/*.json`), champs (`results/fields/*.npz`), figures (`figures/validation/`).

## 4) Contrôle de résolution (grid convergence)

Sur A et D (extrêmes attendus) : `N ∈ {256, 384, 512}`. Estimer incertitude sur `D,P,K,I_Q`.

## 5) Extension (si interaction détectée)

- `Re ∈ {100, 400, 1000}`
- `n ∈ {2, 4, 6, 8}` (couverture plus large)
- `alpha ∈ {0.05, 0.10}` (sensibilité amplitude)

## 6) Validation croisée

1. **Référence B (α=0)** : comparer A,C vs `se_n4_Re400_N256`, `se_n8_Re400_N256` du Projet B. Métriques : `D, K, |u_max|, psi_primary`. Tolérance cible ±2-3% (diff lid plat vs tangent).
2. **FEM** : 1-2 cas (Taylor–Hood P2/P1, `fem_taylor_hood.py` copié). Comparaison qualitative/topologique + scalaires.
3. **Contrôle physique (et non symétrie)** : la cavité entraînée n'est PAS symétrique gauche-droite — la réflexion `x -> -x` transforme la vitesse du lid `(1,0)` en `(-1,0)`, donc les CL ne sont pas invariantes, et le tourbillon principal est décalé du côté aval (sens du lid). On suit donc la position du tourbillon principal comme contrôle quantitatif vs le benchmark de Ghia et al. (1982) (Re=400 : centre `~ (0.5547, 0.6055)`).

## 7) Livrables

- `audit/inventaire_cas.csv`, `audit/synthese_inventaire.json`, `docs/00_synthese_audit.md`
- `geometry/*` + `figures/geometry/check_*.png` + rapport validation
- `solver/*`, `diagnostics.py`, `checks.py`, `cases.py`, `run_case.py`
- `results/scalars/pilot_4cas.csv`, `results/scalars/interaction_IQ.json`
- `docs/04_resultats_pilote.md` (à rédiger après exécution)
- `provenance/ORIGINS.md` (traçabilité)

## 8) Critère d'arrêt avant extension

Poursuivre uniquement si `|I_P| > 2*sigma_{num}(I_P)$ **ou** `|I_D| > 2*sigma_{num}(I_D)$ (interaction non-nulle au seuil d'incertitude). Sinon conclure additivité (H0).

---

## 9) État d'avancement (18/10/2026)

Le plan ci-dessus est le document de cadrage initial. L'état réel :

| étape | état | note |
|---|---|---|
| Campagne pilote 4 cas (A–D) | **fait** | requalifiée : `n ∈ {2,4,6,8}`, `α ∈ {0, 0.10}` |
| Extension §5 | **fait** | `Re ∈ {100, 400, 1000}` × 8 cas = **24 cas à N=256**, 0 échec, 223,9 min |
| Tolérance | **1e-8 adoptée** | test A/B : écart 0,002 % sur `I_D`, gain 21 % d'itérations (note 06) |
| Contrôle de résolution §4 | **partiellement fait** | test N=128→256 sur les 3 Re (note 06) ; **N=512 reporté** (mesuré : ~49 h pour le Tier 1, non compatible avec exploration parallèle) |
| Cartographie `Re` | **en cours** | Re tous les 50 de 150 à 950 à N=128, 20 valeurs au total (note 08) |
| Validation FEM §6.2 | **à faire** | |
| `α ∈ {0.05}` §5 | **non fait** | amplitude intermédiaire non explorée |

Résultats détaillés : notes 05, 06, 07, 08.
