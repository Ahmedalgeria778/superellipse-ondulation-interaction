# 02 — Protocole : quantifier l'interaction (différence de différences)

## 1) Objectif

Mesurer si l'effet d'une ondulation locale `(alpha, Nw)$` sur les grandeurs
énergétiques dépend de la forme globale `n$` de la cavité. Test par
**différence de différences** (interaction factorielle à 2×2 minimal).

## 2) Grandeurs cibles

On se concentre sur les indicateurs énergétiques issus du bilan mécanique :
- `K = (1/2) ∫_Ω |u|^2 dA$` (énergie cinétique, adim par $U_ref^2 H^2$ ou conforme aux conventions)
- `P = ∫_{Γ_{lid}} u_w · (σ n)\,ds$` (puissance fournie par le lid)
- `D = ∫_Ω 2 S:S\,dA$` (dissipation visqueuse totale)

Critère de convergence numérique : `|P-D|/D <= delta$` avec `delta$` faible
(cible typique < 1-2% à N=256, < 0.5% à N=384). On exige ce critère pour tous
les cas retenus dans l'analyse d'interaction.

## 3) Différence simple (effet de l'ondulation pour un n donné)

```
ΔQ(n, alpha) = Q(n, alpha) - Q(n, 0),   Q ∈ {K,P,D}.
```

Interprétation : variation induite par l'ondulation à forme globale fixée.

## 4) Différence de différences (interaction)

Pour deux formes contrastées `n_a, n_b$` :

```
I_Q(n_a, n_b, alpha) = ΔQ(n_b, alpha) - ΔQ(n_a, alpha)
                     = [Q(n_b,alpha)-Q(n_b,0)] - [Q(n_a,alpha)-Q(n_a,0)].
```

- `I_Q ≈ 0$` → effet additif (H0). L'ondulation agit de façon quasi identique
  pour `n_a$` et `n_b$`.
- `I_Q != 0$` → interaction significative (H1). L'effet de l'ondulation **dépend**
  de la forme globale.
- Amplitude relative (facultatif) : `i_Q = I_Q / ( (1/2)(|ΔQ(n_a)|+|ΔQ(n_b)|) )$` si les dénominateurs non nuls.

## 5) Campagne minimale (pilote)

Cas A–D à `Re=400, N=256, lid plat, H=1, a/H=0.5, b/H$ fixé` :

| Cas | n | alpha | Nw | Description |
|---|---:|---:|---:|---|
| A  | 4 | 0.00  | 1  | fond superelliptique lisse (référence n=4) |
| B  | 4 | 0.10  | 1  | n=4 + ondulation modérée |
| C  | 8 | 0.00  | 1  | fond quasi-rectangulaire (référence n=8) |
| D  | 8 | 0.10  | 1  | n=8 + ondulation modérée |

Interaction : `n_a=4, n_b=8, alpha=0.10$`.

Objectif pilote : obtenir `I_P, I_D$` avec incertitude numérique estimée. Si
`|I_Q| > 2*sigma_{num}$`, on retient H1 comme plausible.

## 6) Estimation d'incertitude numérique (minimal)

Approche pragmatique (compatible avec les deux bases) :

1. **Convergence spatiale** : exécuter `N=256, 384, 512$` sur **un sous-ensemble**
   représentatif (au moins A et D, idéalement tous). Estimer erreur par
   Richardson/ordre apparent ou par différence relative `|Q_{2N}-Q_N|/|Q_N|$`.
2. **Bilan** : `r = (P-D)/D$` est un proxy direct d'erreur de fermeture.
3. **Stabilité itérative** : résidu LBM `res < 1e-8$` (seuil cohérent avec Projet A/B).
4. **Croisement FEM** (validation qualitative) : comparer topologie de tourbillons
   et `K,D$` sur 1-2 cas critiques (FEM Taylor–Hood P2/P1 copié).

On rapportera `sigma_{num}(I_Q)$` par propagation (différences) à partir des
incertitudes estimées sur `Q$`.

## 7) Décomposition spatiale (mécanisme)

Pour éclairer H2, on calculera `D_{near-wall}$ vs `D_{core}$ avec **définition unifiée** :

- `d(x, Γ)$ = distance minimale à toute frontière fluide (`Γ = lid ∪ parois ∪ fond effectif y_b$`).
- Seuil `delta_{w} = c * h$ avec `c \in [1,3]$` (maillage), ou `delta_w/H = const$`
  (indépendant de N) pour comparer. À fixer après premier test.
- `D_{nw} = ∫_{Ω : d < delta_w} 2S:S dA$, `D_{core} = D - D_{nw}$.
- Examiner `ΔD_{nw}(n,alpha), ΔD_{core}(n,alpha)$ et leurs interactions.

L'objectif n'est **pas** de reproduire exactement les chiffres Projet B, mais
d'appliquer une **même définition** à tous les cas combinés.

## 8) Conditions strictes d'analyse

- **Même solveur** pour A–D (pas de mélange de codes). Seul le **même jeu de CL,
  normalisation, diagnostics** garantit la validité de `I_Q$`.
- **Pas de régression** : valider `(n=4, α=0)$` vs Projet B à `Re=400,N=256$`
  sur `D, K$` (ordre de grandeur + tendance). Écart toléré raisonnable (1-3%)
  compte tenu de différences (lid plat vs tangent, définitions). **P ne peut
  pas** être comparé directement — c'est justement ce que ce projet mesure.
- **Traçabilité** : enregistrer pour chaque cas `case_id = se_n{n}_w{alpha}_Nw{Nw}_Re{Re}_N{N}` et tous les scalaires dans `results/`.
- **Documentation** : justifier tout écart > seuil dans `docs/03_plan_campagne.md`.
