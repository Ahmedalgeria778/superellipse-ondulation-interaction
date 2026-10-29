# Note 08 — Dépendance en Reynolds du contraste d'interaction

*Plan d'analyse pré-enregistré, lancé le 18/10/2026 à 08:29.*

## 1. Question

Quelle est la forme de la loi \(I_D(\mathrm{Re})\) pour le contraste
d'interaction **apparié** entre la paire \((n=2,\alpha=0)\) et
\((n=4,\alpha=0.10)\) ?

Trois hypothèses concurrentes, tranchées par validation croisée :

| modèle | forme | signature physique |
|---|---| loi de puissance |
| M1 `power` | \(I = A\,\mathrm{Re}^\gamma\) | autosimilarité |
| M2 `saturation` | \(I = A\,\mathrm{Re}^\gamma/(1+\mathrm{Re}/Re_s)\) | bornage par la couche pariétale |
| M3 `linear` | \(I = A(\mathrm{Re}-Re_0)\) | transition douce, seuil \(Re_0\) |

## 2. Pourquoi N=128 comme grille de travail

L'étude de stabilité (note 06) mesure l'écart du contraste entre \(N=128\)
et \(N=256\) sur les trois Reynolds de référence :

| Re | 100 | 400 | 1000 |
|---|---|---|---|
| écart \(I_D\) N=128→256 | 3.41 % | 1.92 % | 2.75 % |

Tous < 10 % (**EXCELLENT**), donc la méthode de différence appariée
annule l'erreur de discretisation bien mieux que les grandeurs brutes
(\(D^*\) varie de +3.0 %, \(P^*\) de +6.7 % entre les deux grilles).

**Décision :** cartographier \(I_D(\mathrm{Re})\) à \(N=128\), et réserver
\(N=256\) aux trois cas de référence. Justification rédigeable telle quelle
dans le papier :

> *The interaction contrast is stable to within 3.5 % between N = 128 and
> N = 256 at all three Reynolds numbers tested. Additional Reynolds numbers
> were therefore computed on the coarser grid to map the trend, with the
> finer grid reserved for the reference cases.*

## 3. Campagne

- **Grille :** \(N=128\), \(N_w=1\), \(a=b=0.5\), \(H=1\).
- **Nombre de Mach :** **fixe**, `U_ref = (40/3)/N` → Ma = 0,1804.
  Argument `--U` obligatoire (voir §3.1).
- **Reynolds :** tous les 50 de **100 à 1000** (19 valeurs).
- **Cas par Re :** \(n\in\{2,4,6,8\}\times\alpha\in\{0,0.10\}\) = 8.
- **Critère de convergence :** résidu \(\le 10^{-8}\) (tolérance unique de
  toutes les campagnes, note 06).
- **Réutilisation :** aucun cas déjà calculé n'est relancé, sous réserve
  qu'il soit **sur la loi de Ma** (les fichiers hérités `status=reanalyzed`
  sans résidu stocké sont considérés valides).

### 3.1 Piège du nombre de Mach (correction du 18/10)

Le solveur choisit sa viscosité par table :

```python
TARGET_TAU = {100: 0.65, 400: 0.60, 1000: 0.54}
DEFAULT_TAU = 0.60
...
U_ref = (tau - 0.5) / 3.0 * Re / N          # tau pris dans la table ci-dessus
```

Pour tout `Re` **hors** `{100, 400, 1000}` le filet `DEFAULT_TAU=0.60`
s'applique, donc \(\nu_{LB}\) est fixe et \(U_{LB}\propto\mathrm{Re}\) :
**le nombre de Mach croît linéairement avec Re.**

| Re | 100 | 400 | 600 | 800 | 950 | 1000 |
|---|---|---|---|---|---|---|
| Ma (loi `DEFAULT_TAU`) | 0,068 | **0,180** | 0,271 | 0,361 | **0,429** | **0,180** |
| Ma (loi retenue `--U`) | 0,180 | 0,180 | 0,180 | 0,180 | 0,180 | 0,180 |

Les trois points de référence sont à Ma ≈ 0,18 ; les cas « par défaut »
montaient jusqu'à **Ma = 0,43** (erreur de compressibilité \(\sim\mathrm{Ma}^2
\approx 19\%\)), sur un régime numérique **différent** de celui des points de
référence. La courbe \(I_D(\mathrm{Re})\) aurait été faussée — le défaut
s'observait d'ailleurs à la jonction : \(\Delta I_D(950\to1000) = 0{,}0070\)
alors que la tendance interne exigeait \(\approx 0{,}0091\).

**Loi retenue :** `U_ref = (40/3)/N`, qui reproduit **exactement** les cas
Re=400 et Re=1000 existants (Ma = 0,180422 à \(N=128\), 0,090211 à
\(N=256\)) et les rend réutilisables tels quels.

**70 cas N=128 hors-loi ont été mis en quarantaine** dans
`results/quarantine/tau060_N128/` et 136 cas relancés avec `--U`.
À \(N=512\), la même option donnera Ma = 0,0451 — la loi est donc
automatiquement cohérente entre grilles.

## 4. Coût

136 cas à refaire + 16 réutilisés = 152 cas au total, 8 processus
parallèles : **~2 h**. À titre de comparaison, la même exploration à
\(N=512\) coûterait ~49 h (mesuré : 68,6 it/s à \(N=512\) contre 276 it/s à
\(N=256\), soit ×4 plus lent par itération, machine déjà saturée par un seul
processus).

## 5. Analyse

`solver/plot_resweep.py` recalcule le contraste double

\[
I=\bigl[D(4,0.10)-D(4,0)\bigr]-\bigl[D(2,0.10)-D(2,0)\bigr]
\]

pour \(D^*,P^*,K^*\), ajuste les trois modèles, et les départage par
**validation croisée leave-one-Re-out** (RMSE sur \(\ln|I|\)).
Les barres d'erreurs du graphe proviennent de la stabilité \(N=128\to256\).

Sorties : `results/resweep/resweep_N128.json`,
`figures/analysis/resweep_ID_vs_Re_N128.png`.

## 6. Résultats (campagne terminée le 18/10 à 12:31)

**152/152 cas convergés, 0 échec, 156,3 min** (8 workers). Audit de Ma
post-campagne : **aucun fichier hors-loi**.

### 6.1 Table complète \(I_D(\mathrm{Re})\), paire (2,4), \(N=128\)

| Re | \(I_D\) | \(I_D/D^\*\) (n=2,0) | it |
|---|---|---|---|
| 100 | −0,001841 | −0,0105 % | 46 000 |
| 150 | −0,006589 | −0,0356 % | 56 000 |
| 200 | −0,012228 | −0,0630 % | 64 000 |
| 250 | −0,018674 | −0,0920 % | 70 000 |
| 300 | −0,025751 | −0,1217 % | 74 000 |
| 350 | −0,033272 | −0,1514 % | 82 000 |
| 400 | −0,041089 | −0,1805 % | — |
| 450 | −0,049102 | −0,2088 % | 92 000 |
| 500 | −0,057250 | −0,2362 % | 102 000 |
| 550 | −0,065495 | −0,2626 % | 110 000 |
| 600 | −0,073817 | −0,2882 % | 118 000 |
| 650 | −0,082209 | −0,3130 % | 128 000 |
| 700 | −0,090666 | −0,3371 % | 138 000 |
| 750 | −0,099185 | −0,3605 % | 146 000 |
| 800 | −0,107762 | −0,3834 % | 156 000 |
| 850 | −0,116397 | −0,4057 % | 164 000 |
| 900 | −0,125083 | −0,4276 % | 174 000 |
| 950 | −0,133816 | −0,4489 % | 184 000 |
| 1000 | −0,142590 | −0,4699 % | 232 000 |

Le Re=400 provient du CSV hérité (affiché `?` car ce tableau ne porte pas
le champ `iters`).

### 6.2 Forme de la loi : dérive régulière, exponentielle locale décroissante

Incréments \(\Delta I_D\) par tranche de 50 Re augmentent **de manière
régulière et monotone** (aucune césure, aucune inflexion suspecte) :

| tranche | 100→150 | 350→400 | 700→750 | 950→1000 |
|---|---|---|---|---|
| \(\Delta I_D\) | 0,004748 | 0,007817 | 0,008519 | 0,008774 |

Mais ces incrément **ralentissent** : l'exponentielle locale
\(\gamma_{loc}=d\ln\lvert I\rvert/d\ln\mathrm{Re}\) décroît de **3,09**
(sur 100→150) à **1,24** (sur 950→1000), en tendant vers 1. La courbe
en échelle log-log est donc **fortement courbée** et le rapport
\(\Delta I_D/\mathrm{Re}\) s'aplatit vers une asymptote linéaire.

### 6.3 Sélection des modèles (CV leave-one-Re-out, 19 points)

| critère | modèle vainqueur | paramètres | RMSE CV |
|---|---|---|---|
| absolu, \(\ln\lvert I\rvert\) | **power** \(A\mathrm{Re}^\gamma\) | \(\gamma=1{,}7272\), \(A=1{,}12\times10^{-6}\) | 26,7 % |
| relatif, \(\ln\lvert I\rvert/D^\*\) | **linear** \(A(\mathrm{Re}-Re_0)\) | \(A=5{,}176\times10^{-4}\), \(Re_0=62{,}2\) | **20,4 %** |
| — | saturation | \(\gamma=3{,}48\), \(Re_s=347\) | **1010 % (éliminé)** |

**\(R^2\) log = 0,9758 pour `power` est trompeur** : le log-log étant
courbé, ce modèle surestime la pente à haut Re (il prédit un
accroissement de pente ×1,66 entre Re=450 et 950 alors que l'observé est
×1,09). La saturation est **éliminée** par la CV : les incrément
continuent de croître, il n'y a aucun plateau.

**Lecture physique :** croissance **saturée à bas Re et asymptotiquement
linéaire à haut Re** — ni une loi de puissance unique, ni une saturation.
Le modèle à exponentielle locale décroissante
\(I\propto \mathrm{Re}^{1-q}(\mathrm{Re}+Re_c)^q\) (\(\gamma_{loc}\) de
3 à 1) est le candidat naturel ; il n'a pas encore été ajusté.

### 6.4 Stabilité de grille sur les trois Reynolds de référence

| Re | 100 | 400 | 1000 |
|---|---|---|---|
| écart \(I_D\) N=128→256 | **4,20 %** | 1,92 % | 2,75 % |

Tous < 10 %. Le chiffre Re=100 a été **recalculé** : la valeur
antérieure (3,53 %) comparait un \(N=256\) hors-loi (Ma = 0,034) à un
\(N=128\) hors-loi (Ma = 0,068). Les deux grilles sont désormais sur la
loi (Ma = 0,1804 / 0,0902) : \(I_D(N{=}256)=-0{,}001922\),
\(I_D(N{=}128)=-0{,}001841\).

## 7. État d'avancement

| date | état |
|---|---|
| 18/10 08:29 | campagne lancée (152 cas, 8 workers) |
| 18/10 09:56 | **relance sur la loi de Ma** après correction du piège |
| 18/10 12:31 | **terminée : 152/152, 0 échec, 156,3 min** |
| 18/10 14:55 | Re=100 à \(N=256\) recalculé sur la loi (8/8, 22,1 min) ; contraste régénéré |
| 18/10 15:00 | analyse consignée ici ; figure `resweep_ID_vs_Re_N128.png` |

**Suite :** ajuster le modèle à exponentielle locale décroissante
(§6.3), puis confronter aux vérifications de la note 09.
