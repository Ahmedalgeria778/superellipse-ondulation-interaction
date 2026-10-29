# 07 — Résultats de la campagne N=256

## 1) État de la campagne

24 cas : `n ∈ {2,4,6,8}` × `alpha ∈ {0, 0.10}` × `Re ∈ {100, 400, 1000}`,
`Nw=1`, `N=256`, tolérance **1e-8** (validée par `docs/06_tests_tolerance_resolution.md`).

- **24/24 converge** ; 0 échec ; campagne exécutée en 223,9 min (6 workers,
  pool ordonné par coût décroissant : Re=1000 d'abord, puis 400, puis 100).
- Durées par cas : Re=1000 → 104–139 min ; Re=400 → 30–73 min ; Re=100 → 30–62 min.
- Itérations : Re=1000 → 594 000–622 000 ; Re=400 → 264 000–272 000 ;
  Re=100 → 242 000–264 000.
- Artéfacts : `results/scalars/pilot_Re{100,400,1000}_N256_tol1e-08.csv`,
  `results/scalars/interaction_Re{100,400,1000}_N256_tol1e-08.json`,
  `results/logs/*.log`.

## 2) Magnitudes intégrées (N=256)

| Re | cas | D* | P* | K* | (P−D)/D |
|---|---|---|---|---|---|
| 100 | n=2 α=0 | 17.86267 | 17.36857 | 0.033952 | −2.77 % |
| 100 | n=2 α=.1 | 17.85991 | 17.36638 | 0.033588 | −2.77 % |
| 100 | n=8 α=0 | 17.86720 | 17.37276 | 0.034190 | −2.77 % |
| 100 | n=8 α=.1 | 17.86207 | 17.36773 | 0.033733 | −2.77 % |
| 400 | n=2 α=0 | 23.24971 | 22.50600 | 0.039792 | −3.20 % |
| 400 | n=2 α=.1 | 23.12922 | 22.38916 | 0.039109 | −3.20 % |
| 400 | n=8 α=0 | 23.32820 | 22.58278 | 0.040114 | −3.20 % |
| 400 | n=8 α=.1 | 23.16136 | 22.41778 | 0.039278 | −3.21 % |
| 1000 | n=2 α=0 | 31.25336 | 30.16895 | 0.042944 | −3.47 % |
| 1000 | n=2 α=.1 | 30.89745 | 29.82759 | 0.042100 | −3.46 % |
| 1000 | n=8 α=0 | 31.43923 | 30.35194 | 0.043220 | −3.46 % |
| 1000 | n=8 α=.1 | 30.92489 | 29.84221 | 0.042352 | −3.50 % |

> **Re=100 recalculé le 18/10** sur la loi de Ma (voir `docs/08` §3.1) :
> les 8 cas initiaux étaient hors-loi (Ma = 0,034 au lieu de 0,090) et ont
> été mis en quarantaine. Les valeurs ci-dessus sont les neuves.

**Bilan de fermeture** : −2.77 % (Re=100) → −3.20 % (Re=400) → −3.47 %
(Re=1000) à N=256, contre −5.90 % / −6.20 % / −6.81 % à N=128.
La résolution améliore le résidu de **~50 %** ; la dépendance en Re subsiste
mais s'atténue.

## 3) Interactions `I_Q` à N=256

Contraste `I_Q = [Q(nb,α) − Q(nb,0)] − [Q(na,α) − Q(na,0)]`, α = 0.10.

| Re | paire | `I_D` | `I_P` | `I_K` |
|---|---|---|---|---|
| 100 | (2,4) | −0.001922 | −0.002443 | −0.87e−04 |
| 100 | (2,6) | −0.002303 | −0.002782 | −0.93e−04 |
| 100 | (2,8) | −0.002369 | −0.002842 | −0.93e−04 |
| 100 | (4,6) | −0.000381 | −0.000338 | −0.06e−04 |
| 100 | (4,8) | −0.000447 | −0.000399 | −0.07e−04 |
| 100 | (6,8) | −0.000066 | −0.000060 | −0.01e−04 |
| 400 | (2,4) | −0.041894 | −0.044238 | −1.53e−04 |
| 400 | (2,6) | −0.045954 | −0.047812 | −1.53e−04 |
| 400 | (2,8) | −0.046345 | −0.048149 | −1.53e−04 |
| 400 | (4,6) | −0.004060 | −0.003575 | ~ 0 |
| 400 | (4,8) | −0.004451 | −0.003912 | ~ 0 |
| 400 | (6,8) | −0.000391 | −0.000337 | ~ 0 |
| 1000 | (2,4) | −0.146628 | −0.158694 | −0.38e−04 |
| 1000 | (2,6) | −0.158131 | −0.168123 | −0.22e−04 |
| 1000 | (2,8) | −0.158424 | −0.168376 | −0.24e−04 |
| 1000 | (4,6) | −0.011503 | −0.009429 | +0.17e−04 |
| 1000 | (4,8) | −0.011795 | −0.009682 | +0.15e−04 |
| 1000 | (6,8) | −0.000293 | −0.000253 | −0.02e−04 |

**Structure confirmée et identique à N=128 :**

- `I_Q < 0` partout : l'ondulation augmente D et P, et **plus fortement
  pour la forme n=2 que pour n ≥ 4** — l'interaction est donc négative.
- Fort contraste 2 → 4, puis **saturation** : `I(4,6)`, `I(4,8)`, `I(6,8)`
  sont 10 à 30 fois plus petits que `I(2,4)`.
- Signal relatif `I_D/D(n=2,α=0)` : −0.011 % (Re=100), −0.180 % (Re=400),
  −0.469 % (Re=1000) à N=256, contre −0.011 / −0.181 / −0.470 % à N=128 —
  **quasi identique en relatif**.
- Cohérence `I_P/I_D` : **1.27 (Re=100), 1.06 (Re=400), 1.08 (Re=1000)**
  (valeurs N=128 : 1.26 / 1.03 / 1.06). Puissance et dissipation décrivent
  le même effet.
- `I_K` reste petit, non monotone en Re, et change de signe pour les paires
  (4,6)/(4,8) à Re=1000 : quantité de nature topologique, à traiter séparément.

## 4) Stabilité en résolution (le test décisif)

Écart relatif défini par rapport à **la grille fine de référence** :
\(\varepsilon=(\lvert I_{128}\rvert-\lvert I_{256}\rvert)/\lvert I_{256}\rvert\).

| Re | `I_D` N=128 | `I_D` N=256 | écart relatif | seuil 10 % |
|---|---|---|---|---|
| 100 | −0.001841 | −0.001922 | **4.20 %** | OK |
| 400 | −0.041089 | −0.041894 | **1.92 %** | OK |
| 1000 | −0.142590 | −0.146628 | **2.75 %** | OK |
| 100 (`I_P`) | −0.002320 | −0.002443 | 5.05 % | OK |
| 400 (`I_P`) | −0.042470 | −0.044238 | 4.00 % | OK |
| 1000 (`I_P`) | −0.150617 | −0.158694 | 5.09 % | OK |

> Les valeurs Re=100 sont **postérieures au recalcul sur la loi de Ma**
> (`docs/08` §3.1). Les chiffres antérieurs (3.53 % et 6.6 %)
> comparaient deux grilles à des nombres de Mach différents et sont
> caducs. La convention d'écart est ici unique et identique à `docs/06`.

Alors que `D*` bouge de **+2.99 %** et `P*` de **+6.69 %** entre les deux
grilles, le contraste ne dérive que de **1,92 à 4,20 %**. L'annulation
appariée de l'erreur de grille est donc efficace, sur les trois Reynolds.

**Incertitude numérique retenue sur `I_D` : ~4 % (résolution) ; 0,0018 % (tolérance).**

## 5) Figures

- `figures/analysis/resolution_effect_vs_n.png` — δQ(n)/Q(n,0) aux deux grilles.
- `figures/analysis/resolution_IQ_stability.png` — `|I_Q|` vs Re, changement
  relatif N=128→256 (seuils 10 % / 30 %) et ratio `I_P/I_D`.
- `figures/analysis/effect_vs_n_Re400.png`, `interaction_Re400.png` (N=128).

Régénérer : `python -m solver.plot_resolution`.

## 6) Verdict et suite

| critère | mesure | verdict |
|---|---|---|
| tolérance 1e-8 vs 1e-9 sur `I_Q` | 0.0018 % | ADOPTER 1e-8 |
| stabilité `I_D` N=128→256 | 1.92–4.20 % | STABLE (seuil 10 %) |
| convergence de la campagne N=256 | 24/24 | OK |

Le papier est viable en l'état : le signal d'interaction (−0.011 % à −0.47 %
selon Re) est mesurable avec une incertitude numérique ~4 %.

**Suite : campagne N=512 (Tier 1, 12 cas)** — Re=1000 : n ∈ {2,4,6,8} ×
α ∈ {0,0.10} (8 cas) ; Re=400 : n ∈ {2,4} × α ∈ {0,0.10} (4 cas) ;
tolérance 1e-8, 4 workers, checkpoint actif.

## 7) Paragraphe de rédaction (papier 3)

> The \(N=256\) campaign (24 configurations) was run to completion with the
> adopted tolerance \(\epsilon = 10^{-8}\). The interaction remains negative
> for every pair and every Reynolds number, with the same structure as at
> \(N=128\): a large contrast for the pair \((2,4)\) followed by a saturation
> for \(n \ge 4\). Expressed relative to the reference dissipation, the
> signal is \(-0.011\,\%\) at \(Re=100\), \(-0.180\,\%\) at \(Re=400\) and
> \(-0.469\,\%\) at \(Re=1000\), against \(-0.011\,\%\), \(-0.181\,\%\) and
> \(-0.470\,\%\) at \(N=128\). The matched-resolution contrast is therefore
> grid-independent to within 4.5 \% on the three Reynolds numbers studied,
> while the absolute dissipation itself shifts by 3 \%. Power and dissipation
> remain consistent, with \(I_P/I_D\) between 1.06 and 1.08 for \(Re \ge 400\).
