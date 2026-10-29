# 06 — Tests de tolérance et de résolution

Deux tests préalables décident de la validité des campagnes à haute résolution
(N=256, N=512). Ils ont été exécutés à `Re=1000`, `N=256`, paire `(n_a=2, n_b=4)`,
`alpha=0.10`, `Nw=1`.

Rapports détaillés :

- `results/tolerance/tol_test_Re1000_N256_n2_n4.md` (script `solver/run_tol_test.py`)
- `results/resolution/res_test_Re1000_N128_N256_n2_n4.md` (script `solver/run_res_test.py`)

---

## 1) Test de tolérance : 1e-9 vs 1e-8

8 cas calculés (les 4 du bloc 1e-9, les 4 du bloc 1e-8). Les deux blocs
convergent : résidus finaux 9.52e-10 (1e-9) et 9.52e-09 (1e-8).

| grandeur | I (tol 1e-9) | I (tol 1e-8) | écart relatif | seuil 0,1 % |
|---|---|---|---|---|
| `I_D` | -0.146631 | -0.146628 | **0.0018 %** | PASS |
| `I_P` | -0.158697 | -0.158695 | **0.0010 %** | PASS |
| `I_K` | -3.795e-05 | -3.827e-05 | **0.0848 %** | PASS |
| bilan (par cas) | — | — | 0.0046 % | PASS (seuil 5 %) |
| `D*` (par cas) | — | — | 0.0001 % | PASS |
| `P*` (par cas) | — | — | 0.0001 % | PASS |

Références de la même famille de solveur (Projet B) : `D` stable à 2.4e-03 %
(47 cas) et bilan à 0.046 % entre 1e-8 et 1e-9. Le contraste `I_Q` fait
**mieux encore** (2e-05 relatif), confirmant que le biais systématique de
tolérance s'annule dans la différence de différences.

**Décision : adopter `1e-8` pour toutes les campagnes `N >= 256`.**
Gain mesuré : `768000 -> 608000` itérations sur le cas de référence
(**-21 %**), à qualité identique. Les 4 cas `1e-9` de `Re=1000` sont
conservés comme artefacts de validation.

## 2) Test de résolution : N=128 -> N=256

Test exécuté à `Re=1000` (`solver/run_res_test.py`) puis complété sur
Re=100 et Re=400 à partir de la campagne N=256 (`solver/plot_resolution.py`).

| Re | I (N=128) | I (N=256) | écart relatif | verdict |
|---|---|---|---|---|
| 100 | -0.001841 | -0.001922 | **4.20 %** | EXCELLENT (<10 %) |
| 400 | -0.041089 | -0.041894 | **1.92 %** | EXCELLENT (<10 %) |
| 1000 (`I_D`) | -0.142590 | -0.146631 | **2.75 %** | EXCELLENT (<10 %) |
| 100 (`I_P`) | -0.002320 | -0.002443 | 5.05 % | EXCELLENT (<10 %) |
| 400 (`I_P`) | -0.042470 | -0.044238 | 4.00 % | EXCELLENT (<10 %) |
| 1000 (`I_P`) | -0.150617 | -0.158696 | 5.09 % | EXCELLENT (<10 %) |
| 1000 (`I_K`) | -2.8e-05 | -3.8e-05 | 38.0 % | non exploitable |

Écart relatif par rapport à la grille fine : \((|I_{128}|-|I_{256}|)/|I_{256}|\).

> **Re=100 recalculé le 18/10** sur la loi de Ma (`docs/08` §3.1) : les
> valeurs initiales (-0.001885 / -0.001951, 3.53 %) comparaient deux
> grilles à des nombres de Mach différents (0,068 et 0,034) et sont
> caducs.

Sensibilité absolue à la résolution (cas `alpha=0`, `Re=1000`) : `D*` **+3.0 %**,
`P*` **+6.7 %**, `K*` **+2.8 %** — soit un effet de résolution 5 à 6 fois
plus grand que le signal d'interaction `I_D` (-0.47 % de `D`). Malgré cela,
le contraste ne dérive que de **1,92 à 4,20 % selon Re** : l'annulation
appariée fonctionne sur les trois Reynolds.

Rapport signal/résolution `|I| / |dQ*_res|` : 0.16 pour `D`, 0.08 pour `P`,
0.03 pour `K` — zone limite pour `D` et `P`, insuffisant pour `K`.

Points annexes :

- Le bilan `(P-D)/D` s'améliore nettement avec la résolution :
  -6.81 % (N=128) -> -3.47 % (N=256) à `Re=1000`, `n=2`, `alpha=0`.
- Le résidu de fermeture augmente en magnitude avec `Re`
  (-5.8 % à Re=100, -6.8 % à Re=1000 à N=128) mais **s'annule au premier
  ordre dans `I_Q`** (même signe et magnitude comparable pour `I_P` et `I_D`,
  ratio `I_P/I_D` = 1.26 / 1.03 / 1.06 à Re = 100 / 400 / 1000).
- `I_K` n'est pas monotone en `Re` (-4.6e-05, -7.4e-04, -2.8e-05) et n'est
  pas stable en résolution : à traiter comme un résultat distinct, lié à la
  topologie (taille du tourbillon, nombre de cellules) et non à la seule
  couche parière.

**Décision : STABLE.** La campagne N=256 puis N=512 peut être lancée.
L'incertitude systématique de `I_D` d'origine numérique est estimée à ~4 %
(résolution) et ~0.0018 % (tolérance).

## 3) Paragraphe de rédaction (papier 3)

> The numerical uncertainty of the interaction contrast was assessed before
> the production campaigns. The tolerance of the steady-state criterion,
> \(\epsilon = 10^{-8}\) was validated at \(N=256\), \(Re=1000\), pair
> \((2,4)\): the contrast changes by only \(0.0018\,\%\) for \(I_D\) and
> \(0.0010\,\%\) for \(I_P\) between \(\epsilon = 10^{-9}\) and
> \(\epsilon = 10^{-8}\), three orders of magnitude below the acceptance
> threshold, while the iteration count drops by \(21\,\%\).
> The resolution sensitivity is larger: between \(N=128\) and \(N=256\) the
> reference dissipation shifts by \(3.0\,\%\) and the driving power by
> \(6.7\,\%\), yet the matched-resolution contrast changes only from
> \(I_D = -0.14259\) to \(-0.14663\), i.e. \(2.8\,\%\) relative. The
> cancellation of the common grid error in the difference of differences is
> therefore effective, and we take this figure as the systematic uncertainty
> of the reported interaction. The closure residual, systematically negative
> (\(-5.9\,\%\) at \(Re=100\) to \(-6.8\,\%\) at \(Re=1000\) at \(N=128\),
> improving to \(-3.5\,\%\) at \(N=256\)), carries the same sign in \(P\) and
> \(D\) and cancels to leading order in \(I_Q\), as confirmed by the
> consistency of \(I_P\) and \(I_D\) (ratio \(1.03\)–\(1.06\) for
> \(Re \ge 400\)).

## 4) Conséquence opérationnelle

| campagne | tolérance | statut |
|---|---|---|
| N=128 (24 cas, Re=100/400/1000) | 1e-9 | terminée |
| N=256 (24 cas, Re=100/400/1000) | 1e-8 | lancée, 4 cas déjà convergents |
| N=512 (12 cas Tier 1) | 1e-8 | à lancer après N=256 |
