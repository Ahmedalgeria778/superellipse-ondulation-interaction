# 12 - Vérification du nombre de Mach et origine de la dérive du contraste à `Re=100`

**Date** : 2026-11-17 (version finale ; la première version du 2026-11-16 portait des chiffres provisoires)
**Statut** : terminé. Le run `n=4`, `N=768`, `Re=100` a convergé le 2026-11-17 à 08:49:26
(résidu `9.889e-9`), ce qui permet de trancher entre « dérive physique » et « artefact de convergence ».
**Données** : `results/scalars/`, `results/fields/`, `results/ma_sweep/`, `results/dAdB_multires_v2.json`,
`results/dAdB_N768_final.json`, `results/richardson_extrap_v2.json`, `results/mach_summary.json`.

---

## 1. Le nombre de Mach pose-t-il un problème ?

### 1.1 Valeurs de production

Le pilotage est `U = (tau - 1/2) Re / (3 N)` avec `tau` fixé par Reynolds
(`TARGET_TAU = {100: 0.65, 400: 0.60, 1000: 0.54}`), donc `Ma = sqrt(3) U`
décroît en `1/N` (scaling diffusif). Les valeurs enregistrées dans les JSON
sont :

| Re | N | U | tau | Ma | `drho/rho` mesuré |
|---:|---:|---:|---:|---:|---:|
| 100 | 128 | 0.039063 | 0.65 | 0.0677 | 2.36e-2 |
| 100 | 256 | 0.019531 | 0.65 | 0.0338 | — |
| 100 | 384 | 0.013021 | 0.65 | 0.0226 | 2.89e-3 |
| 100 | 768 | 0.006510 | 0.65 | 0.0113 | 7.56e-4 |
| 400 | 128 | 0.104167 | 0.60 | 0.1804 | — |
| 400 | 384 | 0.034722 | 0.60 | 0.0601 | — |
| 400 | 512 | 0.026042 | 0.60 | 0.0451 | 3.82e-3 |
| 1000 | 128 | 0.104167 | 0.54 | 0.1804 | 3.73e-2 |
| 1000 | 384 | 0.034722 | 0.54 | 0.0601 | 4.27e-3 |

Sur les grilles de production (`N >= 384`) : `Ma <= 0.0601` et
`drho/rho <= 4.3e-3` (0.43 %). Au `N=768`, `Ma = 0.0113` et
`drho/rho = 7.6e-4`. La compressibilité est donc faible là où les résultats
sont cités. Elle n'est **pas** négligeable sur les grilles grossières
(`N=128` : `Ma = 0.1804`, `drho/rho` jusqu'à 3.7e-2 pour `n=2`), ce qui
justifie de n'utiliser ces grilles que pour le test de résolution.

Les cas `n>=6` à `Re=100`, `N=128` ont été lancés avec `U = 5/48`
(`tau = 0.9`, `Re=100` exact) et atteignent `Ma = 0.1804` avec
`drho/rho = 1.475e-1`. Ils ne sont pas comparables aux cas `n=4` à
`Ma = 0.0677` ; c'est déjà écrit dans le SM (Table SM.X, note).

### 1.2 Test de sensibilité en Mach à `Re`, `N` fixés

Un balayage dédié existe dans `results/ma_sweep/` : à `Re` et `N` fixés, on
fait varier `U` d'un facteur 8 et `tau = 1/2 + 3 U N / Re` en conséquence,
en recalculant **les deux** membres `n=2` et `n=4`. Le contraste
`I_D(2,4)` est donc une différence appariée au même `(U, tau)`.

Résultats (seuls les cas `status=converged` sont utilisés) :

| Re | N | plage de Ma | variation de `D*` | plage de `I_D` (pp) | variation relative |
|---:|---:|---:|---:|---:|---:|
| 100 | 128 | 0.0226 -> 0.1804 | 1.29 % | +0.021502 -> +0.021804 | 1.39 % |
| 100 | 256 | 0.0226 -> 0.1804 | 1.39 % | +0.021877 -> +0.021967 | 0.41 % |
| 400 | 128 | 0.0451 -> 0.3608 | 0.60 % | +0.321025 -> +0.322940 | 0.59 % |
| 400 | 256 | 0.0451 -> 0.3608 | 0.47 % | +0.322701 -> +0.324650 | 0.60 % |
| 1000 | 128 | 0.0601 -> 0.4811 | 0.69 % | +0.576074 -> +0.584966 | 1.54 % |
| 1000 | 256 | 0.1203 -> 0.4811 | 0.23 % | +0.587871 -> +0.596393 | 1.44 % |

Lecture :

- l'intégrale absolue `D*` varie de 0.2 à 1.4 % sur un facteur 8 de `Ma` :
  elle porte bien une erreur de compressibilité sur les grilles grossières ;
- le **contraste apparié** `I_D(2,4)` ne varie que de 0.4 à 1.5 % de sa
  valeur sur le même facteur 8. L'erreur commune se cancelle dans la
  différence ;
- les valeurs du balayage recoupent la production : à `Re=100`, `N=256`,
  `tau=0.65` (`Ma=0.0338`), la production donne `+0.02195 %`, et
  l'interpolation entre `Ma=0.0226` (`+0.021922`) et `Ma=0.0451`
  (`+0.021959`) donne `+0.02194 %`.

### 1.3 Le Mach peut-il expliquer le déplacement `N=384 -> 768` ?

Non. Trois arguments, tous vérifiés sur les données finales :

- à `Re=100`, `I_D` **augmente** avec `Ma` (`+0.021502` à `Ma=0.0226`,
  `+0.021804` à `Ma=0.1804`) ;
- or le raffinement `N=384 -> 768` **diminue** `Ma` (0.0226 -> 0.0113) ;
  l'effet Mach attendu est donc une très légère **baisse** de `I_D`,
  d'au plus quelques `1e-5` pp, alors que l'on observe une **hausse**
  de `5.39e-4` pp ;
- le déplacement observé est donc de signe opposé et plusieurs fois plus
  grand que ce que le Mach peut produire, et il est surtout
  **interprété par la partition `dA/dB`** (section 2), qui montre que le
  mouvement est un résidu de cancellation et non un effet de
  compressibilité.

Rappel de cohérence : le run `N=768` a le **Mach le plus faible** de toute la
campagne (`Ma = 0.0113`), donc s'il y avait un effet Mach, il irait dans le
sens d'un résultat *plus* fiable, pas moins.

**Conclusion Mach** : le nombre de Mach ne menace ni le signe, ni l'ordre de
grandeur, ni la comparaison appariée `I_D(2,4)`. Il explique une partie de
l'erreur sur les intégrales absolues aux grilles `N <= 256`, et c'est
pourquoi seuls les contrastes appariés à `(Re, N, tau)` fixés sont utilisés
pour les conclusions. Le Mach n'est pas la cause du déplacement
`N=384 -> 768`.

---

## 2. Origine du déplacement du contraste à `Re=100`

### 2.1 Les chiffres (tous convergés)

| N | `D*` n=2 | `D*` n=4 | `I_D(2,4)` |
|---:|---:|---:|---:|
| 128 | 17.417928 | 17.421691 | +0.021605 % |
| 256 | 17.774193 | 17.778090 | +0.021924 % |
| 384 | 17.915213 | 17.919115 | +0.021785 % |
| 768 | 18.070223 | 18.074257 | **+0.022323 %** |

`n=2`, `N=768` : convergé, résidu `9.713e-9`, `1 172 000` itérations en un
seul processus (le journal montre la rampe de couvercle complète de 768 pas,
donc aucun redémarrage).
`n=4`, `N=768` : convergé, résidu `9.889e-9`, avec un redémarrage depuis un
checkpoint. Le compteur d'itérations du solveur repart de zéro à chaque
reprise (`ckpt.activate` remet `it = 0` et `load` ne restaure pas le
compteur) : seul le dernier segment est donc enregistré, `378 000`
itérations (`iters_last_process` dans le JSON). Le checkpoint consommé avait
été écrit `400 000` pas après le démarrage de son propre processus, et un
journal antérieur (`.before_restart.log`) montre un processus plus ancien
allant jusqu'à `1.04e6` pas. Le **total cumulé n'est pas enregistré** et
n'est donc cité nulle part.

La variation `N=384 -> 768` vaut :

- en absolu : `+5.39e-4` point de pourcentage (pp) ;
- en relatif : `+2.5 %`.

Les deux intégrales, elles, changent de `+0.865 %` entre `N=384` et `N=768`
(`+0.8652 %` pour `n=2`, `+0.8658 %` pour `n=4`). Autrement dit, l'erreur de
grille commune se cancelle dans la différence à hauteur de ~99.94 %, et le
résidu de cancellation (`5.4e-4` pp, soit `3e-5` de `D*`) reste comparable au
contraste lui-même (`2.2e-5` de `D*`). C'est ce rapport qui rend le contraste
relativement sensible à la grille, tout en le laissant stable en absolu.

### 2.2 Partition exacte `dA/dB` aux deux maillages les plus fins

Calculée avec le code de production (`solver/referee_mech_analysis.py`,
`build_combined_geometry` + `diss_density`), sur les deux champs convergés.
Unités : `1e-4 D*(2,0)`. Source : `results/dAdB_multires_v2.json`, dont la
ligne `N=384` reproduit `partition_bands_N384.json` à `7e-18`.

| bande | N=384 | N=768 | changement |
|---|---:|---:|---:|
| lid | -1.332 | -1.614 | -0.282 |
| bottom | -1.289 | -1.282 | +0.007 |
| side | +7.132 | +7.173 | +0.041 |
| core | -9.325 | -9.006 | +0.319 |
| **dA** | **-4.814** | **-4.729** | **+0.085** |
| **dB** | **+6.992** | **+6.961** | **-0.031** |
| **net** | **+2.179** | **+2.232** | **+0.054** |
| `I_D` (%) | +0.021785 | +0.022323 | +0.000539 |
| `C` | 3.96 | 4.03 | +0.07 |

Deux faits :

1. **Tout est petit et stable** : sur quatre maillages, `dB` reste dans
   `6.89`-`7.01` et `dA` dans `-4.73`-`-4.82`, soit moins de 2 % de variation
   sur chaque terme, alors que leur somme se déplace de 2.5 %.
   La région ajoutée ne dissipe pas davantage en raffinant ; c'est la
   redistribution dans la région commune qui bouge légèrement.
2. `dA` est une somme de bandes qui s'opposent : `side = +7.2` contre
   `core = -9.0` et `lid+bottom = -2.9`. L'indice de cancellation
   `C = sum|bandes| / |dA|` vaut 4.0 sur les quatre grilles
   (`3.9` à `4.0`). C'est cette cancellation qui fait qu'une variation de
   quelques pour cent des bandes déplace `dA` de quelques pour cent et le
   contraste de quelques pour cent.

La séquence `dA` à `Re=100` est `-4.73`, `-4.82`, `-4.81`, `-4.73` : elle
revient à sa valeur `N=128`, sans tendance monotone. Il n'y a donc **pas** de
dérive systématique à expliquer ; le point `N=768` est simplement un
échantillon supplémentaire d'une quantité qui oscille dans une bande de
±2 %.

### 2.3 Le déplacement de +19 % était un artefact de convergence

C'est le point que la première version de cette note ne pouvait pas trancher.
Le checkpoint intermédiaire du champ `n=4` (environ 1.4 M itérations dans le
run, résidu `~6e-7`)
donnait `D*(n=4) = 18.074911` et `I_D = +0.025944 %`, soit `+19 %` au-dessus
de `N=384`. Le run convergé donne `D* = 18.074257` et `I_D = +0.022323 %`,
soit `+2.5 %`.

Entre les deux états, `D*` n'a changé que de `3.6e-5` en relatif, mais
`I_D` a changé de 14 % en relatif. Le facteur d'amplification est
exactement `1/I_D` : une erreur relative `eps` sur **un seul** des deux
champs déplace le contraste en relatif de `eps / 2.23e-4`, soit un facteur
`~4.5e3`. Pour connaître `I_D` à 1 %, il faut donc `eps < 2.2e-6`
(0.00022 %) sur chaque champ — c'est précisément la raison pour laquelle le
run a été poussé jusqu'à `residual < 1e-8` et non arrêté à un plafond
d'itérations. Un plafond à 400 k itérations (`*_CAPPED.json`, `D* = 18.0254`)
aurait donné un contraste faux de plus de 20 %.

**Aucun nombre issu d'un champ non convergé n'est cité comme résultat** dans
le manuscrit ni dans le SM.

### 2.4 Extrapolation de Richardson (indicative)

Avec les quatre niveaux (`N=128,256,384,768`, tous convergés) :

| série | `Q_inf` | `p` | écart au plus fin |
|---|---:|---:|---:|
| `n=2` | 18.302256 | 0.7482 | +1.284 % |
| `n=4` | 18.306518 | 0.7479 | +1.285 % |
| contraste extrapolé | — | — | **+0.023288 %** |

Stabilité au retrait d'un maillage :

| niveaux | `p` (n=2 / n=4) | contraste extrapolé |
|---|---:|---:|
| 128,256,384,768 | 0.748 / 0.748 | +0.0233 % |
| 128,256,384 | 0.699 / 0.700 | +0.0191 % |
| 256,384,768 | 0.820 / 0.818 | +0.0257 % |

Le `p` est presque identique pour `n=2` et `n=4` à chaque sous-ensemble : les
deux intégrales obéissent à la même loi de convergence, et le contraste est
le résidu d'une cancellation, pas une quantité autonome. Le contraste
extrapolé varie d'un facteur `1.35` selon le jeu de maillages : c'est
pourquoi **aucune valeur extrapolée n'est citée comme résultat**, et pourquoi
aucun GCI n'est calculé.

Un ajustement direct de `I_D(N)` (au lieu des deux intégrales) est dégénéré :
l'exposant part aux bornes imposées (`p = 0.2` sur 4 niveaux et sur
`256-768`, `p = 8.0` sur `128-384`) et le contraste extrapolé varie de
`+6 %` à `+7.7 %` au-dessus de la valeur la plus fine. La note de
`results/richardson_extrap_v3_provisional.json` documente cet échec ; il n'en
est rien tiré.

### 2.5 Conclusion

Le déplacement `N=384 -> 768` (+2.5 % en relatif, `+5.4e-4` pp en absolu)
n'est pas :

- un effet de Mach (section 1.3 : signe opposé, amplitude trop grande) ;
- un effet de Reynolds (Re est fixé à 100) ;
- une preuve que le contraste diverge (les quatre valeurs tiennent dans une
  bande de `7.2e-4` pp, soit 3.2 %).

C'est un **résidu de troncature amplifié par la cancellation dans `dA`** :
les intégrales absolues changent de 0.87 % entre `N=384` et `N=768`, leur
différence de 0.02 %, et le résidu de cette différence est du même ordre que
le contraste. Le point `N=768` est une **extension de résolution**, pas un
maillage convergé : le critère `< 0.2 %` entre les deux maillages les plus
fins n'est pas atteint (0.87 % sur `D*`) et n'est revendiqué nulle part.

---

## 3. Ce que cela change dans le manuscrit (état final)

- Les valeurs de production restent celles de `N=384` (`+0.021785 %` à
  `Re=100`), et le point `N=768` convergé est ajouté comme quatrième niveau :
  `+0.022323 %`.
- L'incertitude de grille à `Re=100` est le span sur les maillages
  disponibles : `7.2e-4` pp (3.2 %) sur `N=128-768`. Le span `N=128-384`
  (`3.4e-4` pp, 1.6 %) est cité séparément là où la distinction compte.
- Le texte ne dit plus « stable à 2.4 % » de façon globale : il dit
  « stable à quelques pour cent de sa valeur sur les maillages disponibles »,
  avec les trois spans explicites (7.2e-4, 2.4e-3, 1.4e-2 pp).
- Les conclusions (signe positif, ordre de grandeur, sous-extensivité,
  hiérarchie des transitions) ne changent pas : le contraste reste
  sous-extensif d'un facteur 13 à 364 par rapport au modèle de surface.
- Le déplacement de +19 % du checkpoint intermédiaire est documenté dans le
  SM (§ « High-resolution extension ») comme un artefact de convergence,
  précisément pour montrer pourquoi aucun résultat n'est tiré d'un champ non
  convergé.

---

## 4. État des calculs au 2026-11-17

### Terminés

- LBM `N=768`, `Re=100`, `n=2` et `n=4`, `alpha=0` : **convergés** à
  `tol=1e-8` (résidus `9.71e-9` et `9.89e-9`), sans plafond d'itérations.
- Partition `dA/dB` sur 12 lignes (`Re=100/400/1000` × `N=128...768`),
  validée à `7e-18` sur la ligne `N=384`.
- Richardson 4 niveaux et tests de retrait (`richardson_extrap_v2.json`).
- FEM Taylor-Hood `Re=1000`, `nx=64/128/256` (`n=2,4`) : terminés ;
  `nx=256` donne `D*(n=2)=32.20162`, `D*(n=4)=32.39251`,
  `I_D = +0.593 %`, bilans natifs `+0.0109 %` et `+0.0068 %`.
- FEM `Re=400`, `nx=32` : `n=2` et `n=4` terminés (2026-11-17, 44 et 48
  itérations, 23 s et 25 s ; résidus de point fixe < `1e-10`). Résultats :
  `D* = 25.554` (`n=2`) et `26.599` (`n=4`) sur le champ projeté,
  `25.625` et `25.691` en natif, soit `I_D(FEM, nx=32) = +4.1 %` contre
  `+0.32 %` en LBM à `N=384`. L'écart (facteur 13) confirme que `nx=32` est
  trop grossier pour un contraste : ces deux cas restent hors des tableaux
  et des figures, et seront versés à la prochaine archive de données
  (v1.2.0), pas à la v1.1.0 déjà publiée. Le diagnostic signale par ailleurs
  un Mach local de 1.67 sur ce maillage (cf. `Ma_ok = false` dans les JSON),
  contre 0.0601 en LBM : autre raison de ne pas s'en servir.
- Balayage Mach (`results/ma_sweep/`, tableaux SM `supptab:fitma` et
  `supptab:ma_low_high`).

### En cours

- LBM `N=768`, `Re=400`, `alpha=0`, `n=2` et `n=4` : **lancés le
  2026-11-17 à 11:23**, sans plafond (`tol=1e-8`, `max_iter=1e8`), sous le
  superviseur `solver/watch_long.py` qui relance automatiquement après une
  coupure et s'arrête sur `status = converged` (JSON). Débit mesuré :
  20 000 itérations en 1 075 s, soit `0.0537` s/it ; à ce rythme, 550 k à
  1.2 M itérations donnent **8 à 18 h** par cas, les deux tournant en
  parallèle. Journaux : `results/logs/watch_Re400_N768_n{2,4}.log` et
  `results/logs/se_n{2,4}_w0.00_Nw1_Re400_N768_tol1e-08.log` ; checkpoints
  toutes les 200 000 itérations dans `results/ckpt/`.
- LBM `N=768`, `Re=1000`, `alpha=0`, `n=2` et `n=4` : **lancés le
  2026-11-17 à 11:53**, dans la foulée des cas `Re=400`, avec le même
  protocole sans plafond (`tol=1e-8`, `max_iter=1e8`) et les mêmes
  superviseurs. Au 2026-11-17 11:56, la rampe de couvercle est terminée dans
  les deux journaux (`results/logs/se_n{2,4}_w0.00_Nw1_Re1000_N768_tol1e-08.log`)
  mais le premier point de mesure `it=20000` n'est pas encore atteint : le
  débit et l'ETA ne seront connus qu'après cette première ligne. Quatre
  superviseurs et leurs quatre processus solveurs sont actifs simultanément.
  Ce calcul étend le test de résolution au Reynolds le plus élevé, où le
  span est déjà le plus large (`1.4e-2` pp sur `N=128-384`).

### Non lancés (perspectives, explicitement étiquetés)

- Limite 2D à `Re=1000` : explicitée dans le manuscrit, non quantifiée par un
  calcul 3D.
- Crossover `Re ≈ 346` : exploratoire, présenté comme tel.
