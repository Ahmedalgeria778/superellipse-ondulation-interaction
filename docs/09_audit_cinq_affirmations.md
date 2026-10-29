# Note 09 — Audit des cinq affirmations du texte d'ambition

*Réponse à la demande d'audit : confronter les contrastes \(I_Q\), les
décompositions Wall/Core et les données de convergence aux cinq
affirmations, puis décider lesquelles peuvent être présentées comme des
conclusions démontrées.*

**Principe :** chaque affirmation est jugée sur les fichiers réellement
présents dans le projet (`results/`, `figures/`, `docs/`), jamais sur
l'intention. Trois verdicts possibles :

- **DÉMONTRÉ** — le résultat est mesuré, significatif vis-à-vis de
  l'incertitude, et reproductible à partir des artefacts ;
- **TESTABLE** — le calcul existe ou est prévu, mais n'a pas été fait ;
- **HYPOTHÈSE** — aucune donnée ne le soutient pour l'instant.

---

## 0. Définition des grandeurs \(Q\) (point explicitement demandé)

Source : `provenance/Lid_driven_ondulation/diagnostics.py:9-12`.

| symbole | grandeur | définition |
|---|---|---|
| \(D^\*\) | dissipation visqueuse | \(D_\nu=\int 2\mu\,S{:}S\,dA\), normalisée par \(LHU_0^2\) |
| \(P^\*\) | puissance d'entraînement de la lid | \(P_{lid}=\mu\int U_{lid}(\partial u/\partial y)_{lid}dx\), normalisée |
| \(K^\*\) | énergie cinétique | \(K/(LHU_0^2)\) |
| \((P-D)/D\) | résidu du bilan mécanique | doit tendre vers 0 à l'état stationnaire |

Le contraste est **toujours** une différence de différences appariée :

\[
I_Q(n_a,n_b)=\bigl[Q(n_b,\alpha)-Q(n_b,0)\bigr]-\bigl[Q(n_a,\alpha)-Q(n_a,0)\bigr]
\]

Le texte d'ambition devra nommer \(Q\) explicitement ; « énergétique »
seul est insuffisant devant un relecteur.

---

## 1. Non-additivité des effets géométriques → **DÉMONTRÉ**

### Ce qui existe

- **18 contrastes à \(N=256\)** (6 paires × 3 Re) dans
  `results/scalars/interaction_Re{100,400,1000}_N256_tol1e-08.json`,
  valeurs reprises dans `docs/07:43-62`.
- **19 contrastes à \(N=128\)** (une paire, 19 Re de 100 à 1000) dans
  `results/resweep/resweep_N128.json`, figure
  `figures/analysis/resweep_ID_vs_Re_N128.png`.
- **\(\mathbf{I_Q<0}\) sans exception** : 18/18 à \(N=256\), 19/19 à
  \(N=128\). L'ondulation augmente \(D\) et \(P\), et **plus fortement
  pour \(n=2\) que pour \(n\ge4\)** : les deux effets ne s'additionnent
  pas, ils se compense partiellement.

### Amplitude et significativité

| Re | \(\lvert I_D\rvert\) (N=128) | \(\lvert I_D\rvert\)/\(D^\*\) | incertitude résolution | rapport signal/bruit |
|---|---|---|---|---|
| 100 | 0,001841 | −0,011 % | 4,20 % | **~24** |
| 400 | 0,041089 | −0,180 % | 1,92 % | **~52** |
| 1000 | 0,142590 | −0,469 % | 2,75 % | **~36** |

Le point méthodologique à écrire noir sur blanc :

> le signal d'interaction est **10 à 50 fois plus petit** que l'erreur de
> grille sur \(D^\*\) lui-même (+3,0 %), mais **18 à 52 fois plus grand**
> que l'incertitude résiduelle sur le contraste apparié. C'est
> précisément l'appariement qui rend l'effet mesurable.

### Robustesse

- Même signe et même structure aux **deux grilles**.
- Deux grandeurs physiquement distinctes concordent :
  \(I_P/I_D = 1{,}03\) à \(1{,}24\) (valeurs N=256 : 1,24 / 1,06 / 1,08).
- Tolérance : passer de \(10^{-9}\) à \(10^{-8}\) ne déplace \(I_D\) que de
  **0,0018 %** (`docs/06`).

### Ce qui manque pour être irréprochable

- Incertitude estimée sur **3 Reynolds seulement**, par réplication de
  grille et non par réplication statistique.
- Stabilité à **2 grilles** : pas d'ordre apparent ni d'extrapolation
  Richardson sur \(I_Q\).
- **Révisé le 18/10 :** le cas Re=100 était hors-loi de Mach (Note 08
  §3.1) ; il a été recalculé sur les deux grilles. Stabilité finale :
  **4,20 %** pour `I_D` et 5,05 % pour `I_P`. Les chiffres antérieurs
  (3,53 % et 6,6 %) comparaient deux grilles à des nombres de Mach
  différents et sont caducs.

### Verdict

**DÉMONTRÉ.** La non-additivité est mesurée, du bon signe, significative
vis-à-vis de sa propre incertitude, et confirmée par deux grandeurs et
deux résolutions. C'est la contribution principale défendable.

---

## 2. Médiation par la couche pariétale → **HYPOTHÈSE NON TESTÉE**

### Ce qui existe : aucune décomposition Wall/Core

Recherche exhaustive (`solver/`, `provenance/`, `results/`, `figures/`) :

- **`solver/` : rien.** `checks.py:44` ne produit que \(D^\*,P^\*,K^\*\)
  intégrés sur **tout** le domaine. Aucun masque, aucune bande, aucun
  découpage pariétal.
- **`provenance/` : deux scripts non exploitables**
  - `06_physics/analyse_coeur_paroi.py` — bande \(d(x,\Gamma)\ge d_0\),
    \(d_0=0{,}05\) ; **non exécutable** : importe
    `etude_convergence_paroi`, fichier absent.
  - `06_physics/analyse_etape2_3_coeur.py` — écrit dans
    `results_coeur/`, **répertoire inexistant** : zéro sortie générée.
- **Le cahier des charges existe** (`docs/02:74-85` : seuil
  \(\delta_w=c\,h\), \(c\in[1,3]\), robustesse au choix de \(c\)) mais
  **n'a jamais été implémenté** ; `docs/05:104` le liste en « Suite ».
- Le `D_core` visible dans `audit/inventaire_cas.csv` est **hérité du
  dataset externe B**, pas calculé ici (`docs/00:44` prévient d'ailleurs
  que les deux définitions sont incompatibles).

### Conséquence

L'interprétation H2 (« l'ondulation mord sur la région proche-paroi »,
`docs/05:89-91`) **n'appuyée sur aucune donnée de ce projet**. Dans une
cavité fermée il faut de surcroît parler de *région pariétale* et non de
couche limite au sens classique.

### Verdict

**HYPOTHÈSE.** À implémenter avant toute phrase causale. Priorité 1.

---

## 3. Vérification et validation différentielles → **TESTABLE (partiellement fait)**

### Ce qui existe

| contrôle | valeur | artefact |
|---|---|---|
| tolérance \(10^{-8}\) vs \(10^{-9}\) sur \(I_D\) | 0,0018 % | `results/tolerance/tol_test_*.md` |
| stabilité \(I_D\), \(N=128\to256\) | 1,92 % (Re=400), 2,75 % (Re=1000) | `docs/06`, `figures/analysis/resolution_IQ_stability.json` |
| \(D^\*\) brut \(N=128\to256\) | +2,99 % | idem |
| \(P^\*\) brut \(N=128\to256\) | +6,69 % | idem |
| bilan \((P-D)/D\) à \(N=256\) | −3,07 → −3,47 % | `docs/07:34` |
| position du tourbillon, Re=400 | (0,556 ; 0,606) vs Ghia (0,5547 ; 0,6055) | `docs/05:32-33` |

Le test de stabilité du **contraste** répond exactement à l'objection du
relecteur (« la convergence d'une différence ne découle pas
automatiquement de la convergence de chacune des simulations ») : on
mesure directement \(\Delta I_Q\) au raffinement, pas seulement les
grandeurs absolues.

### Ce qui manque

- **2 grilles** : aucun ordre apparent, aucune extrapolation
  Richardson, aucun \(N=384\). `provenance/.../convergence_order.py`
  existe mais n'a jamais été lancé sur ce projet.
- **Validation FEM non faite** : `fem/fem_taylor_hood.py` est copié et
  fonctionnel, `docs/03:77` le marque « à faire », aucun résultat.
- **Profils de vitesse non comparés à Ghia** : seule la position du
  centre du tourbillon est testée ; les profils \(u(y),v(x)\) ne sont
  que dans `provenance/`.
- **Bilan hors cible** : plan initial \(|P-D|/D<0{,}02\) (`docs/03:33`),
  mesuré **3,1 %** à \(N=256\).
- **Pas de bilan de masse, pas de test automatisé** (recherche négative).

### Verdict

**TESTABLE — partielle et non achevée.** La *méthode* de vérification
différentielle est en place et fonctionne ; la *validation* externe
(FEM, profils de référence) et l'*ordre de convergence* restent à faire.
Ne pas écrire que le projet a « établi une norme méthodologique ».

---

## 4. Jeu de données de référence → **RESSOURCE REPRODUCTIBLE (pas encore un benchmark)**

### Ce qui existe

| ensemble | taille | état |
|---|---|---|
| campagne \(N=128\) | **152 cas** (19 Re × 8) | 152/152 convergés, 0 échec, tol \(10^{-8}\), Ma fixé à 0,1804 |
| campagne \(N=256\) | **24 cas** (3 Re × 8) | 24/24 ; Re=100 recalculé sur la loi de Ma le 18/10 |
| quarantaine | 70 + 9 + 9 fichiers | `results/quarantine/tau060_N{128,256}/`, cas hors-loi de Mach |
| scripts | `run_campaign`, `run_case`, `plot_resweep`, `plot_resolution` | entrées/CLI tracées |
| documentation | `docs/00` → `docs/09` | protocole, décisions, résultats, pièges |

Point fort : le **piège de Mach** (Note 08 §3.1) a été identifié, 70 cas
invalides mis en quarantaine plutôt que supprimés, et la correction est
rétro-compatible avec les cas de référence. C'est de la bonne hygiène de
données et ça se raconte.

### Ce qui manque pour parler de « benchmark »

- Référence **externe indépendante** (FEM, profils Ghia) — absente.
- Tests automatisés, bilan de masse, README d'utilisation, DOI/versionnage.

### Verdict

**RESSOURCE REPRODUCTIBLE — oui, dans l'état. BENCHMARK de référence —
non.** Formulation retenue : *dataset and analysis procedures intended
to support reproducibility*, pas *reference benchmark*.

---

## 5. Saturation de l'effet de forme pour \(n\ge4\) → **TESTABLE (rendements décroissants démontrés, saturation complète non démontrée)**

### Ce qui existe

Effet de l'ondulation \(\Delta D^\*=D^\*(\alpha{=}0{,}10)-D^\*(\alpha{=}0)\),
Re=400, \(N=128\) (`docs/05:19-26`) :

| \(n\) | \(D^\*(\alpha{=}0)\) | \(\Delta D^\*\) | incrément sur \(n\) précédent |
|---|---|---|---|
| 2 | 22,75808 | −0,11671 | — |
| 4 | 22,83148 | −0,15780 | −0,04109 (−35 %) |
| 6 | 22,83430 | −0,16170 | −0,00390 (−2,5 %) |
| 8 | 22,83429 | −0,16205 | −0,00035 (−0,2 %) |

Contrastes correspondants à \(N=256\) (`docs/07:51-56`) :

| paire | \(\lvert I_D\rvert\) | rapport à \(\lvert I(2,4)\rvert\) |
|---|---|---|
| (2,4) | 0,041894 | 1 |
| (4,6) | 0,004060 | **×10** |
| (4,8) | 0,004451 | ×9 |
| (6,8) | 0,000391 | **×107** |

### Comparaison avec l'incertitude (le point décisif)

Chaque incrément **est lui-même** un contraste apparié, mesuré aux deux
grilles. Stabilité \(N=128\to256\) mesurée sur les trois Reynolds :

| contraste | Re=100 | Re=400 | Re=1000 | lecture |
|---|---|---|---|---|
| \(I(2,4)\) | 4,20 % | 1,92 % | 2,75 % | **stable** |
| \(I(4,6)\) | 13,66 % | 3,96 % | 3,92 % | stable pour Re ≥ 400, **marginal à Re=100** |
| \(I(6,8)\) | 11,12 % | 10,54 % | 7,25 % | **à la limite de la résolution** |

Rapport signal/bruit pour Re=400 (l'incertitude est celle du contraste
lui-même, mesurée, pas une hypothèse) :

| incrément | \(\lvert I\rvert\) | incertitude | rapport | verdict |
|---|---|---|---|---|
| \(2\to4\) | 0,04189 | ±0,0008 | **×52** | résolu sans ambiguïté |
| \(4\to6\) | 0,00406 | ±0,00016 | **×25** | **détecté** |
| \(6\to8\) | 0,00039 | ±0,00004 | **×10** | marginal |

### Lecture physique : rendements décroissants, non saturation

Part de chaque incrément dans le changement total \(n=2\to8\)
(\(\Delta D^\*\) de −0,1167 à −0,1621 à Re=400) :

| Re | \(2\to4\) | \(4\to6\) | \(6\to8\) |
|---|---|---|---|
| 100 | 82,6 % | 14,7 % | 2,6 % |
| 400 | **90,6 %** | 8,6 % | 0,8 % |
| 1000 | **92,6 %** | 7,2 % | 0,2 % |

Le saut \(2\to4\) capte **83 à 93 %** de l'effet géométrique total.
Mais l'incrément \(4\to6\) reste **détectable** (×25 l'incertitude à
Re ≥ 400) et représente encore 7 à 15 % du changement : la saturation
**n'est pas complète en \(n=4\)**.

### Verdict

**TESTABLE — rendements décroissants démontrés, saturation complète non
démontrée.**

Formulations défendables :

- *the change from \(n=2\) to \(n=4\) accounts for 83–93 \% of the total
  geometric effect over \(n=2\ldots8\)* — **démontré** (×25 à ×52
  l'incertitude) ;
- *residual changes beyond \(n=4\) remain detectable at \(Re\ge400\) but
  represent less than 9 \% of the total* — **démontré** ;
- *the effect saturates for \(n\ge4\)* — **non soutenu** : à Re=100 la
  stabilité de \(I(4,6)\) n'est que de 13,7 %, et à Re=400 le résidu est
  significatif. Un relecteur refera ce calcul et conclura à une
  surinterprétation.

Il faut en outre signaler que \(I(6,8)\) n'est stable qu'à 7–11 % : ce
dernier incrément n'est pas résolu et ne doit pas être commenté.

---

## 6. Synthèse

| # | affirmation | verdict | condition pour écrire |
|---|---|---|---|
| 1 | non-additivité | **DÉMONTRÉ** | nommer \(Q\) explicitement (Re=100 déjà recalculé sur la loi de Ma) |
| 2 | médiation pariétale | **HYPOTHÈSE** | implémenter la décomposition (`docs/02` §7) |
| 3 | vérification différentielle | **TESTABLE** | FEM + profils Ghia + ordre de convergence + bilan < 2 % |
| 4 | jeu de données de référence | **RESSOURCE** (pas benchmark) | validation externe, tests, README |
| 5 | saturation \(n\ge4\) | **TESTABLE** | dire « rendements décroissants : 83–93 % du changement en 2→4 », pas « saturation pour \(n\ge4\) » |

### Ce qui peut entrer dans l'introduction sans risque

Le paragraphe proposé par la relecture, qui transforme chacune des
affirmations en question testable :

> *Although the effects of global cavity geometry and local wall
> corrugation have both received attention in confined-flow research,
> their combined energetic effects require a systematic assessment. In
> particular, the extent to which these geometric contributions can be
> treated as additive remains to be quantified for the configuration
> considered here. We address this question by introducing an
> interaction contrast that isolates the coupled contribution of global
> shape and local wall perturbation from their individual effects. The
> contrast is evaluated over a Reynolds-number sweep, with matched
> numerical configurations and a resolution-sensitivity assessment
> designed to distinguish small physical differences from
> discretization-induced biases. A complementary wall–core decomposition
> is used to investigate the spatial distribution of the interaction and
> assess whether its Reynolds-number dependence is consistent with a
> near-wall-mediated mechanism. Finally, the superelliptic geometry
> provides a continuous parameterization of global shape, enabling us to
> examine whether the interaction exhibits systematic geometric trends
> or saturation. The resulting numerical dataset and analysis procedures
> are intended to support reproducibility and future data-driven
> investigations of geometrically modified cavity flows.*

**Réserve :** cette version annonce une décomposition wall–core qui
**n'existe pas encore** (§2). Elle n'est rédigeable qu'une fois le code
écrit et exécuté, sinon il faut retirer la phrase « A complementary
wall–core decomposition is used… ».

### Plan d'action prioritaire (ordre de valeur pour le relecteur)

| priorité | vérification | effort | état |
|---|---|---|---|
| 1 | décomposition Wall/Core + robustesse au seuil \(\delta_w\) | code + 8×N réduction énergie | **non commencé** |
| 2 | validation FEM (`fem/fem_taylor_hood.py` déjà copié) | post-traitement | **non commencé** |
| 3 | profils \(u(y),v(x)\) vs Ghia à Re=400/1000 | post-traitement | **non commencé** |
| 4 | convergence différentielle sur **\(I_Q\)** à 3 résolutions (\(N=128,256,384\)) | 3 grilles | **2 grilles faites** |

Tant que 1 et 3 ne sont pas faites, les contributions 2 et 3 du texte
d'ambition doivent rester formulées comme objectifs et non comme résultats.
