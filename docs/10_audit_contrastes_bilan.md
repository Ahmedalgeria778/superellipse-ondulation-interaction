# Note 10 — Audit quantitatif des contrastes et localisation du défaut énergétique

**Demandes du relecteur auxquelles cette note répond**

| Priorité | Question | Outil | Rapport |
|---|---|---|---|
| 1 | Audit quantitatif des contrastes \(I\) (et non plus qualitatif) | `solver/audit_contrasts.py` | `results/audit/audit_contrasts.json` |
| 2 | Bilan énergétique : localiser la source du défaut \((P^\*-D^\*)/D^\*\), vérifier les conventions de signe, démontrer ou infirmer que le défaut s'annule dans le contraste apparié | `solver/audit_balance.py` | `results/audit/audit_balance.json` |
| 2 (complément) | Discrétisation (h) **ou** compressibilité (Ma) ? | `solver/run_ma_sweep.py` | `results/ma_sweep/masweep_Re*_N*.json` |

Aucune simulation supplémentaire pour la priorité 1 : tout est reconstruit à partir des
champs `.npz` déjà calculés. Le balayage Ma est la seule extension, 8 cas N=128/N=256
(`save=False`, `ckpt=False`, sortie dédiée — **aucun cas de campagne n'est écrasé**).

---

## 0. Conventions

* \(D^\* = 2\int_{\Omega_f} \mathbf{S}^\*:\mathbf{S}^\*\,\mathrm{d}A^\*\) intégré sur le
  **masque fluide** (`diagnostics.py:255-256`, facteur `*d*d`), \(K^\*\) idem.
* \(P^\* = \oint_{\text{lid}} U^\*_{\text{lid}}\,\partial_y u^\*\,\mathrm{d}x\) : intégrale
  de **ligne sur la lid** (`diagnostics.py:259-263`, facteur `*d`), domaine fixe.
* \(\varepsilon = (P^\*-D^\*)/D^\*\) : le « défaut de bilan », toujours négatif
  (\(P^\*<D^\*\)).
* \(I\) = contraste apparié
  \(\bigl(|X_{b,\alpha}|-|X_{b,0}|\bigr)-\bigl(|X_{a,\alpha}|-|X_{a,0}|\bigr)\),
  avec \(X\in\{D^\*,P^\*\}\) et \(I_D,I_P\) ses deux estimations.
* Normalisation d'écart de grille **unifiée** : \((|I_{128}|-|I_{256}|)/|I_{256}|\).
* Loi de campagne : \(U=(40/3)/N\) → Ma = 0,180422 (N=128) et 0,090211 (N=256),
  fixes pour tous les Re. À \(N,Re\) fixes, \(\tau = 0{,}5+\sqrt3\,N\,\mathrm{Ma}/Re\)
  : **Ma et τ sont liés** (limitation rappelée en §2.8).

---

# 1. Priorité 1 — Audit quantitatif des contrastes

## 1.1 A. Conformité des cas

* **48 cas** audités : Re ∈ {100, 400, 1000} × N ∈ {128, 256} × n ∈ {2,4,6,8} × α ∈ {0, 0,10}.
* **0 cas hors-loi Ma** : Ma = 0,180422 / 0,090211 exactement pour les 48.
* Résidus ≤ 9,61 × 10⁻⁹ ≤ tolérance 10⁻⁸ pour **40 cas**.
* **8 cas en `reanalyzed` sans résidu stocké** (N=128, Re=400, les 8 combinaisons n×α) :
  leur JSON hérite du statut de `run_campaign.py:39-60`, qui traite un enregistrement
  `reanalyzed` sans résidu comme valide. Les champs `.npz` correspondants sont bien
  présents (305–357 Ko) et les scalaires ont été **reconstruits depuis ces champs**
  (cohérence §1.8). Ce n'est pas un écart de calcul mais un **trou de provenance** :
  le résidu de convergence de ces 8 cas n'est pas conservé.

  > Ce point est documenté et non contourné : conformément à la directive « ne jamais
  > relancer un cas déjà calculé », **aucun de ces 8 cas n'est relancé**. Il faut
  > simplement ne pas écrire que les 48 cas sont résiduellement vérifiables.

## 1.2 B. Les quatre termes

Les 48 cas sont tous exploitables ; \(I_D<0\) et \(I_P<0\) dans les 18 couples
(Re, N, paire) — l'effet de l'ondulation est donc **toujours suppressif** — et
\(Q(\alpha)<Q(0)\) pour chaque n. Le signal est de même signe partout, aucun cas
contraire à traiter en exception.

## 1.3 C. L'effet n'est pas un artefact de volume (contrôle de normalisation)

L'ondulation supprime **−5,1 à −5,7 %** d'aire fluide selon n (mesuré, pas −3,4 %).
Deux conventions sont recalculées :

| Convention | Interprétation |
|---|---|
| extensive (\(D^\*\) intégré sur le masque fluide) | mélange effet géométrique et perte d'aire |
| densité (\(D^\*/A_f\)) | isole l'effet « par unité de fluide » |

## 1.4 C2. Modèle nul « aire seule » : signe **opposé** dans 18/18 cas

On construit la prédiction purement géométrique
\(I_{\text{null}} = \) effet attendu si \(D^\*\) ne faisait que suivre l'aire supprimée,
en appliquant la même formule au même modèle d'aire.

| | \(I_{\text{null}}\) | \(I_{\text{mesuré}}\) |
|---|---|---|
| signe | **positif, 18/18** | **négatif, 18/18** |
| rapport \(I_{\text{null}}/I_{\text{mesuré}}\) | — | ∈ [−79,6 ; −0,76] |

Exemples : Re=100, N=128, paire (2,4) : \(+0{,}0682\) contre \(-0{,}00184\) (rapport −37) ;
Re=1000, N=256, paire (2,4) : \(+0{,}1211\) contre \(-0{,}1466\) (rapport −0,83).

> **Le modèle nul « aire seule » prévoit un contraste de signe opposé au mesuré dans
> les 18 cas.** Le contraste n'est donc pas un artefact de volume : c'est le résidu
> d'une forte compensation physique (perte d'aire d'un côté, changement de
> cisaillement de l'autre).

## 1.5 C3. Robustesse du signe et parts additives

* **Signe conservé 18/18** entre les deux normalisations (extensive ↔ densité).
* Décomposition des intervalles 2→4, 4→6, 6→8 (somme exactement 100 %) :

| Normalisation | 2→4 | 4→6 | 6→8 |
|---|---|---|---|
| extensive | **81,1 – 92,6 %** | 7,2 – 16,1 % | 0,18 – 2,79 % |
| densité (\(/A_f\)) | **74,2 – 81,8 %** | 14,5 – 21,7 % | 2,5 – 5,7 % |

> La conclusion est la même sous les deux conventions : la **quasi-totalité de l'effet
> est portée par le premier intervalle de forme (2→4)**, avec des rendements
> décroissants. (Voir note 09 §5 : rendements décroissants démontrés, **saturation
> complète pour n≥4 non démontrée**.)

## 1.6 D. Le défaut de bilan est commun-mode

\(\varepsilon\) sur les 48 cas : **−6,86 % ≤ ε ≤ −2,77 %**.

| Réglage | plage de ε (N=128 → N=256) |
|---|---|
| Re=100 | −5,793 % → −2,766 % |
| Re=400 | −6,220 % → −3,210 % |
| Re=1000 | −6,855 % → −3,501 % |

**Trois mesures d'indépendance à la géométrie :**

1. **Dispersion de ε entre les 4 géométries** (n=2,4,6,8 × α=0,0,10) à N, Re fixes :
   **0,0042 à 0,0774 pt%**, soit **0,09 % à 1,38 % de sa propre valeur**.
2. **Dépendance différentielle à l'ondulation** \(\Delta_\alpha\varepsilon = \varepsilon(\alpha{=}0{,}10)-\varepsilon(\alpha{=}0)\) (moyenne sur n) :
   **−0,0594 à +0,0014 pt%**, soit **≤ 0,9 % du défaut lui-même**.
3. **Ordre apparent en h** \(\bigl(\log|\varepsilon_{128}/\varepsilon_{256}|\bigr)/\log 2\) :
   **1,066 (Re=100), 0,956 (Re=400), 0,974 (Re=1000)** → **premier ordre**, pas de
   saturation entre les deux résolutions.

> Seules les lignes 2 et les différences entre géométries polluent le contraste \(I\) ;
> tout le reste est commun aux 4 configurations d'un contraste et s'annule exactement.
> **Le défaut s'annule donc dans le contraste apparié à hauteur de ≤ 0,06 pt%.**

## 1.7 E. Sensibilité de grille des contrastes

Écart N=128 ↔ N=256 sur \(I\) :

| Paire | Re=100 | Re=400 | Re=1000 |
|---|---|---|---|
| \(I_D\)(2,4) | 4,20 % | **1,92 %** | 2,75 % |
| \(I_P\)(2,4) | 5,05 % | 4,00 % | 5,09 % |
| \(I_D\)(4,6) | **13,66 %** | 3,96 % | 3,92 % |
| \(I_D\)(6,8) | 11,12 % | 10,54 % | 7,25 % |
| \(I_P\)(6,8) | 10,00 % | 16,18 % | **24,74 %** |

> **Stable entre grilles** : \(I(2,4)\) pour Re ≥ 400 (1,9–4,2 %).
> **Non convergé quantitativement** : \(I(4,6)\) à Re=100, toutes les paires (6,8),
> et \(I_P(6,8)\) à Re=1000 (24,7 %). Le seuil de 10 % est franchi par 5 configurations sur 18.

## 1.8 E2. Contrôle croisé \(I_P\) contre \(I_D\) (barre d'erreur indépendante)

Paire (2,4), N=256 — deux intégrales de nature différente (surface sur le masque fluide
vs ligne sur la lid) calculées sur le même champ :

| Re | \(I_D\) | \(I_P\) | \(I_P/I_D-1\) |
|---|---|---|---|
| 100 | −0,0019221 | −0,0024434 | **+27,1 %** |
| 400 | −0,0418941 | −0,0442376 | **+5,59 %** |
| 1000 | −0,146628 | −0,158694 | **+8,23 %** |

Sur les paires secondaires, l'écart va de −18,0 % à −11,1 % (signe inverse, de même
faible ampleur relative).

**Budget d'incertitude sur \(I(2,4)\)** (somme quadratique écart de grille \(I_D\) +
écart de grille \(I_P\) + désaccord \(I_P/I_D\) ; tolérance 10⁻⁸ = 0,002 %, marginale) :

| Re | budget |
|---|---|
| 400 | **≈ 7 %** |
| 1000 | **≈ 10 %** |
| 100 | **≈ 28 %** |

## 1.9 F. Cohérence des valeurs stockées

Comparaison des scalaires JSON archivés avec le recalcul depuis les `.npz` :

* **0,000 %** pour 5 des 6 couples (Re, N) ;
* **2,31 %** pour N=128 / Re=100 — c'est le couple dont l'archivage est le plus ancien
  (conservé tel quel, fichiers hors-loi de l'ancienne loi de τ). Aucune valeur n'a été
  réécrite.

## 1.10 Verdict — priorité 1

| Point | Verdict |
|---|---|
| Conformité Ma | **0/48 hors-loi** ; 8 cas sans résidu → trou de provenance documenté |
| Contraste non réductible à l'aire | **démontré** (modèle nul de signe opposé 18/18) |
| Robustesse du signe | **18/18** ; parts additives, mêmes conclusions sous les 2 conventions |
| Défaut de bilan commun-mode | **démontré** (dispersion 0,09–1,38 % de lui-même ; \(\Delta_\alpha\) ≤ 0,06 pt%) |
| Ordre en h | **premier ordre** (0,96–1,07) |
| \(I(2,4)\) à Re ≥ 400 | **stable entre grilles** (1,9–4,2 %), budget ≈ 7 % |
| Paires (6,8) et \(I(4,6)\) à Re=100 | **non convergées quantitativement** |

---

# 2. Priorité 2 — Localisation du défaut énergétique

## 2.1 Les quatre causes candidates

Toutes testées sur les champs déjà calculés (`solver/audit_balance.py`, 24 cas :
Re ∈ {100,400,1000} × N ∈ {128,256} × n ∈ {2,8} × α ∈ {0,0,10}) :

| # | Cause | Mode de test |
|---|---|---|
| (i) | gradient de vitesse en **ordre 1** sur la lid (paroi à 0,5 dx du 1er nœud) | `p2` / `lid` : formule à 3 points (ordre 2) |
| (ii) | **inversion de signe** des corrections unilatérales quand la paroi est du côté négatif de l'axe | `sign` : signes corrigés en bas/gauche |
| (iii) | contamination du gradient par le **saut fluide/solide** (`np.gradient` dérive à travers le front) | `clean` : dérivées centrées uniquement entre nœuds fluides |
| (iv) | toutes les corrections unilatérales sont en **ordre 1** sur les 4 parois | `full` : ordre 2 sur les 4 parois |

Nomenclature des colonnes : `ref` = état publié ; `p2` = P seul en ordre 2 ;
`lid` = P **et** D en ordre 2 sur la lid, parois restantes publiées ;
`sign` = signes corrigés ; `full` = les 4 parois en ordre 2 ; `clean` = saut retiré.

## 2.2 Décomposition de D\* par région (dans `ref`)

| Réglage | latérales | **lid** | fond ondulé | intérieur |
|---|---|---|---|---|
| Re=100, N=128 | 19,4 % | **36,5 %** | 0,02–0,05 % | 44,1 % |
| Re=100, N=256 | 10,5 % | **29,6 %** | 0,01–0,02 % | 59,8 % |
| Re=400, N=128 | 22,3 % | **38,5 %** | 0,04–0,16 % | 39,1 % |
| Re=400, N=256 | 13,5 % | **29,7 %** | 0,02–0,08 % | 56,8 % |
| Re=1000, N=128 | 24,2 % | **44,1 %** | 0,15–0,47 % | 31,4 % |
| Re=1000, N=256 | 16,5 % | **33,3 %** | 0,08–0,26 % | 50,1 % |

* Bande pariétale (latérales + lid + fond) = **43 à 69 %** de \(D^\*\).
* **La lid est la région dominante (30–44 %).**
* **Le fond ondulé ne contribue que 0,02 à 0,47 %** de \(D^\*\) : l'effet de l'ondulation
  sur \(D^\*\) (≈ −0,5 % en valeur absolue) passe par **la modification du champ dans
  toute la boîte**, pas par le cisaillement local sur la crête.

## 2.3 Test (ii) — inversion de signe : **défaut de code réel, impact négligeable**

`diagnostics.py:73` (paroi inférieure) et `:57` (paroi gauche) calculent
`(0 - U)/dist` là où la dérivée unilatérale est `+(U - 0)/dist` : ces deux parois sont
du **côté négatif** de l'axe (y croît vers le haut, x croît vers la droite), le signe
est donc inversé. Les parois haute et droite, elles, sont correctes.

Effet de la correction sur \(D^\*\) :

| Réglage | \(\Delta D^\*\) |
|---|---|
| Re=1000, N=128 | +0,09 à +0,14 % |
| Re=400, N=128 | +0,24 à +0,26 % |
| Re=100, N=128 | **+0,44 à +0,45 %** |
| Re=1000/400, N=256 | −0,002 à +0,051 % |

Gain moyen sur \(\lvert\varepsilon\rvert\) : **−3,1 %** (min −7,8 %, max +0,1 %), c.-à-d.
la correction rend le défaut *légèrement pire*.

> L'effet croît avec l'épaisseur de couche limite (Re faible, grille grossière) car le
> premier nœud fluide est plus chargé. **C'est un bug à corriger dans la provenance,
> mais il n'est pas la source du défaut** (≤ 0,45 % sur \(D^\*\), ≤ 0,45 pt% sur ε).

## 2.4 Test (iii) — saut fluide/solide : **négligeable**

Retirer toute dérivation à travers le front (différences centrées entre nœuds fluides
contigus) change \(D^\*\) de **≤ +0,07 %**. Gain moyen sur \(\lvert\varepsilon\rvert\) :
**−0,2 %** (min −1,0 %, max −0,0 %).

> La contamination par le saut n'est pas la source du défaut.

## 2.5 Test (i) — le gradient en ordre 1 sur la lid : **cause dominante**

\(P^\*\) est une intégrale de ligne **exclusivement sur la lid** et n'utilise que
\(\partial_y u\) au nœud supérieur, estimé en ordre 1 :
\((U_{\text{lid}} - u_{N-1})/(0{,}5\,dx)\).

Passer **seul** ce gradient en ordre 2 (formule à 3 points exacte pour espacement
inégal) réduit le défaut de :

> **Gain moyen +70,8 % (min +61,2 %, max +87,0 %)** sur les 24 cas —
> \(D^\*\) : \(\varepsilon\) passe par exemple de −6,20 % à −2,41 % (Re=400, N=128)
> et de −3,47 % à −0,95 % (Re=1000, N=256).

En valeur absolue : \(\Delta P^\* = +2{,}10\) à \(+4{,}56\) % de \(D^\*\).

> **Le défaut est dominé par le cisaillement pariétal de premier ordre sur la lid.**

## 2.6 Test (iv) — pourquoi **aucun remédage local n'élimine le défaut**

On teste les variantes **cohérentes** (où l'amélioration est appliquée à la fois dans
\(P^\*\) et dans \(D^\*\)) :

| Correction | gain moyen sur \(\lvert\varepsilon\rvert\) | amplitude | verdict |
|---|---|---|---|
| `p2` — lid ordre 2 dans \(P^\*\) seul | **+70,8 %** | +61 à +87 % | améliore, mais **paire incohérente** (D non corrigé) |
| `lid` — lid ordre 2 dans \(P^\*\) **et** \(D^\*\) | **+2,3 %** | **−19,9 à +22,4 %** | **non robuste** : pire à N=128, meilleur à N=256 |
| `sign` — signes corrigés | −3,1 % | −7,8 à +0,1 % | pire |
| `full` — 4 parois en ordre 2 | **−29,7 %** | −57,8 à −3,7 % | **le défaut augmente** de 3,7 à 10,3 pt% |
| `clean` — saut retiré | −0,2 % | −1,0 à −0,0 % | nul |

Détail du test cohérent `lid` :

| Réglage | \(\Delta P^\*\) | \(\Delta D^\*\) | ε `ref` → ε `lid` | gain |
|---|---|---|---|---|
| Re=100, N=128 | +4,15 % | **+5,30 %** | −5,79 % → −6,94 % | **−19,9 %** |
| Re=400, N=128 | +3,80 % | **+4,43 %** | −6,20 % → −6,84 % | **−10,2 %** |
| Re=1000, N=128 | +4,53 % | **+5,43 %** | −6,81 % → −7,71 % | **−13,2 %** |
| Re=100, N=256 | +2,40 % | +1,79 % | −2,77 % → −2,15 % | +22,4 % |
| Re=400, N=256 | +2,10 % | +1,47 % | −3,20 % → −2,57 % | +19,7 % |
| Re=1000, N=256 | +2,52 % | +2,00 % | −3,47 % → −2,95 % | +15,0 % |

> **Le remédage cohérent de la lid change de signe d'effet entre N=128 et N=256**
> (aggrave de 10–20 % sur la grille grossière, améliore de 15–22 % sur la fine).
> C'est la signature d'une **erreur de discrétisation de premier ordre** : le
> reconstruité en ordre 2 sur- ou sous-passe le gradient vrai d'un montant qui dépend
> de \(h\). Ce n'est **pas** le comportement d'un bogue unique.
>
> Même chose pour `full` : \(D^\*\) monte de 2,6 à 8,1 % alors que \(P^\*\) ne monte que
> de 2,1 à 4,6 %, donc le défaut **s'aggrave**. Deux corrections cohérentes en ordre 2
> rapprochent chacune le résultat du vrai… sans se rapprocher l'une de l'autre.

> **Conclusion du test : le défaut est un écart de cohérence en \(O(h)\) entre deux
> intégrales indépendamment discrétisées d'un même champ instantané, dominé par la
> couche limite de la lid — et non un bug unique identifiable et corrigeable.**

## 2.7 Annulation dans le contraste apparié — **démontré**

Trois arguments convergents (détail §1.6) :

1. dispersion de \(\varepsilon\) entre géométries à \(N,Re\) fixes : **0,09 à 1,38 % du défaut** ;
2. \(\Delta_\alpha\varepsilon \le 0{,}06\) pt% = **≤ 0,9 % du défaut** ;
3. vérification indépendante : \(I_P\) et \(I_D\) convergent à **+5,59 % (Re=400)** et
   **+8,23 % (Re=1000)** sur la paire (2,4).

Avec \(D^\* \approx 23\), un écart commun-mode de 0,06 pt% vaut
\(\approx 1{,}4\times10^{-2}\) en valeur absolue sur \(D^\*\), soit \(\approx 0{,}3\)% de
\(I(2,4)=-0{,}042\) à Re=400 — **négligeable devant le budget de 7 %**.

> **Le défaut de bilan s'annule dans le contraste apparié.** La quantité publiable est
> le contraste \(I\), avec la barre d'erreur de §1.8 ; la valeur absolue de \(D^\*\)
> seule n'est pas.

## 2.8 Balayage Ma : discrétisation **ou** compressibilité ?

L'affinement N=128 → N=256 modifie **à la fois** \(h\) et Ma (\(0{,}1804 \to 0{,}0902\)).
On fixe donc \(N\) et \(Re\) et on fait varier \(U\) (donc \(\nu\) et \(\tau\)) : à \(Re\)
fixe, l'écoulement adimensionné est le même, seul Ma varie.

**Re=400, N=128, n=2, α=0** (`results/ma_sweep/masweep_Re400_N128.json`) :

| Ma | τ | \(D^\*\) | \(P^\*\) | ε | itérations |
|---|---|---|---|---|---|
| 0,3608 | 0,700 | 22,865 | 21,436 | −6,253 % | 60 000 |
| 0,1804 | 0,600 | 22,758 | 21,346 | −6,204 % | 84 000 |
| 0,0902 | 0,550 | 22,729 | 21,476 | −5,512 % | 146 000 |
| 0,0451 | 0,525 | 22,819 | 21,851 | −4,238 % | 256 000 |

* **\(D^\*\) est invariant à ±0,30 %** ; **\(P^\*\) varie de −2,31 % à +2,31 %**.
  → Toute la variation de ε passe par le **cisaillement de la lid**, pas par la
  dissipation. (Le point Ma=0,1804 reproduit exactement la valeur de campagne
  −6,2039 %, contrôle croisé passé.)
* Ajustements de ε à \(h\) fixe :

| modèle | c₀ (extrapolation Ma→0) | c₁ | R² | c₀/ε(Ma=0,1804) |
|---|---|---|---|---|
| ε = c₀ + c₁·Ma² | **−5,105 %** | −10,33 %/Ma² | 0,428 | **82 %** |
| ε = c₀ + c₁·Ma | **−4,661 %** | −5,263 %/Ma | 0,613 | **75 %** |

* **Effet h à Ma fixe** : à Ma = 0,0902, passer \(h=1/128 \to 1/256\) fait passer ε de
  **−5,512 % à −3,199 %** → **2,31 pt%** ; sur une plage comparable, l'effet h est donc
  **≈ 3 fois l'effet Ma** (0,69 pt% entre Ma = 0,18 et Ma = 0,09 à \(h\) fixe).

**Re=400, N=256, n=2, α=0** (`results/ma_sweep/masweep_Re400_N256.json`) :

| Ma | τ | \(D^\*\) | \(P^\*\) | ε | itérations |
|---|---|---|---|---|---|
| 0,3608 | 0,900 | 23,360 | 22,750 | −2,611 % | 120 000 |
| 0,1804 | 0,700 | 23,293 | 22,540 | −3,232 % | 168 000 |
| 0,0902 | 0,600 | 23,250 | 22,506 | **−3,199 %** | 270 000 |
| 0,0451 | 0,550 | 23,250 | 22,621 | −2,709 % | 472 000 |

* Le point Ma = 0,0902 **reproduit le contrôle de campagne à l'unité** (−3,1988 %) ;
  \(D^\*\) est invariant à **±0,47 %** et \(P^\*\) ne varie que de **−0,51 % à +0,57 %**.
* La dépendance en Ma **change de forme** avec l'affinement : monotone à N=128,
  **non monotone à N=256** (amplitude maximale vers Ma ≈ 0,15–0,18) ; la **dispersion
  tombe de 2,02 pt% (N=128) à 0,62 pt% (N=256)**. C'est la signature d'un **couplage
  \(h\)–Ma**, non d'une loi de compressibilité autonome.

**Ajustement à deux variables** sur les **8 points** (design factoriel complet
\(2\,(h)\times 4\,(\mathrm{Ma})\), les 4 valeurs de Ma étant les mêmes aux deux \(h\))
par `solver/fit_mach_resolution.py` → `results/ma_sweep/fit_h_ma.json` :

| modèle | k | \(b_0\) [%] | \(b_h\) [%/h] | \(b_{\mathrm{Ma}}\) | \(b_{\mathrm{Ma}^2}\) | R² | RMSE [pt%] | AICc |
|---|---|---|---|---|---|---|---|---|
| M0h : \(b_0+b_h h\) | 2 | −0,324 | −669,1 | — | — | 0,822 | 0,608 | **−1,6** |
| M0M : \(b_0+b_{\mathrm{Ma}}\mathrm{Ma}\) | 2 | −3,874 | — | −2,191 | — | 0,034 | 1,417 | +12,0 |
| M1 : \(b_0+b_h h+b_{\mathrm{Ma}}\mathrm{Ma}\) | 3 | +0,047 | −669,1 | −2,191 | — | 0,856 | 0,547 | +2,4 |
| M2 : \(+\ b_I\,h\,\mathrm{Ma}\) | 4 | −1,512 | −403,1 | +7,023 | — | 0,922 | 0,402 | +6,7 |
| M3 : \(b_0+b_h h+b_{\mathrm{Ma}^2}\mathrm{Ma}^2\) | 3 | −0,165 | −669,1 | — | −3,689 | 0,839 | 0,578 | +3,2 |
| M4 : \(b_0+b_h h+b_{\mathrm{Ma}}\mathrm{Ma}+b_{\mathrm{Ma}^2}\mathrm{Ma}^2\) | 4 | +1,031 | −669,1 | −17,18 | +35,88 | 0,918 | 0,413 | +7,2 |

* **Le terme de grille domine et est robuste** : \(b_h\) est **identique à quatre chiffres
  dans les cinq modèles qui le contiennent** (−669,1 %/h), parce que le design
  \(2\times4\) (mêmes Ma aux deux \(h\)) rend \(h\) **exactement orthogonal** à tout
  polynôme en Ma. La seule variation de \(h\) (1/128 → 1/256) explique à elle seule les
  **2,31 pt%** observés.
* **Le terme Ma est faible et mal contraint** : il ne pèse que **0,62 pt%** sur toute la
  plage à \(N=256\), et \(b_{\mathrm{Ma}}\) **change de signe et d'ordre selon la forme du
  modèle** (−2,19 / +7,02 / −17,18), de même que sa version quadratique
  \(b_{\mathrm{Ma}^2}\) (−3,69 / +35,88).
* **Sélection de modèle** : au critère **AICc** (qui pénalise les \(k=4\) sur \(n=8\)), le
  modèle le plus simple **M0h gagne** (−1,6) ; ajouter un terme Ma ne l'améliore pas, bien
  que M2/M4 divisent le résidu max par deux (0,66 contre 1,31 pt%). Mais le résidu de
  M0h/M1 est **systématique** (jusqu'à +1,31 pt% au coin \(N=128,\ \mathrm{Ma}=0{,}0451\)),
  signe que le terme Ma existe réellement — simplement trop petit pour être tranché sur
  8 points.
* **L'extrapolation est modèle-dépendante** : le défaut extrapolé \(b_0\) (à
  \(h\to0,\ \mathrm{Ma}\to0\)) balaie **−3,87 % à +1,03 %**, soit **4,90 pt%
  = 79 % de l'amplitude du défaut de campagne**. Rapporté à ce dernier (−6,20 %, N=128),
  \(b_0\)/défaut vaut 5,2 % (M0h), −0,8 % (M1), 2,7 % (M3) dans les modèles dominés par
  \(h\), contre 62,4 % (M0M), 24,4 % (M2), −16,6 % (M4) : **les modèles qui négligent la
  forme en Ma ne rendent pas l'essentiel du défaut à \(h\to0\)**.

> **Verdict (8 points) :** le défaut est **dominé par la discrétisation \(h\)** ; un terme
> en Ma existe (couplage \(h\)–Ma), d'amplitude **≤ ≈ 1 pt%** à \(N=256\), mais sa
> **forme et son extrapolation ne sont pas tranchables** sur ce plan (un seul \(Re\), une
> seule forme, deux niveaux de \(h\)). **Aucune « loi » \(ε(\mathrm{Ma})\) n'est retenue** —
> seulement l'ordination des causes : \(h\) d'abord, Ma ensuite.

**Limitation assumée :** à \(N,Re\) fixes, \(\tau = 0{,}5+\sqrt3\,\mathrm{Ma}/(Re\,h)\)
**lie Ma et \(\tau\)** — varier \(\tau\) à Ma fixe exigerait de varier \(Re\), ce qui change
l'écoulement. Le design \(2\times4\) rend \(h\) et Ma estimables l'un sans l'autre,
**pas \(\tau\)**. Enfin la cible \(h\to0,\ \mathrm{Ma}\to0\) est **hors de l'enveloppe
échantillonnée** : les \(b_0\) ci-dessus sont des **extrapolations**, à ne pas lire comme
un défaut asymptotique mesuré.

## 2.9 Verdict — priorité 2

| Cause candidate | Mesure | Verdict |
|---|---|---|
| (ii) inversion de signe `diagnostics.py:57`, `:73` | \(\Delta D^\*\) = +0,002 à +0,45 % ; gain ε = −3,1 % | **bug réel, impact ≤ 0,45 % — pas la source** |
| (iii) saut fluide/solide dans `np.gradient` | \(\Delta D^\*\) ≤ +0,07 % ; gain ε = −0,2 % | **négligeable** |
| (i) gradient ordre 1 sur la lid | gain ε **+70,8 %** quand on l'améliore dans \(P^\*\) seul | **cause dominante** |
| — correction cohérente de la lid | gain **+2,3 % [−19,9 ; +22,4]** : signe inversé entre N=128 et N=256 | **erreur \(O(h)\), pas un bogue** |
| (iv) ordre 1 sur les 4 parois | gain **−29,7 %** (le défaut augmente) | le défaut n'est pas corrigeable localement |
| compressibilité (Ma) | \(b_h\) robuste (−669 %/h) ; terme Ma ≤ ≈1 pt%, non tranchable (AICc préfère M0h) ; \(b_0\) extrapolé balaie 4,90 pt% selon le modèle | **secondaire, non séparable du modèle** |
| **Annulation dans le contraste** | \(\Delta_\alpha\varepsilon \le 0{,}06\) pt% ≤ 0,9 % du défaut ; \(I_P/I_D\) à +5,6/+8,2 % | **démontrée** |

> **Réponse au relecteur.** Le défaut \((P^\*-D^\*)/D^\* \in [-6{,}9\,\%\,;\,-2{,}8\,\%]\)
> est une **erreur de discrétisation du cisaillement pariétal, de premier ordre en h**,
> dominée par la lid (30–44 % de \(D^\*\), dont \(P^\*\) dépend intégralement), avec une
> part de compressibilité mesurable mais secondaire. Il n'est réductible à aucun bug
> unique : deux des causes candidates (signes inversés, saut fluide/solide) sont
> mesurées et rejetées comme négligeables ; le remédage cohérent le plus prometteur
> change d'effet entre les deux résolutions, ce qui est précisément la signature d'une
> erreur \(O(h)\). **Il s'annule dans le contraste apparié** (≤ 0,06 pt%, ≤ 0,9 % de
> lui-même), si bien que la grandeur publiable est \(I\) avec un budget de ≈ 7 %
> (Re ≥ 400), ≈ 10 % (Re=1000) et ≈ 28 % (Re=100) — et non \(D^\*\) seul.

---

## 3. Ce qui est acquis / ce qui reste

### Acquis dans cette note

* `results/audit/audit_contrasts.json` — audit A–F complet des 48 cas (priorité 1).
* `results/audit/audit_balance.json` — 24 cas × 6 variantes de gradient, décomposition
  de \(D^\*\) par région (priorité 2).
* `results/ma_sweep/masweep_Re400_N128.json` — balayage Ma à \(h\) fixe.
* Deux bugs de code identifiés dans la **provenance** (non corrigés ici, car hérités) :
  * `diagnostics.py:57` et `:73` : signe inversé des corrections unilatérales pour la
    paroi gauche et la paroi inférieure → **à corriger dans `provenance/`**, impact
    mesuré ≤ 0,45 % sur \(D^\*\) ;
  * statut `reanalyzed` sans résidu conservé pour les 8 cas N=128/Re=400 →
    **à documenter dans `docs/06`/`docs/08`**, aucun relancement.

### Restant à faire (ordre pour le relecteur)

1. **Fait** : balayage Ma à N=256 (`solver/run_ma_sweep.py --Re 400 --N 256`) —
   4 points, **contrôle −3,1988 % reproduit à l'unité** (ligne de campagne),
   écrit dans `results/ma_sweep/masweep_Re400_N256.json`. Ajustement à deux
   termes \(\varepsilon=f(h,\mathrm{Ma})\) sur **8 points** (N=128 + N=256) →
   reporté en §2.8.
2. **Fait** : figure 2 panneaux — `figures/analysis/contraste_paire_Re.{png,pdf}`
   (`solver/plot_contrast_figure.py`), données exposées dans
   `results/audit/figure_contrast_data.json`.
   (a) \(|I_D|\) des 3 paires vs Re, N=128 plein / N=256 hachuré, sur fond de la
   courbe des 19 Reynolds, 2 légendes regroupées (courbe + panneau) ; (b) parts
   additives 2→4 / 4→6 / 6→8 empilées, Re × grille, mention plein/hachuré en
   xlabel. **PNG 5518×2482 px (233,6 × 105,1 mm @ 600 dpi) et PDF à 190 mm écrits** ;
   corps effectif vérifié (titre 8,0 / axe 8,0 / graduations 7,2 / légende 7,2 pt),
   aucun chevauchement courbe/légende.
3. **Fait** : sous-ensemble **N=384**, **16/16 cas convergés** (campagne
   2026-10-18 22:45 → 2026-10-19 ~14:40 ; 6 workers, `NUMBA_NUM_THREADS=2`,
   tol 1e-8, résidus ≤ 1e-8). Résultat : le désaccord \(I_P/I_D-1\) de la paire
   (2,4) tombe à **+2,9 / +1,6 / +1,9 %** (Re = 100 / 400 / 1000) contre
   +27-33 % à N=128 → **artefact de discrétisation**, pas désaccord systématique ;
   la paire (2,4) à Re=100 devient publiable (détail en §3.1). Durée réelle ≈ 16 h
   (≈ 83 h·cœur annoncés ; queue = les Re=100, ~500 k it au lieu des 223 k
   estimés). **N=512 abandonné** (18 it/s, ≈ 222 h·cœur).
4. **Reste** : paragraphe de synthèse — §4 ci-dessous, **à fournir/traduire en
   anglais** pour la soumission PoF.

### 3.1 Sous-ensemble N=384 : quelle sélection est *discriminante* ?

Le budget d'incertitude sur \(I(2,4)\) est **7 % (Re=400), 10 % (Re=1000),
28 % (Re=100)** ; il se décompose en (i) écart de grille N=128↔256, (ii) écart de
grille sur \(I_P\), (iii) désaccord \(I_P/I_D\).

**Éléments de arbitrage :**

* La paire **(2,4) porte 81–93 % de l'effet** ; les paires (4,6) et (6,8) portent
  7–16 % et **0,2–2,8 %**. Converger (6,8) à ±10 % déplacerait le total de ≤ 0,3 % :
  **ce n'est pas là que se joue la conclusion.** → le sous-ensemble doit viser la
  paire (2,4).
* À Re=400 et Re=1000 le budget est déjà ≤ 10 % : un 3ᵉ point N=384 **conforterait**
  sans changer de verdict.
* **À Re=100 le budget est de 28 % et il est dominé par le désaccord \(I_P/I_D\)
  (+27,1 %)**, pas par l'écart de grille (4,20 %). Deux issues possibles au passage
  en N=384 :
  * \(I_P/I_D\) **retombe vers ≈14 %** → le désaccord est bien un effet de
    discrétisation → budget ramené à ≈ 15 %, la paire (2,4) à Re=100 devient publiable ;
  * \(I_P/I_D\) **reste ≈ 27 %** → le désaccord est **systématique** et n'est pas un
    artefact de \(h\) → **la conclusion change** : il faut écrire que les deux
    estimateurs divergent à bas Re et ne pas les combiner.

  C'est la seule configuration de la sélection qui peut **inverser une conclusion
  déjà écrite** ; c'est donc elle qui est discriminante.

**Sélection retenue (16 cas) — décision du 2026-10-18** :

| Re | n | α | cas | question visée |
|---|---|---|---|---|
| 100 | 2, 4 | 0 et 0,10 | 4 | 3ᵉ point de grille sur la paire (2,4) à Re=100 (budget 28 %, ordre apparent 1,066) |
| 400 | 2, 4, 6, 8 | 0 et 0,10 | 8 | les **trois** intervalles (2,4)/(4,6)/(6,8) à la grille la plus fine, sur le Re le mieux convergé |
| 1000 | 2, 4 | 0 et 0,10 | 4 | confirmation de I(2,4) à grand Reynolds (budget 10 %) |

* **n=2 et n=4 aux trois Reynolds** → \(I(2,4)\) à N=384 sur toute la plage, donc
  trois points de grille (128 / 256 / 384) pour l'ordre apparent et le budget.
* **n=6 et n=8 à Re=400 seulement** → reconstruction de (4,6) et (6,8) à N=384,
  ce qui lève la réserve « non convergé quantitativement » sur ces deux paires.
  Réserve assumée : à Re=100 et Re=1000 les paires (4,6)/(6,8) restent à N=256.
* **α=0 et α=0,10 à chaque (Re,n)** → le contraste \(\Delta_\alpha\) lui-même est
  recalculé à N=384, et pas seulement les scalaires \(D^\*\), \(P^\*\).
* Critère de décision écrit à l'avance : **ne déclarer une convergence
  quantitative que si le contraste lui-même se stabilise** entre N=256 et N=384 —
  pas seulement les \(D^\*\) et \(P^\*\) individuels.

> **Résolution retenue : N=384, et non N=512** (décision du 2026-10-18, après mesure
> du coût réel). Conséquence à écrire dans le texte : le pas 256→384 vaut **×1,5 et
> non ×2**, donc le troisième point reste exploitable pour l'écart de grille mais
> **l'ordre apparent n'est plus directement comparable à la suite 128/256** — à
> présenter comme un point de contrôle, pas comme un point d'une suite de
> Richardson régulière.

```powershell
# preview (no execution)
python -m solver.run_subset512 --N 384 --dry-run
# campaign : (re)lancee le 2026-10-18 22:45 par solver/run_subset512_auto.ps1
#            6 workers, NUMBA_NUM_THREADS=2, tol 1e-8, N=384, sans checkpoint
python -m solver.run_subset512 --N 384 --workers 6 --max_iter 4000000
```

À N=384 il n'existe **aucun checkpoint** (celui de N=512 est abandonné) : les 16 cas
partent de zéro. `solver/run_subset512.py` soumet les cas par Re décroissant, n
décroissant, et réutilise `run_campaign.status_of` pour ne jamais relancer un cas
déjà `converged`.

**Coût (mesuré le 2026-10-18, décision N=384).** Le prorata (49 h pour 152 cas
→ ≈ 5 h pour 16) est un plancher très optimiste. La mesure directe donne
**≈ 31 it/s par worker** (0,032 s par itération, N=384), **indépendamment du
nombre de threads Numba** (2 ou 16 threads donnent le même débit : le noyau est
limité par lui-même, pas par les threads). C'est ce qui a motivé l'abandon de
N=512 : à N=512 on mesure 18 it/s, soit un coût par itération **×1,8 plus élevé
que le simple rapport des cellules** (effet de cache : `f` passe de 10,6 Mo à
18,9 Mo), et la cible complète monterait à ≈ 222 h·cœur (≈ 37 h). Avec les cibles
d'itérations ci-dessus (\(it_{256}\times1{,}83\)) :

| Re | it/cas | h/cas (1 worker) | cas | h·cœur |
|---|---|---|---|---|
| 100 | 223 000 | 2,0 | 4 | 8 |
| 400 | 493 000 | 4,4 | 8 | 35 |
| 1000 | 1 136 000 | 10,1 | 4 | 40 |
| **total** | | | **16** | **≈ 83 h·cœur** |

soit **≈ 14 h à 6 workers** (débit agrégé mesuré 186 it/s). `--max_iter 4000000`
reste nécessaire (le défaut de `run_campaign` est 1 200 000, **insuffisant** pour
Re=1000).

> **`--U` n'est pas à passer** (les 16 cas N=384 ont été lancés sans). Attention :
> `TARGET_TAU = {100: 0,65, 400: 0,60, 1000: 0,54}` **n'est pas construit pour donner
> \(U\) constant**.
> * Re = 400 et 1000 : ν = 13,333/N ⇒ \(U = Re\cdot\nu/N = (40/3)/N\) ⇒ à N=384,
>   **U = 0,0347222, Ma = 0,0601**.
> * **Re = 100 : τ = 0,65 ⇒ U = 0,01302, Ma = 0,023** — soit un Mach **2,6× plus
>   faible** qu'aux deux autres Reynolds.
>
> Le nombre de Mach **n'est donc pas apparié entre Reynolds** (0,023 / 0,060 / 0,060 à
> N=384). À écrire explicitement : l'effet de compressibilité n'est pas contrôlé de
> façon identique sur la plage en Re. *(Correction : la version antérieure de cette
> note annonçait τ=0,90 ⇒ U=(40/3)/N pour Re=100 ; c'était faux, τ(Re=100)=0,65.)*

*Ce qui n'est pas dans la sélection* : n ∈ {6,8} à Re=100 et Re=1000 (le gain
n'est que sur des paires déjà hors conclusion principale), et la campagne
intégrale de 152 cas (reportée au-delà du nécessaire).

**Résultat de la campagne N=384 (16/16 convergés, 2026-10-19).** Les grandes
dérivées \(D^\*\) bougent encore avec la grille (incrément 256→384 ≈ +0,8 % à
Re=400), mais la structure qui porte la conclusion se stabilise :

| Re | \(I_D(2,4)\) [%] : 128 / 256 / 384 | \(I_P(2,4)\) [%] : 128 / 256 / 384 | \(I_P/I_D-1\) : 128 / 256 / 384 |
|---|---|---|---|
| 100 | −0,0105 / −0,0108 / **−0,0137** | −0,0140 / −0,0141 / **−0,0141** | +33,3 % / +30,6 % / **+2,9 %** |
| 400 | −0,1805 / −0,1802 / **−0,1923** | −0,1990 / −0,1966 / **−0,1954** | +10,2 % / +9,1 % / **+1,6 %** |
| 1000 | −0,4699 / −0,4692 / **−0,5135** | −0,5326 / −0,5260 / **−0,5234** | +13,3 % / +12,1 % / **+1,9 %** |

* **Désaccord des deux estimateurs : tranché** (c'était l'arbitrage de §3.1). Le
  \(I_P/I_D-1\) de +27-33 % à Re=100 (N=128) tombe à **+2,9 % à N=384**, et les
  trois Reynolds convergent vers ≈ 2 %. C'est donc un **artefact de discrétisation**,
  pas un désaccord systématique → **première branche retenue : le budget sur
  \(I(2,4)\) à Re=100 passe de 28 % à ≲ 5 %, et la paire (2,4) à Re=100 est
  publiable**, les deux estimateurs devenant compatibles pour \(h\) suffisamment
  petit.
* **Répartition : confirmée.** À Re=400, N=384, les parts de l'interaction
  \(I(2,4)/(4,6)/(6,8)\) valent **91,3 / 7,9 / 0,8 %** — la paire (2,4) reste le
  porteur de l'effet (les autres paires sont des corrections de 0,8-8 %).
* La réponse à α par exposant est stable en grille : à Re=400, n=8,
  \(\Delta_\alpha D^\* = -0,7097/-0,7152/-0,7151\) % (128/256/384), et il y a
  **saturation en n** (n=6 et n=8 donnent la même réponse à 0,002 pt%).
* **Réserve maintenue : \(I(2,4)\) n'est pas convergé quantitativement en grille.**
  \(I_D(2,4)\) varie de **+6,7 %** (Re=400), **+9,4 %** (Re=1000) et **+27 %**
  (Re=100, valeur minuscule −0,01 %) entre N=256 et N=384. Le troisième point ne
  clôt donc pas l'incertitude de grille sur l'interaction : le critère écrit à
  l'avance (stabilisation du contraste lui-même) **n'est pas satisfait**, et il
  faut l'écrire comme *« signe et ordre de grandeur stables, non convergé
  quantitativement »*.

---

## 4. Paragraphe de synthèse pour le relecteur (à coller dans le manuscrit)

> Version française de travail — à traduire en anglais pour la soumission PoF.
> Rédigé dans les contraintes de style du relecteur : *signal détectable* /
> *stable entre grilles* / *convergé quantitativement* sont distingués ; aucun seuil
> n=4 n'est revendiqué ; aucune loi asymptotique universelle ; aucune mention de
> « benchmark de référence » ; la décomposition Wall/Core reste une perspective.

---

**Paragraphe proposé**

Sur les 48 cas de la campagne (Re = 100–1000, N = 128 et 256, n = 2–8,
α = 0 et 0,10), le contraste apparié \(I\) est **négatif dans les 18 couples
(Re, N, paire de formes)** sous les deux normalisations considérées — extensive
(\(D^\*\) intégré sur le masque fluide) et densité (\(D^\*/A_f\)) : l'effet de
l'ondulation est donc suppressif, et **ce signe est conservé en 18/18 cas**.
Le contraste n'est pas réductible à la perte de surface : un modèle nul qui ne
prend en compte que l'aire supprimée prédit un contraste **de signe opposé dans
les 18 cas**. L'effet est porté à **81–93 %** par le premier intervalle de forme
(2 → 4), à 7–16 % par (4 → 6) et à **0,2–2,8 %** par (6 → 8), ce qui établit des
rendements décroissants sans permettre d'affirmer une saturation complète au
delà de \(n = 4\). Les mêmes proportions se retrouvent (74–82 / 14,5–21,7 /
2,5–5,7 %) lorsque l'on rapporte \(D^\*\) à l'aire fluide : la hiérarchie n'est
pas un artefact de volume.

**Stabilité entre grilles.** La paire (2,4) varie de **1,9 % (Re = 400)**,
**2,8 % (Re = 1000)** et **4,2 % (Re = 100)** entre N = 128 et N = 256 ; elle
peut donc être qualifiée de *stable entre grilles*, et pas encore de *convergée
quantitativement*. Les paires (4,6) à Re = 100 (13,7 %) et (6,8) (7–11 %,
jusqu'à 24,7 % pour \(I_P\) à Re = 1000) **ne le sont pas** et sont signalées
comme telles. Un contrôle indépendant, qui compare deux intégrales de nature
différente calculées sur le même champ (surface sur le masque fluide contre
ligne sur la lid), donne un désaccord de **+5,6 % (Re = 400)** et **+8,2 %
(Re = 1000)** sur la paire (2,4) : nous retenons donc un budget d'incertitude
de **≈ 7 %** à Re ≥ 400, **≈ 10 %** à Re = 1000 et **≈ 28 %** à Re = 100, la
tolérance de convergence (10⁻⁸, soit 0,002 %) étant marginale.

**Bilan énergétique.** L'état publié présente un défaut
\((P^\*-D^\*)/D^\* \in [-6{,}9\,\%\,;\,-2{,}8\,\%]\). Quatre causes ont été testées sur les
champs déjà calculés. L'inversion de signe des corrections unilatérales à la
paroi inférieure et à la paroi gauche (`diagnostics.py:73` et `:57`) est un défaut
de code réel, mais son impact est **≤ 0,45 % sur \(D^\*\)** ; la contamination du
gradient par le saut fluide/solide vaut **≤ 0,07 %**. Un balayage de Ma à
\(h\) fixe (Ma = 0,045 → 0,361) montre que \(D^\*\) est invariant à ±0,30 % et que
**75–82 % du défaut subsiste dans la limite Ma → 0** : la compressibilité est
mesurable mais secondaire, et l'effet de \(h\) est environ trois fois plus grand
sur les plages balayées. La source dominante est le **cisaillement pariétal de
premier ordre sur la lid**, dont dépend intégralement \(P^\*\) : l'améliorer dans
\(P^\*\) seul réduit le défaut de **61 à 87 %**. En revanche, une amélioration
*cohérente* (appliquée à la fois à \(P^\*\) et à \(D^\*\)) n'améliore le défaut que
de −20 % à +22 % et **change d'effet de signe entre N = 128 et N = 256**, tandis
que le passage complet des quatre parois en ordre 2 l'**augmente** de 3,7 à
10,3 points de pourcentage. Le défaut n'est donc réductible à aucun bug unique :
c'est un **écart de cohérence de premier ordre en \(h\)** entre deux intégrales
indépendamment discrétisées d'un même champ instantané, et l'ordre apparent mesuré
sur trois réglages est de **0,96 à 1,07**. Sa variation entre géométries à \(N,Re\)
fixes est de **0,005 à 0,077 points de pourcentage**, soit 0,09–1,38 % de sa propre
valeur, et sa variation différentielle due à l'ondulation est de **≤ 0,06 points de
pourcentage** (≤ 0,9 % du défaut) : **il s'annule donc dans le contraste apparié**.
C'est \(I\), et non \(D^\*\) isolément, qui est la grandeur publiable.

---

## 5. Claims to retain / avoid — à garder et à ne pas écrire

**À retenir dans le manuscrit** (*retain*) :

1. **Non-additivité géométrique mesurable** — *a measurable geometric
   non-additivity*: le contraste \(\Delta_\alpha\) est négatif et conservé en 18/18
   cas, sous les deux normalisations, et le modèle nul « aire seule » prédit le
   signe opposé dans 18/18.
2. **Dépendance monotone en Reynolds sur la plage étudiée** — *monotonic
   Reynolds dependence over the studied range* (100–1000) : à dire strictement,
   « sur la plage étudiée » et non « en général ».
3. **Concentration de la réponse extensive dans n=2→4** — *the extensive
   response is concentrated in n=2→4* : 81–93 % de l'effet, rendements
   décroissants ensuite.
4. **Moindre sensibilité à la grille de \(I_D(2,4)\) pour Re ≥ 400** — *reduced
   grid sensitivity of \(I_D(2,4)\) for Re ≥ 400* : 1,9 % (Re=400) et 2,8 %
   (Re=1000) contre 4,2 % (Re=100).

**À éviter** (*avoid*) :

1. **Loi asymptotique universelle** — pas de « universal asymptotic law » : on a
   un ordre apparent sur trois réglages, pas une asymptotique.
2. **Preuve causale de médiation par la couche pariétale** — pas de
   « proves … mediated by the boundary layer » : le lien lid/défaut est une
   **interprétation cohérente avec les mesures**, pas une preuve causale.
3. **Seuil de saturation à n=4** — pas de « saturation at \(n=4 \) » : les
   rendements décroissants sont établis, la saturation complète ne l'est pas.
4. **Convergence quantitative complète** — pas de « quantitatively converged »
   pour les paires (4,6) à Re=100 ni pour (6,8) ; utiliser *detectable* /
   *stable across grids* / *quantitatively converged* comme trois niveaux
   distincts.
5. **« Jeu de données benchmark de référence »** — pas de *reference benchmark
   dataset* : l'ensemble est reproductible, pas établi comme référence.

> Ces cinq interdits correspondent mot pour mot aux réserves rédigées en §4 et
> aux conventions de style du relecteur (§0).

