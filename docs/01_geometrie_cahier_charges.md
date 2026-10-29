# 01 — Cahier des charges géométrique (Option A validée)

L'option retenue est : **fond superelliptique + boîte + couvercle plat horizontal**.
Domaine : `D = { (x,y) ∈ R^2 : |x| <= a, y_b(x; n, alpha, Nw) <= y <= H }`.
Parois : latérales verticales fixes `x = ±a` ; fond courbe mobile dans la superposition (`y_b`), lid plat mobile `y = H`.

## 1) Fond de base — superellipse

```
y_b0(x; n) = -b [1 - (|x|/a)^n]^(1/n),   x ∈ [-a, a],   b > 0.
```

Propriétés :
- Épaule à `y=0` en `x=±a` : `y_b0(±a;n)=0` pour tout `n>=1`.
- Profondeur minimale au centre : `y_b0(0;n)= -b`.
- Forme : `n=2` demi-ellipse/demi-cercle (si a=b) ; `n>2` coins plus "plats" à base, fond plus concave en pente ; `n->∞` limite rectangulaire `[ -a,a ] × [-b,0]`.

Paramètres de référence choisis (cohérents avec Projet B) : `a = L/2`, `b` choisi pour épaisseur/rapport (valeur fixée dans campagne : `b/H` constant, e.g. `b=0.5` ou tel que profondeur cohérente ; à préciser dans `geometry/configs.py`).

## 2) Ondulation locale du fond

L'ondulation `eta(x)` est superposée **au fond superelliptique** uniquement :

```
y_b(x; n, alpha, Nw) = y_b0(x; n) + eta(x; alpha, Nw).
```

Critère de continuité/fermeture imposé : **`eta(±a; alpha,Nw) = 0`** pour que le fond rejoint naturellement les parois latérales verticales à `y=0, x=±a` (pas de discontinuité artificielle).

Forme choisie (reprise exacte de Projet A, normalisée sur largeur `L=2a`) :

```
s(x) = (x + a) / (2a)  ∈ [0,1],   x ∈ [-a, a]
eta(x; alpha, Nw) = (alpha/2) [1 - cos(2 pi Nw * s(x))].
```

Donc :
- `eta(-a)= (alpha/2)(1 - cos(0)) = 0`
- `eta(+a)= (alpha/2)(1 - cos(2 pi Nw)) = 0`
- Amplitude crête-à-crête liée à `alpha` (adimensionnée par `H` dans Projet A ; ici aussi par `H`).
- `alpha >= 0`, entier `Nw >= 1`.

Cette forme préserve les bords et évite d'introduire une marche artificielle.

## 3) Domaine complet et CL

- `x ∈ [-a, a]` (largeur `L=2a`)
- `y ∈ [y_b(x), H]` (hauteur verticale totale fixe `H` entre fond effectif et lid plat)
- **Couvercle plat** `y = H`, `x\in[-a,a]` : vitesse imposée `u_w = (U_ref, 0)` (lid entraîné uniformément). Choix **délibéré** (différent du lid tangent Projet B) pour isoler l'interaction forme+ondulation avec forçage identique.
- **Parois latérales verticales** `x = ±a`, `y \in [0, H]` (segment au-dessus de l'épaule `y_b(±a)=0`) : fixes (no-slip). **Fond** `y = y_b(x)$` : no-slip.
- **Géométrie "boîte + relief"** : le volume fluide est une cavité rectangulaire globale `[ -a,a ]×[0,H]` amputée du relief sous `y_b(x)$` (et complétée par l'ondulation). Cette construction facilite le maillage cartésien + traitement Bouzidi (frontière immergée/courbe) tel que pratiqué dans A et B.

## 4) Paramétrage normalisé

Adimensionnement cohérent (Projet A) :
- Longueurs par `H` (ou par `a$` selon besoin). On fixe `H = 1` dans les unités adimensionnelles de calcul.
- Vitesses par `U_ref$`.
- `Re = U_ref H / ν$`.

Choix de campagne pilote (validés) :
- `a/H = 0.5$` (cavité carrée globale `L×H = 1×1$` si on centre ; ou `L=1, H=1$` avec `x\in[0,1]$`). À fixer explicitement dans `configs.py` (convention axe centrée ou translatée). La convention centrée `x\in[-0.5,0.5], H=1$` est naturelle pour la superellipse symétrique.
- `b/H$` : profondeur relative du fond superelliptique (fixe sur tous les cas). Valeur typique `0.5$` (fond à mi-hauteur) ou `b$` tel que la surface fluide reste comparable entre n. À choisir pour minimiser effet de volume (voir §5).
- `n \in \{4,8\}$` (cas pilotes), étendre si interaction significative.
- `(alpha,Nw) \in \{(0,1),(0.10,1)\}$`.
- `Re \in \{400\}$` (pilote), étendre {100,400,1000}.

## 5) Contraintes de validité

1. **Continuité C^0** du fond en `x=±a$` (imposée par `eta(±a)=0$, `y_b0(±a)=0$). Dérivabilité : `dy_b/dx (±a)$` peut être finie ou non selon `n$` (superellipse à coins plus ou moins nets) — c'est physique (forme globale).
2. **Non-recouvrement** : exiger `y_b(x) < H$` sur tout l'intervalle (toujours vrai si `b>=0, alpha$ raisonnable).
3. **Monotonie/forme admissible** : `alpha$` limité pour éviter auto-intersection locale (ondulation trop forte par rapport à la pente locale). Vérification géométrique dans `check_geometry.py` (minimum de `y_b(x)- ?$`, pente, courbure).
4. **Comparabilité entre cas** : pour l'interaction `I_Q$, il est préférable d'avoir **surface fluide `A_f = ∫_{-a}^{a} (H - y_b(x)) dx$** raisonnablement proche entre `(n,alpha)$` pour un même `n$` (effet de volume contrôlé). On pourra documenter `A_f(n,alpha)$` (calcul analytique/numérique exact) et, le cas échéant, discuter de l'effet de confinement. Objectif : variation < ~2% si possible, sinon corriger interprétation.
5. **Traitement LBM (Bouzidi)** : frontière immergée `y_b(x)$` doit être bien résolue (distance aux noeuds fluides/solides) ; contrôle de résolution `N=256,384,512$`.

## 6) Implémentation (modules)

- `geometry/superellipse.py` — fonction `y_b0(x,n,a,b)$`, dérivées premières (facultatives), aire/longueur (réutilisé Projet B).
- `geometry/wall_functions.py` — `eta(x,alpha,Nw,a)$`, normalisation, vérifs.
- `geometry/combined_geometry.py` — `y_b(x,n,alpha,Nw,a,b)$`, `dy_b_dx(x,...)$` (numérique stable), enveloppe, masques.
- `geometry/configs.py` — définitions des 4 configs pilotes + paramètres fixes (`a,b,H,U_ref,Re`).
- `geometry/check_geometry.py` — validation (graphiques + tests numériques : `y_b(±a)=0$, continuité, aire fluide, pente max, distance minimale). Génère figures dans `figures/geometry/`.
