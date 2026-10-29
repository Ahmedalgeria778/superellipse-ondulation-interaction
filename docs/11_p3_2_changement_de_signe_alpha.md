# 11 — P3.2 : le contraste (2,4) change de signe avec l'amplitude à Re=1000

Date : 2026-10-28. Campagne P3.2 terminée (4/4 cas convergés, N=384, Re=1000,
tol=1e-8). Ce résultat est **hors périmètre de la soumission initiale** ; il est
consigné ici pour ne pas être perdu, pas dans le manuscrit.

## Valeurs

Contraste (2,4) à Re=1000, N=384, tous les cas `status=converged` :

| α | D*_n2 | D*_n4 | I_D [%] | I_P [%] | δ_IP/ID [%] | I_Af [%] | contraste par aire [%] |
|---|---|---|---|---|---|---|---|
| 0.00 | 31.54390 | 31.73028 | **+0.5909** | +0.6022 | 1.91 | +7.936 | −6.805 |
| 0.05 | 31.39118 | 31.49260 | **+0.3231** | +0.3279 | 1.50 | +8.145 | −7.233 |
| 0.10 | 31.19820 | 31.22261 | **+0.0782** | +0.0797 | 1.87 | +8.388 | −7.667 |
| 0.20 | 30.82923 | 30.79350 | **−0.1159** | −0.1264 | 9.10 | +8.942 | −8.314 |

Itérations : 496 k / 1 206 k (α=0.05), 536 k / 1 118 k (α=0.20) ; α=0 et 0.10
proviennent de la campagne N=384 déjà publiée.

## Lecture

- `I_D(2,4)` décroît strictement avec α et **s'annule entre α=0.10 et α=0.20**.
  Interpolation linéaire locale : α* ≈ 0.14. Le pas n'est pas linéaire en α
  (−0.268 sur [0, 0.05], −0.245 sur [0.05, 0.10], −0.097 par 0.05 sur
  [0.10, 0.20]) : la décroissance ralentit, donc la vraie racine est plutôt
  au-delà de 0.14 — seul l'encadrement [0.10, 0.20] est solide.
- `I_P` suit `I_D` de près ; δ_IP/ID reste à 1.5–1.9 % tant que |I_D| est
  franchement non nul, et monte à 9.1 % à α=0.20 parce que |I_D| y est petit
  (0.116 %) pour un écart absolu quasi constant (~0.010 pp). C'est exactement
  la mise en garde de `sections/02_formulation.tex:102` sur les petits
  contrastes : rien d'anormal, mais à citer avec précaution.
- Le contraste par aire (`D*/A_f`) devient plus négatif avec α
  (−6.80 → −8.31 %), donc la normalisation extensive ne compense pas l'effet.
- Physique : l'ondulation ne fait pas que « réduire » le contraste (2,4), elle
  l'inverse. Combiné au fait déjà publié que (4,6) et (6,8) s'inversent dès
  α=0.1, cela suggère que l'actionnement de paroi inverse l'ordre énergétique
  des bas exposants au-delà d'une amplitude seuil.

## Vérification de cohérence avec le manuscrit (faite le 2026-10-28)

Aucune contradiction. Points exacts :

1. `sections/04_results.tex:172` : « at Re=1000 it falls from +0.591% to
   +0.078% » — comparaison α=0 → 0.1 seulement, toujours vraie.
2. `sections/04_results.tex:174` : « while (2,4) remains positive (+0.0782%) »
   — phrase portée explicitement par « At Re=1000 and α=0.1 », donc vraie.
   Elle ne doit pas être lue comme une propriété pour tout α ; la formulation a
   été resserrée depuis (voir « Portée ajoutée au corps du texte »).
3. Table SM.X (`supplementary/contrasts_table.tex`) : **aucune** ligne
   α=0.05 ni α=0.20 (27 lignes à 0.00, 27 à 0.10). Aucune valeur figée ne
   peut donc être contredite.
4. `sections/04_results.tex:84` : la hiérarchie |I(2,4)|≫|I(4,6)|≫|I(6,8)| est
   annoncée « at both amplitudes examined », soit 0 et 0.1 ; (4,6) et (6,8)
   n'existent pas à α=0.20, la phrase n'est pas testable au-delà.

## Portée ajoutée au corps du texte (faite le 2026-10-28)

`04_results.tex:174` disait « while (2,4) remains positive ». Remplacé par
« while $(2,4)$ is still positive **at this amplitude** ($+0.0782\%$) ».
Choix délibéré de ne pas écrire « see below for larger amplitudes » : le
manuscrit n'examine que α=0 et 0.1, donc un renvoi vers un « plus loin » aurait
pointé dans le vide. La portée est maintenant portée par la phrase elle-même.
Recherche « monoton » et « remains positive / positive at / all amplitudes » :
aucune autre affirmation générique en α dans `sections/` ni `supplementary/`.

Compilation après modification : `_bmain.pdf` = 26 pages, 0 erreur,
0 référence non définie, 0 overfull, aucun rerun demandé ; copié en `main.pdf`.

## Identité des fichiers (à revérifier avant toute mise à jour de SM.X)

- SM.X canonique du papier superellipse, et unique dans ce projet :
  `paper/supplementary/contrasts_table.tex`
  (label `supptab:contrasts`, 54 lignes, α=0.00 et 0.10 seulement).
- Aucun autre manuscrit (ni papier PINN-cavité) n'est présent dans ce
  workspace ; `ls -d paper*` ne renvoie que `paper_pof_new` et deux archives
  `.zip` d'anciennes versions.
- Piège de build déjà rencontré : le `.bbl` local s'appelle `_bmain.bbl`, donc
  la séquence est `pdflatex -jobname=_bmain main.tex`, `bibtex _bmain`,
  `pdflatex` ×2, puis copie de `_bmain.pdf` vers `main.pdf`. Sans le `bibtex
  _bmain`, toutes les citations sortent en `[?]` et le PDF tombe à 25 pages
  (la bibliographie disparaît), ce qui ressemble à une régression de contenu.

## À faire dans une passe ultérieure (pas maintenant)

- Si α est étendu au-delà de 0.1 : la fin de `04_results.tex:174`
  (« the sign reversal affects the trailing transitions ») devient incomplète,
  le (2,4) s'inversant aussi pour α ≳ 0.14 ; et la phrase « is still positive
  at this amplitude » devra être suivie de la valeur à l'amplitude supérieure.
- Table SM.XVI (interaction en α) reste hors soumission, à produire depuis
  `solver/audit_contrasts.py` — et non depuis les contrastes SM.X.
- Si une figure α est ajoutée, tracer I_D(2,4) contre α à Re=1000 avec les
  quatre points ci-dessus et marquer l'encadrement de la racine.

## Provenance

`results/scalars/se_n2_w0.05_Nw1_Re1000_N384_tol1e-08.json`,
`..._n4_...0.05...`, `..._n2_...0.20...`, `..._n4_...0.20...` ;
logs correspondants dans `results/logs/`. Chiffres reproductibles via
`python -m solver.extract_converged` (contrastes) une fois la table étendue.
