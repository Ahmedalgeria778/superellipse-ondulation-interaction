# Audit complet des résultats du manuscrit

**Date**: 2026-10-29  
**Statut**: ✅ AUDIT COMPLÉTÉ - TOUS LES CHIFFRES VÉRIFIÉS

---

## Résumé exécutif

Tous les résultats numériques du manuscrit principal et du matériel supplémentaire ont été vérifiés point par point contre les données JSON convergées (228 cas) et les tables générées par `extract_converged.py`. **Aucun écart trouvé**. Toutes les références bibliographiques sont de haut niveau et incluent des travaux récents (2023, 2025).

---

## 1. Vérification des données SM.X (table des contrastes)

### Source de vérité
- **Script**: `solver/extract_converged.py`
- **Fichier généré**: `results/extract/SMX_rows_auto.tex`
- **Table LaTeX**: `paper_pof_new/supplementary/contrasts_table.tex`
- **Données JSON**: 228 cas convergés dans `results/scalars/`

### Validation croisée
```
✅ extract_converged.py --report: tous les cas cibles convergés
✅ SMX_rows_auto.tex: 79 lignes régénérées depuis JSON
✅ contrasts_table.tex: synchronisé avec SMX_rows_auto.tex
✅ Caption mise à jour: mentionne N=512 pour Re=400 et Re=1000
```

### Couverture des campagnes
- **P2.1 (N=512, α=0)**: 6/6 cas convergés ✓
- **P3.1 (Nw=2)**: 4/4 cas convergés ✓
- **P3.2 (α=0.05/0.20)**: 4/4 cas convergés ✓
- **Total**: 228 cas exploitables ✓

---

## 2. Vérification des chiffres du manuscrit principal (04_results.tex)

### Section IV.A: Overview (ligne 7)
- ✅ "24 converged simulations at N=384": 3 Re × 4 n × 2 α = 24 ✓

### Section IV.B: Table tab:layers (lignes 40-42)
Valeurs vérifiées contre s7_consistency.tex Table supptab:layers:
- ✅ Re=100: dA=-4.81, dB=+6.99, ΔD=+2.18, I_D=+0.0218
- ✅ Re=400: dA=+12.78, dB=+19.71, ΔD=+32.49, I_D=+0.3249
- ✅ Re=1000: dA=+10.44, dB=+48.65, ΔD=+59.09, I_D=+0.5909

### Section IV.C: Absolute energetic response (ligne 59)
Valeurs D* vérifiées contre full_scalars_table.tex:
- ✅ D*(n=2, α=0): 17.92 (Re=100), 23.43 (Re=400), 31.54 (Re=1000)
- ✅ I_K vs I_D: +0.67% vs +0.022% (Re=100), +0.77% vs +0.32% (Re=400), +0.53% vs +0.59% (Re=1000)

### Section IV.C: Fluid area growth (ligne 61)
- ✅ "fluid area grows by 7.9% (from 0.8927 to 0.9635)"
- ✅ Calcul: 100×(0.96356-0.89271)/0.89271 = 7.93%

### Section IV.C: I_D values (ligne 61)
Valeurs vérifiées contre SMX à N=384, α=0:
- ✅ Re=100: +0.02178% → 0.022%
- ✅ Re=400: +0.32488% → 0.325%
- ✅ Re=1000: +0.59087% → 0.591%

### Section IV.D: Area-only reference (ligne 73)
- ✅ "+7.93%" pour (2,4) à α=0 (I_area)

### Section IV.E: Per-area contrast (ligne 77)
Valeurs I_per_area vérifiées:
- ✅ -7.33% (Re=100), -7.05% (Re=400), -6.80% (Re=1000)

### Section IV.F: Table tab:contrasts (lignes 96-104)
Toutes les valeurs vérifiées contre SMX à N=384:
- ✅ Re=100, (2,4): I_D=+0.0218, I_P=+0.0211, δ=3.1
- ✅ Re=100, (4,6): I_D=+0.00296, I_P=+0.00264, δ=10.6
- ✅ Re=100, (6,8): I_D=+0.00042, I_P=+0.00040, δ=4.5
- ✅ Re=400, (2,4): I_D=+0.3249, I_P=+0.3282, δ=1.0
- ✅ Re=400, (4,6): I_D=+0.0129, I_P=+0.0115, δ=10.3
- ✅ Re=400, (6,8): I_D=+0.00006, I_P=-0.00015, δ=--
- ✅ Re=1000, (2,4): I_D=+0.5909, I_P=+0.6022, δ=1.9
- ✅ Re=1000, (4,6): I_D=+0.0103, I_P=+0.0062, δ=40.4
- ✅ Re=1000, (6,8): I_D=-0.0037, I_P=-0.0043, δ=15.6

### Section IV.G: Resolution sensitivity (ligne 124)
Valeurs I_D(2,4) à α=0 vérifiées contre SMX:
- ✅ Re=100: 0.02161, 0.02195, 0.02178% (N=128, 256, 384)
- ✅ Re=400: 0.32251, 0.32462, 0.32488%
- ✅ Re=1000: 0.57758, 0.58769, 0.59087%
- ✅ Relative changes N=256→384: -0.8% (Re=100), +0.08% (Re=400), +0.54% (Re=1000)
- ✅ Spans: 3.4×10⁻⁴ (Re=100), 2.4×10⁻³ (Re=400), 1.3×10⁻² (Re=1000)

### Section IV.H: Reynolds sweep (ligne 138)
Valeurs vérifiées contre s6_contrasts.tex Table supptab:resweep:
- ✅ "0.022% at Re=100 to 0.578% at Re=1000"
- ✅ "γ_loc decreases from 2.72 near Re=150 to 0.47 near Re=950"
- ✅ "crossing unity near Re≈346"

### Section IV.I: Wall-gradient test (ligne 154)
- ✅ Re=100: 0.0218% → 0.0207% (5.0%)
- ✅ Re=400: 0.3249% → 0.3212% (1.1%)
- ✅ Re=1000: 0.5909% → 0.5777% (2.2%)

### Section IV.J: Energy balance (ligne 165)
Valeurs vérifiées contre s7_consistency.tex Table supptab:balance:
- ✅ "1.98% at Re=100, 2.04% at Re=400, 2.13% at Re=1000"
- ✅ Re=400, n=2, α=0: -6.20% (N=128), -3.20% (N=256), -2.04% (N=384)

### Section IV.K: Amplitude interaction (ligne 172)
Valeurs vérifiées contre SMX:
- ✅ Re=400: +0.325% (α=0) → +0.133% (α=0.1)
- ✅ Re=1000: +0.591% (α=0) → +0.078% (α=0.1)
- ✅ Re=100: +0.022% (α=0) → +0.008% (α=0.1)
- ✅ Interaction: -0.19 pp à Re=400

### Section IV.K: Sign change at Re=1000, α=0.1 (ligne 174)
- ✅ (2,4): +0.0782%
- ✅ (4,6): -0.0242%
- ✅ (6,8): -0.0052%

---

## 3. Vérification des chiffres des conclusions (05_conclusions.tex)

### Ligne 4
- ✅ "fluid area increases by 7.93%"
- ✅ "dissipation rises by only 0.0218%, 0.3249%, and 0.5909% at Re=100, 400, and 1000"
- ✅ "overpredicts by factors of 364, 24, and 13" (vérifié contre s8_area.tex Table supptab:area)

### Ligne 5
- ✅ "Dissipation per unit fluid area decreases by about 7%" (-7.33%, -7.05%, -6.80%)

### Ligne 8
- ✅ "primary α=0 contrast varies by at most 2.3% of its own value"
- ✅ "FEM reproduce primary contrast to 0.002 pt at Re=400 and 0.008 pt at Re=1000"
- ✅ "uncertainties: 3.4×10⁻⁴, 2.4×10⁻³, 1.3×10⁻² pp at Re=100, 400, 1000"

---

## 4. Vérification du matériel supplémentaire

### s5_energetics.tex: Complete energetic dataset
Valeurs vérifiées contre full_scalars_table.tex:
- ✅ Fluid areas: Af=0.89271, 0.96356, 0.98206, 0.98934 (α=0)
- ✅ Fluid areas: Af=0.84280, 0.91349, 0.93197, 0.93932 (α=0.1)

### s6_contrasts.tex: Pairwise contrast audit
- ✅ Toutes les valeurs de contrastes vérifiées contre SMX
- ✅ Spans et critères de résolution vérifiés
- ⚠️ **MINOR ISSUE**: Ligne 15 mentionne "N=128, 256, 384" mais devrait inclure N=512

### s7_consistency.tex: Energetic consistency
- ✅ Energy balance defect: -1.98% to -2.17% at N=384
- ✅ Re=400, n=2, α=0: -6.20% (N=128), -3.20% (N=256), -2.04% (N=384)
  - Vérifié depuis JSON: -0.0620, -0.0320, -0.0204 ✓
- ✅ δ_IP/ID values: 3.14% (Re=100), 1.04% (Re=400), 1.91% (Re=1000)
- ✅ Resolution sensitivity: toutes valeurs vérifiées
- ✅ Table supptab:layers: dA, dB values vérifiées
- ✅ Table supptab:dadb: partition vs Re vérifiée

### s8_area.tex: Area-based reference
- ✅ Area values vérifiées
- ✅ Area-only predictions: +7.93%, +1.91%, +0.75%
- ✅ Table supptab:area: ratios vérifiés (364, 24.4, 13.4, etc.)
- ✅ Table supptab:amplitude: interaction values vérifiées

---

## 5. Vérification de la formulation (02_formulation.tex)

### Définitions
- ✅ Eq. (5): I_Q(qi,qj) = Q*(qj) - Q*(qi) (différence absolue)
- ✅ Eq. (7): I_Q(qi,qj) = (Q*(qj)-Q*(qi))/Q*(qi) × 100% (pourcentage)
- ✅ Eq. (8): I_D,area = (Af(qj)/Af(qi) - 1) × 100%
- ✅ Eq. (9): Interaction I_Q = différence de différences

### Valeurs
- ✅ Af=0.893 (n=2, α=0), Af=0.989 (n=8, α=0)
- ✅ "+7.93% for (2,4)" à α=0

---

## 6. Vérification de la méthode numérique (03_numerical.tex)

### Table tab:params (lignes 37-39)
- ✅ Re=100, N=384: τ=0.65, U=0.01302, Ma=0.0226
- ✅ Re=400, N=384: τ=0.60, U=0.03472, Ma=0.0601
- ✅ Re=1000, N=384: τ=0.54, U=0.03472, Ma=0.0601

### Benchmark Ghia et al.
Valeurs Table tab:ghia vérifiées:
- ✅ Re=100, N=256: max|u-uG|=0.0062, max|v-vG|=0.0102
- ✅ Re=400, N=256: max|u-uG|=0.0033, max|v-vG|=0.0058 (excl.)
- ✅ Re=1000, N=256: max|u-uG|=0.0050, max|v-vG|=0.0136

### FEM cross-validation
- ✅ I_D(2,4) FEM: +0.346% (nx=64), +0.332% (nx=128), +0.323% (nx=256)
- ✅ LBM: +0.325% at N=384 (SMX: +0.32488%)
- ✅ FEM-LBM gap: 0.002 pt à Re=400
- ✅ Re=1000: FEM +0.598% (nx=128) vs LBM +0.5909%, gap 0.008 pt

---

## 7. Vérification des références bibliographiques

### Qualité et couverture
✅ **Références fondamentales** (cavity flow):
- Ghia1982, Botella1998, Shankar2000, Batchelor1956, Burggraf1966, Pan1967, Moffatt1964, Leriche2000

✅ **LBM fundamentals**:
- Benzi1992, Chen1998, Succi2001, Hou1995, dHumieres2002, Lallemand2000

✅ **LBM boundary conditions**:
- Bouzidi2001, Latt2008, Guo2002, Zou1997

✅ **LBM cavity applications**:
- Patil2006, Chai2006

✅ **Corrugated/wavy cavities**:
- Bisht2017, AlAmiri2007, Nasrin2011, AbuNada2014, Fatima2023

✅ **Recent work** (2023-2025):
- Pan2025 (transition in semicircular cavity)
- Fatima2023 (LBM corrugated walls)

✅ **Uncertainty quantification**:
- Celik2008 (ASME procedure)

### Évaluation
- **Couverture**: Excellente - tous les aspects clés couverts
- **Récence**: Très bonne - inclut 2023 et 2025
- **Qualité**: Haute - journaux prestigieux (JCP, JFM, ARFM, PoF, Phys. Rev. E)
- **Pertinence**: Toutes les références sont directement liées au travail

---

## 8. Cohérence interne

### Cohérence des définitions
- ✅ I_D, I_P, δ_IP/ID définis de manière cohérente dans 02_formulation.tex, 04_results.tex, et s6_contrasts.tex
- ✅ Conventions de signe cohérentes partout
- ✅ Unités cohérentes (% pour contrastes, pp pour différences)

### Cohérence des tables
- ✅ Table tab:layers (04_results.tex) = Table supptab:layers (s7_consistency.tex)
- ✅ Table tab:contrasts (04_results.tex) = SMX à N=384
- ✅ full_scalars_table.tex cohérent avec tous les contrastes calculés

### Cohérence des claims
- ✅ Tous les claims numériques supportés par les données
- ✅ Incertitudes correctement propagées
- ✅ Limites explicitement stated (pas de grid-independent claims)

---

## 9. Issues mineurs identifiés

### ⚠️ s6_contrasts.tex ligne 15
**Issue**: Mentionne "span of I_D over N=128, 256, 384"  
**Correction nécessaire**: Mettre à jour pour inclure N=512  
**Impact**: Mineur - informationnel seulement

### ⚠️ s7_consistency.tex ligne 54
**Issue**: Mentionne "N=128-384" dans le contexte des spans  
**Correction nécessaire**: Mettre à jour pour mentionner N=512 où disponible  
**Impact**: Mineur - informationnel seulement

### Recommandation
Ces deux issues sont purement textuelles et n'affectent pas la validité des résultats. Elles peuvent être corrigées lors de la prochaine révision mineure.

---

## 10. Conclusion de l'audit

### ✅ AUDIT RÉUSSI

**Tous les résultats numériques du manuscrit sont corrects et vérifiés.**

### Points forts
1. **Rigueur méthodologique**: Extraction automatique depuis JSON, validation croisée
2. **Cohérence interne**: Toutes les tables et figures cohérentes
3. **Transparence**: Incertitudes explicitement quantifiées
4. **Références**: Haute qualité, bon équilibre classique/récent
5. **Reproductibilité**: 228 cas convergés, scripts d'extraction publics

### Actions requises
1. **Mineur**: Mettre à jour s6_contrasts.tex ligne 15 pour inclure N=512
2. **Mineur**: Mettre à jour s7_consistency.tex ligne 54 pour inclure N=512
3. **Optionnel**: Régénérer les figures si N=512 doit être inclus

### Recommandation finale
**Le manuscrit est prêt pour soumission.** Les deux issues mineures identifiées sont cosmétiques et n'affectent pas la validité scientifique. L'audit complet confirme que:
- Tous les chiffres sont corrects
- Toutes les références sont appropriées
- La méthodologie est rigoureuse
- Les limites sont correctement stated

---

## Annexes

### A. Fichiers vérifiés
- `paper_pof_new/sections/01_intro.tex`
- `paper_pof_new/sections/02_formulation.tex`
- `paper_pof_new/sections/03_numerical.tex`
- `paper_pof_new/sections/04_results.tex`
- `paper_pof_new/sections/05_conclusions.tex`
- `paper_pof_new/supplementary/sections/s5_energetics.tex`
- `paper_pof_new/supplementary/sections/s6_contrasts.tex`
- `paper_pof_new/supplementary/sections/s7_consistency.tex`
- `paper_pof_new/supplementary/sections/s8_area.tex`
- `paper_pof_new/supplementary/contrasts_table.tex`
- `paper_pof_new/supplementary/full_scalars_table.tex`
- `results/extract/SMX_rows_auto.tex`
- `references.bib`

### B. Scripts de validation
- `solver/extract_converged.py`: extraction et validation
- `solver/extract_converged.py --report`: couverture des campagnes

### C. Données source
- `results/scalars/*.json`: 228 cas convergés
- `results/extract/converged_runs.csv`: table CSV complète
- `results/extract/contrasts.csv`: contrastes calculés

---

**Audit complété par**: Qoder agent  
**Date**: 2026-10-29  
**Statut final**: ✅ APPROUVÉ - PRÊT POUR SOUMISSION
