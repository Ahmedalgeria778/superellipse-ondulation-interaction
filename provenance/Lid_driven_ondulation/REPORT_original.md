# Study Report — Energy Budgets of a Lid-Driven Cavity with a Corrugated Floor

**Subject.** Verified quantitative study of the mechanical-energy budgets (injected
power `P*`, viscous dissipation `D*`, kinetic energy `K*`) of a cavity driven by a
moving lid, with a **corrugated bottom wall** and a **controlled lid profile**. Two
effects are distinguished: (A) the geometrical effect of the corrugation, and
(B) the kinematic effect of the velocity profile imposed on the lid.

**Data status.** Complete campaign: `results/scalar_table.json` = **78 rows, all
`converged` with `checks_ok=True`**. Resolution-convergence study
N=128/256/384/512 completed (`results/convergence/resolution.csv`, 18/18).
Protocol B completed (3/3). Computations stopped on 2026-10-16.

---

## 1. Method and LBM validation

### 1.1 Solver
- **D2Q9 MRT**, `c_s^2 = 1/3`, shear relaxation `omega = 1/tau`,
  viscosity `nu = c_s^2 (tau - 1/2)`, `Ma = sqrt(3) * U_ref`.
- Chapman–Enskog coefficient of the non-equilibrium tensor
  (`Pi^neq = -2 rho c_s^2 tau * S`, with **tau** and not tau-1/2): validated by a
  node-by-node fit on the flat cavity at Re=100, N=64; measured slope **-0.4355**
  vs `-2 rho c_s^2 tau = -0.4333` (0.5%), against -0.100 for `-2 rho nu`.
- **Boundary conditions**: the lid is imposed through wall populations; the
  corrugated bottom wall uses **Bouzidi interpolation**; corner correction
  `q_frac = 0.5`.

### 1.2 Geometries and lid profiles
- Corrugated floor parametrised by `(h/H, Nw)`; flat reference `F00`. Grid
  N in {128, 256, 384, 512}. Corrugated cases studied here: **G06** (h/H=0.05,
  Nw=2), **G12** (h/H=0.075, Nw=4), **G13** (h/H=0.10, Nw=1).
- **Uniform profile**: `U_lid/U_ref = 1` everywhere (mean = 1.000).
- **Regularised profile**: smoothstep `3t^2 - 2t^3` over a width `eps = 0.02` at
  both corners (exact mean = 1 - eps = **0.980**, minimum ~ 0 at the corners).
- **`smoothstep_mean1` profile** (protocol B): same shape renormalised to
  **mean = 1.000** (maximum = 1.0211).

### 1.3 Dimensionless quantities
```
K* = 0.5 * int |u*|^2 dA*                 (kinetic energy)
D* = 2 * int S*:S* dA*                    (dissipation, tensor S by finite differences)
P* = int U*_lid (du*/dy*)|lid dx*         (lid power)
balance = (P* - D*) / D*                  (residual; -> 0 at steady state)
```

### 1.4 Validation checks
- Reference flat cavity (`V01/V02/V03`, uniform lid): included in the campaign,
  automatic checks `checks_ok=True` for all **78 rows**.
- Mach check: `Ma_max_local <= 0.10` flagged per case. At `N=128, Ma=0.180`
  (`consistent` variant) the incompressible limit is exceeded, so those points
  are treated with caution.

### 1.5 Reynolds protocols and numerical variants
- `Re = U_ref * N / nu` (recomputed at each grid -> **`Re_eff` invariant**).
- **`consistent`**: `Ma` grows with `N` (Re=400: `Ma` = 0.180/0.090/0.060/0.045
  for N=128/256/384/512).
- **`ma_fixed`**: `U_ref = 0.05208` fixed, `Ma = 0.090`; `tau` = 0.55/0.60/0.65/0.70
  according to N.
- The two variants coincide at N=256 (identical parameters).

---

## 2. Absolute convergence

### 2.1 Quantities vs resolution (Re=400)

| Case | N | `P*` | `D*` | `K*` | balance |
|---|---|---|---|---|---|
| F00 (flat) | 128 | 21.4196 | 22.8341 | 0.03944 | -6.20% |
| F00 | 256 | 22.5825 | 23.3280 | 0.04012 | -3.20% |
| F00 | 512 | 23.2613 | 23.6047 | 0.04045 | -1.46% |
| G06 | 128 | 21.3834 | 22.7983 | 0.03932 | -6.21% |
| G06 | 256 | 22.5448 | 23.2904 | 0.03999 | -3.20% |
| G06 | 384 | 22.9930 | 23.4727 | 0.04021 | -2.04% |
| G06 | 512 | 23.2227 | 23.5662 | 0.04032 | -1.46% |
| G12 | 128 | 21.3214 | 22.7362 | 0.03900 | -6.22% |
| G12 | 256 | 22.4806 | 23.2266 | 0.03966 | -3.21% |
| G12 | 384 | 22.9278 | 23.4074 | 0.03988 | -2.05% |
| G12 | 512 | 23.1572 | 23.5013 | 0.03999 | -1.46% |

### 2.2 Observed orders
- **`K*`**: order ~ **1.0**.
- **`P*`**: observed order **0.73–0.78** (sub-second order).
- **`D*`**: observed order **0.82–0.84**.
- **Budget closure**: order ~ **1.0** — `F00`: -6.19 -> -3.20 -> -1.46;
  `G06`: -6.21 -> -3.20 -> -2.04 -> -1.46; `G12`: -6.22 -> -3.21 -> -2.05 -> -1.46.

### 2.3 Richardson extrapolation (indicative)
- `F00-R400`: `P*inf = 24.21` (residual -3.93%), `D*inf = 23.96` (residual -1.47%).
- `G06/G12-R400`: `P*inf ~ 24.2–24.3` but **residual ~ -5.25%** (extrapolation
  unreliable at Re=400).
- **Conclusion**: the absolute convergence of `P*` **is not reached** (gap ~ 4% at
  the extrapolation, extrapolation uncertainty >= 1%). The N=512 absolute values
  are used as the best direct estimate; **the conclusions rest on the trends and
  the differential quantities, not on the `P*inf` limits.**

---

## 3. Differential effect of the corrugation (contribution A)

`dQ = 100 * (Q_wavy - Q_flat) / Q_flat`, at the same grid, same profile, same
criteria.

### 3.1 Stability with resolution (Re=400)

| Case | N | d`P*` | d`D*` | d`K*` | vortices |
|---|---|---|---|---|---|
| G06 | 128 | -0.169% | -0.156% | -0.308% | 3 / 3 |
| G06 | 256 | -0.167% | -0.161% | -0.316% | 3 / 3 |
| G06 | 512 | **-0.166%** | -0.163% | -0.321% | 3 / 3 |
| G12 | 128 | -0.458% | -0.429% | -1.117% | 5 / 3 |
| G12 | 256 | -0.451% | -0.435% | -1.129% | 4 / 3 |
| G12 | 512 | **-0.447%** | -0.438% | -1.141% | 4 / 3 |

**Drift of d`P*` from N=128 to 512: 0.003 pp (G06) and 0.011 pp (G12)** — i.e.
~1.8% / 2.5% of the measured effect. The effect is therefore **stable over a x4
change in resolution**, despite an absolute budget residual of ~1.5%: an
**empirical cancellation of the common error** in the wavy-minus-flat difference.
Richardson extrapolation of d`P*` for G12: p ~ 0.81 -> d`P*inf` ~ -0.442%
(residual 1.18%), consistent with the N=512 value (-0.447%).

### 3.2 Dependence on Re and geometry (N=256)

| Case | Re | d`P*` | d`D*` | d`K*` |
|---|---|---|---|---|
| G06 | 100 | -0.013% | -0.013% | -0.265% |
| G12 | 100 | -0.024% | -0.022% | -0.765% |
| G06 | 400 | -0.167% | -0.161% | -0.316% |
| G12 | 400 | -0.451% | -0.435% | -1.129% |
| G06 | 1000 | -0.343% | -0.324% | -0.233% |
| G12 | 1000 | -1.049% | -1.001% | -1.062% |
| G13 | 400 | -0.730% | -0.715% | -2.085% |
| G13 | 1000 | -1.678% | -1.634% | -2.018% |

- The effect **grows with Re and with the amplitude `h`** (G13 at h/H=0.10 is the
  strongest).
- At fixed amplitude, **`Nw=1` (G13)** gives the most pronounced effect.

### 3.3 Cross-variant check (N=512)
Difference `consistent` vs `ma_fixed` (different `Ma`/`tau` parameters):

| Case | d`P*` | d`D*` | d`K*` |
|---|---|---|---|
| F00-R400 | 0.028% | 0.052% | -0.015% |
| G06-R400 | 0.028% | 0.053% | -0.013% |
| G12-R400 | 0.027% | 0.053% | -0.012% |

-> differences **<= 0.053%**, identical across the three cases: the measurement
does not depend on the numerical parameter set.

### 3.4 Physical robustness (topology)
- `G06`: **3 vortices at all resolutions** (stable).
- `G12`: **5 (N=128) -> 4 (N=256, 384, 512)** — the weakest secondary vortex falls
  at the detection threshold; to be refined (see §5).

---

## 4. Comparison of lid profiles (contribution B)

Re = 400, N = 256. Three profiles: **U** = uniform, **A** = regularised (same
*maximum* speed as U), **B** = `smoothstep_mean1` (same *mean* speed as U).

| Case | quantity | U | A | B | B/U-1 | A/U-1 |
|---|---|---|---|---|---|---|
| F00 | `K*` | 0.04016 | 0.04012 | 0.04186 | **+4.24%** | -0.10% |
| F00 | `P*` | 28.2330 | 22.5825 | 23.6453 | **-16.25%** | -20.01% |
| F00 | `D*` | 32.3100 | 23.3280 | 24.4286 | **-24.39%** | -27.80% |
| G06 | `P*` | 28.1952 | 22.5448 | 23.6050 | **-16.28%** | -20.04% |
| G13 | `P*` | 28.0681 | 22.4176 | 23.4685 | **-16.39%** | -20.13% |

### 4.1 Operational decomposition
- Effect of the **spatial distribution** (B vs U, equal mean): **-16.3% on `P*`**,
  **-24.4% on `D*`**.
- Effect of the **mean speed** (A vs B): **-3.75 pp on `P*`** — remarkably stable
  (-3.765 / -3.760 / -3.744 pp).

This decomposition (~80% distribution, ~20% mean) is **operational and depends on
the chosen normalisation**, not a universal law.

### 4.2 Geometrical cross-check
d`P*` under protocol B: G06 = **-0.170%**, G13 = **-0.748%** — compatible with the
campaign (N=512: -0.166; N=256: -0.730). **Renormalising the lid does not modify
the geometrical effect.**

### 4.3 Main caveat — budget closure
| lid | balance (P*-D*)/D* |
|---|---|
| uniform (V02-R400) | **-12.6%** |
| regularised (A) | -3.20% |
| `smoothstep_mean1` (B) | -3.21% |

The uniform profile has a velocity jump at **both lid corners**; the smoothed
profile avoids it. **A markedly larger budget residual must be interpreted as a
difference in numerical closure quality, not as a difference in physical
efficiency between profiles.** The `P*`/`D*` budgets (§4.1) must be read with this
explicit caveat.

### 4.4 Definition of the reference Reynolds number
Under protocol B the **mean** speed is 1 but the **maximum** speed of the
renormalised profile reaches 1.0211. **If `Re` is defined on the mean speed it is
preserved; if it is defined on the maximum speed it is not.** Any comparison of
Reynolds number between protocols must state this convention.

---

## 5. Limitations and conclusions

### 5.1 Established results
1. **Geometrical effect (A): weak but stable.** At Re=400 the corrugation reduces
   `P*` by **-0.17% (G06)** to **-0.45% (G12)** and `K*` by -0.32% to -1.14%.
   This effect **varies by <= 0.011 pp between N=128 and 512** and is
   **independent of the numerical variant** (<= 0.053%). The effect grows with Re
   and `h`.
2. **Lid-profile effect (B): dominant.** At equal mean speed, the spatial
   distribution of the velocity reduces `P*` by **-16.3%** and `D*` by **-24.4%**;
   renormalising the mean adds ~ **-3.75 pp**.

### 5.2 Defensible formulation
> At Re=400, the relative differences of power, dissipation, and kinetic energy
> between corrugated and flat geometries show a **weak drift between N=128, 256,
> and 512**. The observed stability, despite the evolution of the absolute values
> and the energy-budget residuals, indicates a **strong cancellation of common
> numerical errors** in the differential quantities. The results support the
> existence of a **weak but measurable energy effect** of the corrugation.

### 5.3 Limitations
- **Incomplete absolute convergence**: `P*` still moves by ~4% at the
  extrapolation; the `P*inf`/`D*inf` limits are indicative, not conclusive. The
  scope of the results is strictly **Re=400** (except as noted in §3.2).
- **Energy budget**: the residual is larger for the uniform lid (-12.6%) than for
  the smooth profiles (-3.2%); this must be documented explicitly.
- **G12 topology**: transition 5 -> 4 vortices with resolution (secondary vortex
  at the threshold) -> robustness to be refined on the **local fields**, not only
  on the integrals.
- **Differential stability alone does not prove the absence of any numerical
  error**; it shows that the effect hardly varies over the tested resolutions.

### 5.4 Complementary validations

**Non-equilibrium dissipation `D*_fneq`** (`Pi^neq` moments; F00-R400 / G06-R400
recomputed at N=256 with saved distributions; `D*_grad u` recalled on the same
fluid mask):

| Case (N=256) | `D*_grad u` | `D*_fneq` | ratio `fneq/grad u` | `D*_fneq` (dev) | balance `grad u` | balance `fneq` |
|---|---|---|---|---|---|---|
| F00-R400 | 23.3280 | 23.2923 | **0.9985** | 23.2655 (0.9973) | -3.20% | **-3.05%** |
| G06-R400 | 23.2904 | 23.2540 | **0.9984** | 23.2273 (0.9973) | -3.20% | **-3.05%** |

- `D*_fneq` reproduces `D*_grad u` to **0.15%** (0.27% for the deviatoric part):
  the measured dissipation **is not an artefact of the finite-difference scheme**.
- The budget residual changes from -3.20% to **-3.05%** with the non-equilibrium
  estimator: the `P ~ D` gap **persists**, so it cannot be attributed to the
  measurement of `D*`.
- The **differential** effect is robust to the choice of estimator:
  d`D*_fneq`(G06-F00) = -0.165% vs d`D*_grad u` = -0.161%.
- `D_star_fneq_error = None` for both cases; responses and `checks_ok=True`.

**Reserved leads (not launched)**:
- (i) `F00-R400` at `N=384`, `consistent` variant, for a `d_Q(384)` at four
  resolutions (not essential).
- (ii) `N=1024` campaign (`F00/G06/G12-R400`, `consistent` variant, script
  `src/convergence_study_n1024.py`) **reserved for a possible revision**. Launched
  then **stopped**: estimated cost ~7 days for a limited gain (gap to the `P*`
  limit reduced from ~3.9% to ~2.2–2.4%, absolute convergence remaining
  incomplete). Set aside in favour of the **stability of the differential effect**,
  already established at N=128/256/512.

### 5.5 Verdict
The data justify a manuscript structured around **two complementary phenomena**: a
**weak but stable energy effect of the corrugated floor**, and a **much stronger
influence of the velocity profile imposed on the lid**. The priority is no longer
an extended campaign, but the **consolidation of the data**, **transparency about
the energy budget**, and a writing that distinguishes established results from
still-conditional interpretations.

---

### Data provenance
| File | Content |
|---|---|
| `results/scalar_table.json` | 78 validated rows (`converged`, `checks_ok=True`) |
| `results/convergence/resolution.csv` | convergence N=128/256/384/512, 18/18 |
| `results/analysis/lid_profile.csv` | lid profiles (means, deficit) |
| `results/analysis/lid_sensitivity.csv` | regularised vs uniform decomposition |
| `results/analysis/wavy_effect.csv` | d_Q wavy-flat per resolution and Re |
| `results/analysis/convergence_order.csv` | observed orders + Richardson extrapolation |
| `results/analysis/protocol_b_compare.csv` | protocols A/B (F00, G06, G13) |
| `results_fneq/` | `D*_fneq` recomputations F00/G06-R400 N=256 (completed) |
