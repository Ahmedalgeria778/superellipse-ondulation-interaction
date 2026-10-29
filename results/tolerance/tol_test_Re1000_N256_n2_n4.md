# Test de sensibilite de I_Q a la tolérance de convergence

Genere par `solver/run_tol_test.py` le 2026-10-18 00:51:20.

Configuration : Re = 1000, N = 256, Nw = 1, paire (n_a=2, n_b=4), alpha = [0.0, 0.1].

Tolérances comparees : **1e-09** (referentiel) et **1e-08**.

## 1. Cas calcules

| tol | cas | status | iters | residual | D* | P* | K* | bilan | t (s) |
|---|---|---|---|---|---|---|---|---|---|
| 1e-09 | n=2 alpha=0.00 | converged | 768000 | 9.52e-10 | 31.253399 | 30.168938 | 0.042944 | -3.4699% | 2614 |
| 1e-09 | n=2 alpha=0.10 | converged | 758000 | 9.56e-10 | 30.897478 | 29.827579 | 0.042100 | -3.4627% | 9914 |
| 1e-09 | n=4 alpha=0.00 | converged | 760000 | 9.55e-10 | 31.437075 | 30.351523 | 0.043166 | -3.4531% | 10088 |
| 1e-09 | n=4 alpha=0.10 | converged | 772000 | 9.54e-10 | 30.934522 | 29.851469 | 0.042284 | -3.5011% | 9976 |
| 1e-08 | n=2 alpha=0.00 | converged | 608000 | 9.52e-09 | 31.253361 | 30.168949 | 0.042944 | -3.4697% | 8110 |
| 1e-08 | n=2 alpha=0.10 | converged | 608000 | 9.64e-09 | 30.897445 | 29.827594 | 0.042100 | -3.4626% | 8206 |
| 1e-08 | n=4 alpha=0.00 | converged | 594000 | 9.56e-09 | 31.437034 | 30.351534 | 0.043165 | -3.4529% | 8309 |
| 1e-08 | n=4 alpha=0.10 | converged | 620000 | 9.63e-09 | 30.934489 | 29.851484 | 0.042283 | -3.5010% | 8575 |

## 2. Interactions I_Q (difference de differences)

| grandeur | I (tol=1e-09) | I (tol=1e-08) | ecart absolu | ecart relatif | seuil 0.1 % | seuil 1 % |
|---|---|---|---|---|---|---|
| `D_star` | -0.146631 | -0.146628 | +2.644e-06 | 0.0018 % | OK | OK |
| `P_star` | -0.158696 | -0.158694 | +1.655e-06 | 0.0010 % | OK | OK |
| `K_star` | -0.000038 | -0.000038 | +3.251e-08 | 0.0848 % | OK | OK |

Composantes intermediaires (utilisees pour lire ou vient l'ecart) :

| grandeur | delta(n=2) @1e-09 | delta(n=4) @1e-09 | delta(n=2) @1e-08 | delta(n=4) @1e-08 |
|---|---|---|---|---|
| `D_star` | -0.355922 | -0.502553 | -0.355916 | -0.502544 |
| `P_star` | -0.341359 | -0.500055 | -0.341356 | -0.500050 |
| `K_star` | -0.000844 | -0.000882 | -0.000844 | -0.000882 |

## 3. Bilan énergétique (P - D)/D

| cas | bilan @1e-09 | bilan @1e-08 | variation relative | seuil 5 % |
|---|---|---|---|---|
| n=2 alpha=0.00 | -3.469900% | -3.469743% | 0.0045 % | OK |
| n=2 alpha=0.10 | -3.462739% | -3.462589% | 0.0043 % | OK |
| n=4 alpha=0.00 | -3.453093% | -3.452933% | 0.0046 % | OK |
| n=4 alpha=0.10 | -3.501116% | -3.500965% | 0.0043 % | OK |

## 4. Verdicts

| critère | mesure | seuil | verdict |
|---|---|---|---|
| `I_D` : ecart relatif entre 1e-09 et 1e-08 | 0.0018 % | 0.1 % / 1 % | PASS<0.1%% |
| `I_P` : ecart relatif | 0.0010 % | 0.1 % / 1 % | PASS<0.1%% |
| `I_K` : ecart relatif | 0.0848 % | 0.1 % / 1 % | PASS<0.1%% |
| bilan : ecart relatif maximal | 0.0046 % | 5 % | PASS |
| D* : ecart relatif maximal (par cas) | 0.0001 % | 0.1 % | PASS |
| P* : ecart relatif maximal (par cas) | 0.0001 % | 0.1 % | PASS |

**Decision : ADOPTER_1e-8** -- la tolérance 1e-08 est validee pour les campagnes N=256 et N=512.

## 5. Etudes de tolérance deja disponibles (appuis)

| etude | localisation | resultats |
|---|---|---|
| Projet B -- audit_stabilite_D.py | `superellipse_LBM_PI_operator/06_physics/results_audit_stabilite_D/stabilite_D.md` | D : mediane 1.005e-04 %, max 2.392e-03 % (47 cas, ckpt vs final) |
| Projet B -- canonique_ghia_re1000.py | `superellipse_LBM_PI_operator/06_physics/results_canonique_ghia/canonique_ghia_Re1000.md` | fermeture : 17.535 % (1e-8) -> 17.527 % (1e-9) = 0.046 % relatif |

Aucune de ces etudes ne porte sur un contraste `I_Q` : c'est l'objet de ce test.

## 6. Paragraphe de rédaction (papier 3)

> The tolerance of the steady-state criterion is set to \(\epsilon = 10^{-9}\) at \(N = 128\).
> The sensitivity of the integrated quantities to \(\epsilon\) was assessed in two independent studies of the present solver family: the viscous dissipation \(D\) varies by at most \(2.4\times 10^{-3}\,\%\) over 47 configurations
> when the run is stopped 20\,000--24\,000 iterations earlier, and the closure residual \((P-D)/D\) changes by \(0.046\,\%\) between \(\epsilon = 10^{-8}\) and \(\epsilon = 10^{-9}\)
> at \(Re=1000\), \(N=256\). For the interaction contrast \(I_D\), which is a difference of differences and is therefore not covered by those two tests, we performed a
> dedicated experiment at \(Re=1000\), \(N=256\), pair \((2,4)\): \(I_D = -0.146631\) at \(\epsilon = 10^{-9}\) and \(I_D = -0.146628\) at \(\epsilon = 10^{-8}\), i.e. a relative change of 0.0018 %,
> well below the measured signal itself and below our 1 % acceptance threshold. All campaigns at \(N \ge 256\) were therefore run with \(\epsilon = 10^{-8}\).
