# Test de stabilite de I_Q en resolution (N=128 -> N=256)

Genere par `solver/run_res_test.py` le 2026-10-18 00:54:46.

Configuration : Re = 1000, Nw = 1, paire (n_a=2, n_b=4), alpha = 0.1, tol = 1e-09.

Question : l'erreur de resolution, qui deplace D* de +2.99 % entre 128 et 256 pour le cas (n=2, alpha=0), s'annule-t-elle dans le contraste I_Q ?

## 1. Cas calcules

| N | cas | status | iters | residual | D* | P* | K* | bilan |
|---|---|---|---|---|---|---|---|---|
| 128 | n=2 alpha=0.00 | converged | 232000 | 7.99e-10 | 30.345879 | 28.278210 | 0.041786 | -6.8137% |
| 128 | n=2 alpha=0.10 | converged | 226000 | 8.62e-10 | 30.003394 | 27.953686 | 0.040992 | -6.8316% |
| 128 | n=4 alpha=0.00 | converged | 230000 | 8.39e-10 | 30.521151 | 28.452428 | 0.041982 | -6.7780% |
| 128 | n=4 alpha=0.10 | converged | 230000 | 8.58e-10 | 30.036076 | 27.977288 | 0.041161 | -6.8544% |
| 256 | n=2 alpha=0.00 | converged | 768000 | 9.52e-10 | 31.253399 | 30.168938 | 0.042944 | -3.4699% |
| 256 | n=2 alpha=0.10 | converged | 758000 | 9.56e-10 | 30.897478 | 29.827579 | 0.042100 | -3.4627% |
| 256 | n=4 alpha=0.00 | converged | 760000 | 9.55e-10 | 31.437075 | 30.351523 | 0.043166 | -3.4531% |
| 256 | n=4 alpha=0.10 | converged | 772000 | 9.54e-10 | 30.934522 | 29.851469 | 0.042284 | -3.5011% |

## 2. Sensibilite absolue a la resolution

| cas | dD*/D* | dP*/P* | dK*/K* |
|---|---|---|---|
| n=2 alpha=0.00 | 2.991 % | 6.686 % | 2.773 % |
| n=2 alpha=0.10 | 2.980 % | 6.704 % | 2.704 % |
| n=4 alpha=0.00 | 3.001 % | 6.675 % | 2.819 % |
| n=4 alpha=0.10 | 2.991 % | 6.699 % | 2.728 % |

## 3. Interactions I_Q aux deux resolutions

| grandeur | I (N=128) | I (N=256) | ecart absolu | ecart relatif | verdict |
|---|---|---|---|---|---|
| `D_star` | -0.142590 | -0.146631 | -4.041e-03 | 2.834 % | EXCELLENT (<10%) |
| `P_star` | -0.150617 | -0.158696 | -8.079e-03 | 5.364 % | EXCELLENT (<10%) |
| `K_star` | -0.000028 | -0.000038 | -1.055e-05 | 37.996 % | PROBLEME (>30%) |

## 4. Composantes intermediaires

| grandeur | delta(n=2) @N=128 | delta(n=4) @N=128 | delta(n=2) @N=256 | delta(n=4) @N=256 |
|---|---|---|---|---|
| `D_star` | -0.342485 | -0.485075 | -0.355922 | -0.502553 |
| `P_star` | -0.324523 | -0.475140 | -0.341359 | -0.500055 |
| `K_star` | -0.000794 | -0.000821 | -0.000844 | -0.000882 |

## 5. Rapport signal / resolution

D* bouge de 2.991 % entre les deux resolutions (cas n=2, alpha=0). Le signal d'interaction vaut -0.146631 a N=256.

| grandeur | |I(N=256)| / |dQ*_res| | lecture |
|---|---|---|
| `D_star` | 0.16 | zone limite : annulation partielle attendue |
| `P_star` | 0.08 | zone limite : annulation partielle attendue |
| `K_star` | 0.03 | signal petit devant le bruit de resolution |

## 6. Verdict global

| critere | mesure | seuil | verdict |
|---|---|---|---|
| `I_D` : ecart relatif N=128 -> N=256 | 2.834 % | 10 % / 30 % | EXCELLENT (<10%) |
| `I_P` : ecart relatif | 5.364 % | 10 % / 30 % | EXCELLENT (<10%) |
| `I_K` : ecart relatif | 37.996 % | 10 % / 30 % | PROBLEME (>30%) |

**Decision : STABLE.** I_D est stable a 2.834 % entre N=128 et N=256 : le signal d'interaction est robuste face a la resolution. Les campagnes N=256 puis N=512 peuvent etre lancees telles qu'elles sont planifiees.

## 7. Paragraphe de redaction (papier 3)

> The interaction contrast must be shown to be free of resolution
> contamination. We therefore recomputed the full \((2,4)\) pair at
> \(Re = 1000\) for \(N = 128\) and \(N = 256\). The dissipation of the reference
> configuration \((n=2,\alpha=0)\) shifts by 2.99 \% between the two grids,
> yet the contrast changes only from \(I_D = -0.14259\) at \(N=128\) to
> \(I_D = -0.14663\) at \(N=256\), i.e. 2.83 \% of relative change (excellent (<10%)).
> The matched-resolution cancellation is therefore effective to within that
> figure, which is used as the systematic uncertainty of the reported
> interaction.
