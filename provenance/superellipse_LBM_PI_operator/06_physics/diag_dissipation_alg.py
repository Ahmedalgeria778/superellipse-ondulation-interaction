# -*- coding: utf-8 -*-
"""
06_physics/diag_dissipation_alg.py
==================================

TEST DECISIF DE LA CONVENTION ENERGETIQUE (Reserve C).

Probleme : le tableau V2 rapportait P/D = 1.80 (Re=100) et 1.99 (Re=1000) sur la
cavite de Ghia N=256, et la production rapportait P/Phi = 4.27 ... 1.47.  Deux
interpretations s'affrontent :

  (H1) la convention Φ = 0.5 U (N-1) D_code est correcte, et l'ecart est un
       defaut de resolution (masque gradient_valid_mask qui exclut la 1re couche
       fluide contre la paroi) ;
  (H2) la convention Φ = 0.5 U (N-1) D_code est FAUSSE d'un facteur, et le
       "vrai" bilan est P ≈ 2 D_code (ce que l'ecart apparent ~2 suggère).

On ne peut pas trancher en comparant deux formules : il faut une MESURE, sans
convention.  LBM dissipe l'energie par l'operateur de collision, et cette
quantite est definie SANS AMBIGUITE :

    A_alg = somme_x somme_q ( f_post[q,x] - f[q,x] ) (e_q + u(x))^2 / 2

C'est exactement l'energie retiree du systeme par la collision, par pas de
temps, en unites lattice.  Aucune hypothese sur S:S, sur le facteur 2, sur le
masque, sur le signe.  En regime stationnaire l'energie ne peut ni s'accumuler
ni disparaitre, donc A_alg DOIT egaler la puissance injectee par les parois.

Ce script compare donc, sur des etats deja sauvegardes (checkpoints de
production, aucune re-simulation) :

    A_alg                      dissipation reellement dissipee par l'algo
    P_MEM                      puissance de paroi (MEM Ladd), code de production
    Phi_code = 0.5 U (N-1) D   convention actuelle de analyse_etape4_puissance.py
    2 D_code                   convention "D_h, Phi = 2 D_h"

Verdict : la convention correcte est celle qui reproduce A_alg.

Ecrit dans 06_physics/results_puissance/diag_dissipation_alg.{md,csv}.
Ne modifie AUCUNE production, AUCUN resultat existant.
"""
from __future__ import annotations

import csv
import sys
import time
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
_ROOT = _HERE.parent
for _p in (_ROOT / "02_LBM", _ROOT / "01_geometry", _HERE):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

from superellipse import Superellipse                     # noqa: E402
from solver_curved import CurvedLBM                       # noqa: E402
from mrt_collision import EX, EY, collide, macroscopic, s_diag  # noqa: E402
from geometry_links import build_link_arrays               # noqa: E402
from etude_convergence_paroi import gradient_valid_mask   # noqa: E402
from valider_bilan_puissance import (U, TOL, grad_D,      # noqa: E402
                                     mem_stream_force)

OUT = _HERE / "results_puissance"


def find_ckpt(N):
    """Cherche un checkpoint PRODUCTION de n=8, Re=1000 a la resolution N."""
    pat = f"se_n8p0_Re1000_N{N}"
    best = None
    for sub in ("data_conv256", "data_convHR", "data_corner", "data_conv"):
        d = _ROOT / "05_dataset" / sub
        if not d.is_dir():
            continue
        for p in sorted(d.glob(pat + "*")):
            if p.suffix != ".npz":
                continue
            # preference au checkpoint (contient f), sinon le cas final
            rank = 0 if "ckpt" in p.name else 1
            if best is None or rank < best[0]:
                best = (rank, p)
    return None if best is None else best[1]


def alg_dissipation(f, f_post, rho, ux, uy, fluid):
    """Deux mesures, sans convention, de la dissipation par la collision.

    A_ghost : somme_x somme_q (f_post - f) (e_q+u)^2 / 2
              (fonctionnelle d'energie que preserve le streaming LBM)

    A_phys  : - somme_x 1/2 rho ( |u'|^2 - |u|^2 )
              (energie CINETIQUE PHYSIQUE, seule comparable a 2 rho nu int S:S)
    """
    tot = 0.0
    for q in range(9):
        uq2 = (EX[q] + ux) ** 2 + (EY[q] + uy) ** 2
        tot += float(np.sum((f_post[q] - f[q]) * uq2)) * 0.5
    rp, up, vp = macroscopic(f_post)
    dE = 0.0
    for (a, b) in ((ux, uy), (up, vp)):
        dE += float(np.sum(rho[fluid] * (a[fluid] ** 2 + b[fluid] ** 2)))
    e_before = 0.5 * float(np.sum(rho[fluid] * (ux[fluid] ** 2 + uy[fluid] ** 2)))
    e_after = 0.5 * float(np.sum(rho[fluid] * (up[fluid] ** 2 + vp[fluid] ** 2)))
    return tot, e_before - e_after, e_before


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    md = [
        "# Test decisif de la convention energetique (Reserve C)",
        "",
        "`A_alg` = energie retiree par la collision, mesuree sur l'etat",
        "reellement calcule : `somme_x somme_q (f_post - f) (e_q+u)^2 / 2`.",
        "Aucune convention supposee. En regime stationnaire `A_alg` doit",
        "egaler la puissance de paroi `P`.", "",
        "| cas | N | it | A_phys | A_ghost | P_MEM | Phi_code | "
        "Phi_sans_masque | 2 D_code | P/A | Ag/A | Phi/A | Phi_nom/A | 2D/A |",
        "|---|---|---|---|---|---|---|---|---|---|---|---|---|",
    ]
    for N in (128, 256, 512):
        p = find_ckpt(N)
        if p is None:
            print(f"[diag] N={N} : aucun checkpoint trouve", flush=True)
            continue
        t0 = time.time()
        with np.load(p) as z:
            f = np.asarray(z["f"], float)
            it = int(np.asarray(z["it"])) if "it" in z.files else -1
        geom = Superellipse(n=8.0)
        Re = 1000.0
        la = build_link_arrays(geom, N, U_lid=U, lid_mode="upper")
        S = s_diag(Re, N, L_char=N - 1, U_ref=U)
        fluid = la.fluid.astype(bool)
        X = la.grid.X
        dx = float(abs(X[0, 1] - X[0, 0]))

        rho, ux, uy = macroscopic(f)
        f_post = collide(f, S, fluid)
        A_ghost, A_phys, E_phys = alg_dissipation(f, f_post, rho, ux, uy, fluid)

        _, _, P = mem_stream_force(f_post, la, rho, None)
        D_code, valid = grad_D(ux, uy, fluid, dx, Re)
        Phi_code = 0.5 * U * (N - 1) * D_code
        D2 = 2.0 * D_code
        # D sans masque valide (toutes les cellules fluides, grille allegee)
        gx = np.full(ux.shape, np.nan); gy = np.full(ux.shape, np.nan)
        vx = np.full(ux.shape, np.nan); vy = np.full(ux.shape, np.nan)
        gx[1:-1, 1:-1] = (ux[1:-1, 2:] - ux[1:-1, :-2]) / (2 * dx)
        vx[1:-1, 1:-1] = (uy[1:-1, 2:] - uy[1:-1, :-2]) / (2 * dx)
        gy[1:-1, 1:-1] = (ux[2:, 1:-1] - ux[:-2, 1:-1]) / (2 * dx)
        vy[1:-1, 1:-1] = (uy[2:, 1:-1] - uy[:-2, 1:-1]) / (2 * dx)
        inb = np.zeros_like(fluid, bool)
        inb[1:-1, 1:-1] = True
        m_all = inb & fluid
        phi_all = np.where(m_all, (2 * gx ** 2 + 2 * vy ** 2 + (gy + vx) ** 2) / Re, 0.0)
        D_all = 2.0 * float(np.sum(phi_all)) * dx * dx
        Phi_all = 0.5 * U * (N - 1) * D_all

        rows.append({"N": N, "it": it, "source": str(p.relative_to(_ROOT)),
                     "A_ghost": A_ghost, "A_phys": A_phys, "E_phys": E_phys,
                     "P_MEM": P, "Phi_code": Phi_code, "D_code": D_code,
                     "Phi_sans_masque": Phi_all,
                     "deux_D_code": D2,
                     "P_sur_A_phys": P / A_phys,
                     "A_ghost_sur_A_phys": A_ghost / A_phys,
                     "Phi_sur_A_phys": Phi_code / A_phys,
                     "Phi_sans_masque_sur_A_phys": Phi_all / A_phys,
                     "deuxD_sur_A_phys": D2 / A_phys,
                     "n_fluid": int(fluid.sum()),
                     "n_masque": int(valid.sum())})
        md.append(
            f"| n=8 Re=1000 | {N} | {it} | {A_phys:+.6e} | {A_ghost:+.6e} | "
            f"{P:+.6e} | {Phi_code:+.6e} | {Phi_all:+.6e} | {D2:+.6e} | "
            f"{P / A_phys:.4f} | {A_ghost / A_phys:.4f} | "
            f"{Phi_code / A_phys:.4f} | {Phi_all / A_phys:.4f} | "
            f"{D2 / A_phys:.4f} |")
        print(f"[diag] N={N} it={it}\n"
              f"        A_phys={A_phys:+.6e}  A_ghost={A_ghost:+.6e}  "
              f"E_phys={E_phys:+.6e}\n"
              f"        P_MEM={P:+.6e}  Phi_code={Phi_code:+.6e}  "
              f"Phi_sans_masque={Phi_all:+.6e}  2D={D2:+.6e}\n"
              f"        P/A_phys={P / A_phys:.4f}  A_ghost/A_phys="
              f"{A_ghost / A_phys:.4f}  Phi/A_phys={Phi_code / A_phys:.4f}  "
              f"Phi_nomasq/A_phys={Phi_all / A_phys:.4f}  "
              f"2D/A_phys={D2 / A_phys:.4f}  ({time.time() - t0:.0f}s)", flush=True)

    md += [
        "",
        "## Lecture",
        "",
        "- `A_phys` = -1/2 somme rho (|u'|^2 - |u|^2) : dissipation mesuree sur",
        "  l'energie cinematique physique. C'est la seule reference comparable a",
        "  la dissipation mecanique 2 rho nu int S:S.",
        "- `A_ghost` = somme (f_post - f) (e_q+u)^2/2 : version LBM complete",
        "  (fonctionnelle que le streaming preserve).",
        "- `Phi_sans_masque` = 0.5 U (N-1) D calcule sur TOUTES les cellules",
        "  fluides (sans `gradient_valid_mask`), pour isoler l'effet du masque.",
        "- Une colonne quotient valant 1 identifie la convention correcte ; on",
        "  retient celle qui vaut 1 aux DEUX resolutions.",
    ]
    (OUT / "diag_dissipation_alg.md").write_text("\n".join(md), encoding="utf-8")
    if rows:
        with (OUT / "diag_dissipation_alg.csv").open("w", newline="",
                                                      encoding="utf-8") as fp:
            w = csv.DictWriter(fp, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)
    print(f"[diag] ecrit -> {OUT / 'diag_dissipation_alg.md'}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())