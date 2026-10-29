# -*- coding: utf-8 -*-
"""
solver/audit_balance.py
=======================
Priorite 2 du relecteur : d'ou vient le defaut d'equilibre energetique
(P* - D*)/D* mesure a 3-7 % ?

Quatre causes candidates, toutes testees sur les CHAMPS DEJA CALCULES
(aucune simulation) :

  (i)  E_p    gradient de vitesse impose sur la lid en ordre 1
              (paroi a 0.5*dx du premier noeud) -> dans P* et dans S12 ;
  (ii) E_sign inversion de signe des corrections unilaterales quand la
              paroi est du cote NEGATIF de l'axe (paroi inferieure :
              (0-U)/dist au lieu de +U/dist ; paroi gauche : idem) ;
              cf. diagnostics.py:57 et :73 ;
  (iii) E_fd  contamination du gradient par le saut solide/fluide :
              np.gradient derive a travers le front (noeuds solides a 0),
              et seules dudy/dvdy (bas) et dudx/dvdx (cotes) sont corriges ;
  (iv)  E_mur couche de paroi : toutes les corrections unilaterales sont
              de premier ordre.

Modes de gradient implementes :
    ref      reproduction exacte de diagnostics._gradients
    signfix  ref + signes corriges en bas et a gauche
    lidonly  ref + ordre 2 sur la LID SEULEMENT (P et D coherents)
    ord2     formule a 3 points (ordre 2) sur les QUATRE parois
    clean    differences centrees uniquement entre noeuds fluides
             (aucun derive a travers le front fluide/solide)

Le bilan est evalue pour chaque croisement (P, D) :
    eps_ref   = P_ref    - D_ref        etat publie
    eps_p2    = P_ord2   - D_ref        seule la lid passe en ordre 2 (dans P)
    eps_lid   = P_ord2   - D_lidonly    lid en ordre 2 de facon coherente
    eps_sign  = P_ref    - D_signfix    seul le signe est corrige
    eps_full  = P_ord2   - D_ord2       passage complet en ordre 2
    eps_clean = P_ref    - D_clean      seul le saut solide/fluide est retire

Usage :
  python -m solver.audit_balance
  python -m solver.audit_balance --Re 400,1000 --N 128,256 --n-list 2,8
"""

from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
_PROV_A = os.path.join(_ROOT, "provenance", "Lid_driven_ondulation")
for _p in (_HERE, _PROV_A):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from solver.combined_geometry import build_combined_geometry  # noqa: E402
from solver.reanalyze_pilot import lid_profile_regularized    # noqa: E402

FIELDS_DIR = os.path.join(_ROOT, "results", "fields")
AUDIT_DIR = os.path.join(_ROOT, "results", "audit")


# --------------------------------------------------------------------------- #
# derivees unilaterales
# --------------------------------------------------------------------------- #
def _wall3(u0, u1, u2, h, d):
    """du/dt en t=0 (paroi), points echantillones a t = 0, h, h+d.

    Interpolation quadratique exacte (ordre 2) pour espacement inegal.
    Le signe retourne est celui de la derivee par rapport a t, t croissant
    depuis la paroi vers l'interieur ; l'appelant applique le signe propre
    a l'axe consider.
    """
    A = u1 - u0
    B = u2 - u0
    b = (B - A * (h + d) / h) / (d * (h + d))
    return (A - b * h * h) / h


def _grad_clean(F, d, axis, fluid):
    """Gradient sans saut solide/fluide (segments contigus de noeuds fluides).

    Difference centree ordre 2 a l'interieur d'un segment, unilaterale
    ordre 1 aux extremites. Contrairement a np.gradient, aucun derive n'est
    pris a travers la frontiere fluide/solide.
    """
    n = F.shape[axis]
    Fm = np.moveaxis(np.asarray(F, dtype=float), axis, -1)
    fm = np.moveaxis(np.asarray(fluid, dtype=bool), axis, -1)
    Gm = np.zeros_like(Fm)

    for idx in np.ndindex(fm.shape[:-1]):
        lf, lF = fm[idx], Fm[idx]
        g = Gm[idx]
        i = 0
        while i < n:
            if not lf[i]:
                i += 1
                continue
            j = i
            while j + 1 < n and lf[j + 1]:
                j += 1
            if j > i:
                g[i] = (lF[i + 1] - lF[i]) / d
                g[j] = (lF[j] - lF[j - 1]) / d
                if j - i >= 2:
                    g[i + 1:j] = (lF[i + 2:j + 1] - lF[i:j - 1]) / (2.0 * d)
            i = j + 1
    return np.moveaxis(Gm, -1, axis)


def gradients_variant(U, V, geom, U_lid_norm, mode="ref"):
    """Gradient de vitesse selon le mode ; meme convention que
    diagnostics._gradients (les memes entrees sont ecrasees)."""
    N = geom["N"]
    d = geom["dx"]
    fluid = geom["fluid"].astype(bool)
    half = 0.5 * d
    three = 3.0 * d

    if mode == "clean":
        dudx = _grad_clean(U, d, 1, fluid)
        dudy = _grad_clean(U, d, 0, fluid)
        dvdx = _grad_clean(V, d, 1, fluid)
        dvdy = _grad_clean(V, d, 0, fluid)
    else:
        dudx = np.gradient(U, d, axis=1)
        dudy = np.gradient(U, d, axis=0)
        dvdx = np.gradient(V, d, axis=1)
        dvdy = np.gradient(V, d, axis=0)

    # --- paroi superieure (lid) : y croit vers le haut, paroi au-dessus ---
    for j in range(N):
        if not fluid[N - 1, j]:
            continue
        if mode in ("ord2", "lidonly") and fluid[N - 2, j]:
            dudy[N - 1, j] = -_wall3(U_lid_norm[j], U[N - 1, j],
                                     U[N - 2, j], half, d)
            dvdy[N - 1, j] = -_wall3(0.0, V[N - 1, j], V[N - 2, j], half, d)
        else:
            dudy[N - 1, j] = (U_lid_norm[j] - U[N - 1, j]) / half
            dvdy[N - 1, j] = (0.0 - V[N - 1, j]) / half

    # --- parois laterales (x croit vers la droite) ---
    for i in range(N):
        if fluid[i, 0]:
            if mode == "ord2" and fluid[i, 1]:
                dudx[i, 0] = _wall3(0.0, U[i, 0], U[i, 1], half, d)
                dvdx[i, 0] = _wall3(0.0, V[i, 0], V[i, 1], half, d)
            elif mode == "signfix":
                dudx[i, 0] = U[i, 0] / half
                dvdx[i, 0] = V[i, 0] / half
            else:
                dudx[i, 0] = (0.0 - U[i, 0]) / half
                dvdx[i, 0] = (0.0 - V[i, 0]) / half
        if fluid[i, N - 1]:
            if mode == "ord2" and fluid[i, N - 2]:
                dudx[i, N - 1] = -_wall3(0.0, U[i, N - 1], U[i, N - 2],
                                         half, d)
                dvdx[i, N - 1] = -_wall3(0.0, V[i, N - 1], V[i, N - 2],
                                         half, d)
            else:
                dudx[i, N - 1] = (0.0 - U[i, N - 1]) / half
                dvdx[i, N - 1] = (0.0 - V[i, N - 1]) / half

    # --- paroi inferieure ondulee : y croit vers le haut, paroi en dessous ---
    yb = geom["yb"]
    for j in range(N):
        col = np.nonzero(fluid[:, j])[0]
        if col.size == 0:
            continue
        i0 = int(col[0])
        h = (i0 + 0.5) * d - float(yb[0, j])
        if h <= 1e-12:
            continue
        if mode == "ord2" and i0 + 1 < N and fluid[i0 + 1, j]:
            dudy[i0, j] = _wall3(0.0, U[i0, j], U[i0 + 1, j], h, d)
            dvdy[i0, j] = _wall3(0.0, V[i0, j], V[i0 + 1, j], h, d)
        elif mode == "signfix":
            dudy[i0, j] = U[i0, j] / h
            dvdy[i0, j] = V[i0, j] / h
        else:
            dudy[i0, j] = (0.0 - U[i0, j]) / h
            dvdy[i0, j] = (0.0 - V[i0, j]) / h

    return dudx, dudy, dvdx, dvdy


# --------------------------------------------------------------------------- #
def _Dstar(dudx, dudy, dvdx, dvdy, fluid, d):
    S11, S22, S12 = dudx, dvdy, 0.5 * (dudy + dvdx)
    return float(np.sum(2.0 * (S11 ** 2 + S22 ** 2 + 2.0 * S12 ** 2)[fluid])
                 * d * d)


def _Pstar(dudy, U_lid_norm, fluid, d, N):
    P = 0.0
    for j in range(N):
        if fluid[N - 1, j]:
            P += U_lid_norm[j] * dudy[N - 1, j]
    return float(P * d)


def _regions(dudx, dudy, dvdx, dvdy, fluid, d, N):
    """Decomposition disjointe de D* : laterales / lid / plafond / interieur."""
    S11, S22, S12 = dudx, dvdy, 0.5 * (dudy + dvdx)
    cell = 2.0 * (S11 ** 2 + S22 ** 2 + 2.0 * S12 ** 2) * d * d

    j = np.arange(N)[None, :]
    i = np.arange(N)[:, None]
    has = fluid.any(axis=0)
    first = np.where(has, np.argmax(fluid, axis=0), N)

    m_side = fluid & ((j <= 1) | (j >= N - 2))
    m_lid = fluid & ~m_side & (i >= N - 2)
    m_bot = fluid & ~m_side & ~m_lid & (i <= first[None, :] + 1)
    m_int = fluid & ~m_side & ~m_lid & ~m_bot
    out = {"side": float(cell[m_side].sum()), "lid": float(cell[m_lid].sum()),
           "bottom": float(cell[m_bot].sum()), "interior": float(cell[m_int].sum())}
    out["wall_band"] = out["side"] + out["lid"] + out["bottom"]
    return out


# --------------------------------------------------------------------------- #
def load_cases(Re_list, N_list, n_list, alpha_list):
    out = {}
    for Re in Re_list:
        for N in N_list:
            for n in n_list:
                for al in alpha_list:
                    base = "se_n%g_w%.2f_Nw1_Re%d_N%d" % (n, al, Re, N)
                    for suf in ("_tol1e-08", ""):
                        path = os.path.join(FIELDS_DIR, base + suf + ".npz")
                        if os.path.exists(path):
                            out[(Re, N, n, al)] = path
                            break
    return out


def analyse(path, eps_lid=0.02):
    z = np.load(path)
    N = int(z["N"])
    n = float(z["n"])
    alpha = float(z["alpha"])
    Re = float(z["Re"])
    ux, uy = z["ux"], z["uy"]

    geom = build_combined_geometry(N, n, alpha, 1, 0.5, 0.5, 1.0)
    dx = geom["dx"]
    fluid = geom["fluid"].astype(bool)
    U = np.where(fluid, ux, 0.0)
    V = np.where(fluid, uy, 0.0)
    U_lid_norm = lid_profile_regularized(geom["xc"], eps_lid)

    G = {m: gradients_variant(U, V, geom, U_lid_norm, m)
         for m in ("ref", "signfix", "lidonly", "ord2", "clean")}
    D = {m: _Dstar(g[0], g[1], g[2], g[3], fluid, dx) for m, g in G.items()}
    P = {m: _Pstar(g[1], U_lid_norm, fluid, dx, N) for m, g in G.items()}
    # P* ne depend que de dudy au niveau de la lid : il est identique pour
    # ref/signfix/clean, et egal a la variante ordre 2 pour lidonly/ord2.
    P_ref = P["ref"]
    P_ord2 = P["ord2"]

    reg = _regions(*G["ref"], fluid, dx, N)

    e = {
        "ref": (P_ref - D["ref"]) / D["ref"] * 100.0,
        "p2": (P_ord2 - D["ref"]) / D["ref"] * 100.0,
        "lid": (P_ord2 - D["lidonly"]) / D["ref"] * 100.0,
        "sign": (P_ref - D["signfix"]) / D["ref"] * 100.0,
        "full": (P_ord2 - D["ord2"]) / D["ref"] * 100.0,
        "clean": (P_ref - D["clean"]) / D["ref"] * 100.0,
    }
    return {
        "Re": Re, "N": N, "n": n, "alpha": alpha,
        "D": D, "P_ref": P_ref, "P_ord2": P_ord2, "eps": e,
        "regions": reg,
        "wall_share_pct": reg["wall_band"] / D["ref"] * 100.0,
        "dD_sign_pct": (D["signfix"] - D["ref"]) / D["ref"] * 100.0,
        "dD_clean_pct": (D["clean"] - D["ref"]) / D["ref"] * 100.0,
        "dD_lidonly_pct": (D["lidonly"] - D["ref"]) / D["ref"] * 100.0,
        "dD_ord2_pct": (D["ord2"] - D["ref"]) / D["ref"] * 100.0,
        "dP_ord2_pct": (P_ord2 - P_ref) / D["ref"] * 100.0,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--Re", type=str, default="400,1000")
    ap.add_argument("--N", type=str, default="128,256")
    ap.add_argument("--n-list", type=str, default="2,8")
    ap.add_argument("--alpha-list", type=str, default="0,0.10")
    args = ap.parse_args()

    Re_list = [int(x) for x in args.Re.split(",")]
    N_list = [int(x) for x in args.N.split(",")]
    n_list = [float(x) for x in args.n_list.split(",")]
    alpha_list = [float(x) for x in args.alpha_list.split(",")]

    cases = load_cases(Re_list, N_list, n_list, alpha_list)
    if not cases:
        print("aucun champ trouve")
        return 1

    print("=" * 106)
    print("AUDIT DU BILAN ENERGETIQUE   (P* - D*) / D*")
    print("=" * 106)
    print("  ref    = etat publie (P et D en ordre 1 sur les parois)")
    print("  p2     = P en ordre 2 sur la lid, D publie")
    print("  lid    = P et D en ordre 2 sur la lid, parois restantes publiees")
    print("  sign   = signes corriges en bas/gauche, P publie")
    print("  full   = P et D en ordre 2 sur les quatre parois")
    print("  clean  = D sans derive a travers le front fluide/solide")
    print()
    print("  Re    N     n    alpha |   ref      p2      lid      sign     full"
          "     clean   | correction dominante")
    print("  " + "-" * 106)

    report = {"cases": [], "by_source": {}}
    rows = []
    for key in sorted(cases, key=lambda k: (k[0], k[1], k[2], k[3])):
        r = analyse(cases[key])
        rows.append(r)
        e = r["eps"]
        base = abs(e["ref"])
        gains = {}
        for lab, v in (("ordre-2 lid", "p2"), ("lid coherente", "lid"),
                       ("signe corrig", "sign"),
                       ("ordre-2 x4 parois", "full"),
                       ("front fluide/solide", "clean")):
            gains[lab] = (base - abs(e[v])) / base * 100.0 if base else 0.0
        best = max(gains, key=gains.get)
        r["gains_pct"] = gains
        r["source_dominante"] = best
        print("  %-5d %-5d %-4g %-6.2f |%+.3f%% %+.3f%% %+.3f%% %+.3f%% %+.3f%% %+.3f%%"
              "  | %s (%+.0f%%)"
              % (r["Re"], r["N"], r["n"], r["alpha"],
                 e["ref"], e["p2"], e["lid"], e["sign"], e["full"], e["clean"],
                 best, gains[best]))
        report["cases"].append(r)

    print()
    print("=" * 106)
    print("GAIN MOYEN EN |defaut| QUAND ON APPLIQUE CHAQUE CORRECTION")
    print("=" * 106)
    for lab in ("ordre-2 lid", "lid coherente", "signe corrig",
                "ordre-2 x4 parois", "front fluide/solide"):
        g = [r["gains_pct"][lab] for r in rows]
        print("  %-22s  moyen %+7.1f%%   (min %+7.1f%%, max %+7.1f%%)"
              % (lab, float(np.mean(g)), float(np.min(g)), float(np.max(g))))
        report["by_source"][lab] = {"mean": float(np.mean(g)),
                                    "min": float(np.min(g)),
                                    "max": float(np.max(g))}

    print()
    print("Poids des couches parietales dans D* (laterales + lid + plafond) :")
    print("  Re    N     n    alpha   bande de paroi   interieur")
    for r in rows:
        reg = r["regions"]
        print("  %-5d %-5d %-4g %-6.2f    %6.2f %%         %6.2f %%"
              % (r["Re"], r["N"], r["n"], r["alpha"],
                 reg["wall_band"] / r["D"]["ref"] * 100.0,
                 reg["interior"] / r["D"]["ref"] * 100.0))

    print()
    print("Contributions absolues des variantes de D* (en % de D_ref) :")
    print("  Re    N     n    alpha   dD(signe)  dD(clean)  dD(lid)"
          "   dD(ordre2)  dP(lid)")
    for r in rows:
        print("  %-5d %-5d %-4g %-6.2f   %+.3f     %+.3f     %+.3f      %+.3f"
              "      %+.3f"
              % (r["Re"], r["N"], r["n"], r["alpha"],
                 r["dD_sign_pct"], r["dD_clean_pct"], r["dD_lidonly_pct"],
                 r["dD_ord2_pct"], r["dP_ord2_pct"]))

    # --- differentialle alpha : ce qui peut polluer le contraste ------------
    print()
    print("DEPENDANCE DIFFERENTIELLE A L'ONDULATION (alpha=0.10 vs 0) :")
    print("  Seules ces differences polluent le contraste I ; tout le reste")
    print("  est commun aux 4 configurations et s'annule.")
    print("  Re    N     n    d(eps_ref)   d(eps_full)   d(D_ref)")
    idx = {(r["Re"], r["N"], r["n"], r["alpha"]): r for r in rows}
    for Re in Re_list:
        for N in N_list:
            for n in n_list:
                a, b = idx.get((Re, N, n, 0.0)), idx.get((Re, N, n, 0.10))
                if not a or not b:
                    continue
                print("  %-5d %-5d %-4g  %+.4f pt%%   %+.4f pt%%    %+.4f %%"
                      % (Re, N, n,
                         b["eps"]["ref"] - a["eps"]["ref"],
                         b["eps"]["full"] - a["eps"]["full"],
                         (b["D"]["ref"] - a["D"]["ref"]) / a["D"]["ref"] * 100.0))

    os.makedirs(AUDIT_DIR, exist_ok=True)
    out = os.path.join(AUDIT_DIR, "audit_balance.json")
    with open(out, "w", encoding="utf-8") as fh:
        json.dump(report, fh, indent=2, ensure_ascii=False)
    print("\n[ok] rapport = %s" % out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
