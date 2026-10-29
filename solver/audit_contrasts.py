# -*- coding: utf-8 -*-
"""
solver/audit_contrasts.py
=========================
Audit quantitatif des contrastes d'interaction, a partir des JSON bruts.

Verifie, sans aucune nouvelle simulation :

  A. conformite des cas (Ma sur la loi, residu <= tol, statut) ;
  B. les quatre termes de chaque contraste et leurs signes ;
  C. l'influence de la normalisation : D* est une integrale SUR LE MASQUE
     FLUIDE (diagnostics.py:255), or l'aire fluide change entre alpha=0 et
     alpha=0.10 (ecart de ~3.4 %). P* est une integrale lineaire sur la lid,
     dont le domaine est inchange. On compare donc les contrastes en
     grandeur extensive et en grandeur normalisee par l'aire.
  D. le bilan (P-D)/D : valeur, independance vis-a-vis de la geometrie,
     ordre de convergence apparent en N ;
  E. l'incertitude differentielle N=128 -> N=256 pour chaque grandeur ;
  F. coherence avec les fichiers interaction_*.json stockes.

Sorties : tracee stdout + results/audit/audit_contrasts.json
"""

from __future__ import annotations

import argparse
import glob
import json
import math
import os

SCALARS_DIR = os.path.join(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))), "results", "scalars")
AUDIT_DIR = os.path.join(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))), "results", "audit")

NS = (2.0, 4.0, 6.0, 8.0)
ALS = (0.0, 0.10)
PAIRS = ((2.0, 4.0), (4.0, 6.0), (6.0, 8.0))
# loi de vitesse de maquette retenue (note 08, section 3.1)
U_LAW_NUM = 40.0 / 3.0


def load_case(N, Re, n, al, tol=1e-8):
    """Charge le JSON d'un cas.

    Le suffixe ``_tol1e-08`` est privilegie ; a defaut on retombe sur le nom
    legacy (aucun suffixe), qui correspond aux cas initiaux resolus a 1e-9.
    """
    names = ["se_n%g_w%.2f_Nw1_Re%d_N%d%s.json"
             % (n, al, Re, N, s) for s in ("_tol1e-08", "")]
    for name in names:
        path = os.path.join(SCALARS_DIR, name)
        if os.path.exists(path):
            with open(path, encoding="utf-8") as fh:
                return json.load(fh), path
    return None, os.path.join(SCALARS_DIR, names[0])


def load_table(N, Re, tol=1e-8):
    out, missing = {}, []
    for n in NS:
        for al in ALS:
            rec, path = load_case(N, Re, n, al, tol)
            if rec is None:
                missing.append(os.path.basename(path))
            else:
                out[(n, al)] = rec
    return out, missing


def law_u(N):
    return U_LAW_NUM / N


def law_ma(N):
    return math.sqrt(3.0) * law_u(N)


def contrast(tab, key, na, nb, renormalize_area=False):
    """I = [Q(nb,a) - Q(nb,0)] - [Q(na,a) - Q(na,0)]."""

    def q(n, al):
        rec = tab[(n, al)]
        val = rec[key]
        if renormalize_area:
            val = val / rec["fluid_area"]
        return val

    da = q(na, 0.10) - q(na, 0.0)
    db = q(nb, 0.10) - q(nb, 0.0)
    return {
        "Q_na_0": q(na, 0.0), "Q_na_a": q(na, 0.10),
        "Q_nb_0": q(nb, 0.0), "Q_nb_a": q(nb, 0.10),
        "delta_a": da, "delta_b": db,
        "I": db - da,
        "I_rel_a_pct": 100.0 * (db - da) / abs(q(na, 0.0)) if q(na, 0.0) else None,
    }


def section(title):
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--Re", type=str, default="100,400,1000")
    ap.add_argument("--tol", type=float, default=1e-8)
    args = ap.parse_args()
    res_list = [float(x) for x in args.Re.split(",")]
    N_LO, N_HI = 128, 256

    report = {"Re": res_list, "tol": args.tol, "conformity": {},
              "normalization": {}, "balance": {}, "uncertainty": {},
              "consistency": {}}

    # ------------------------------------------------------------------ A
    section("A. Conformite des cas (loi de Ma + residu)")
    print("Ma attendue : N=128 -> %.6f | N=256 -> %.6f\n"
          % (law_ma(N_LO), law_ma(N_HI)))
    print("  N     Re    n    a      Ma      residu     statut   verdict")
    tabs = {}
    for N in (N_LO, N_HI):
        for Re in res_list:
            tab, missing = load_table(N, Re, args.tol)
            tabs[(N, Re)] = tab
            if missing:
                print("  %-5d %-5d      cas manquants: %s" % (N, Re, ", ".join(missing)))
            for (n, al), rec in sorted(tab.items()):
                ma, res = rec.get("Ma"), rec.get("residual")
                bad = []
                if abs(ma - law_ma(N)) > 1e-4:
                    bad.append("Ma hors-loi")
                if res is None:
                    bad.append("residu absent")
                elif res > args.tol * 1.0001:
                    bad.append("residu > tol")
                verdict = "OK" if not bad else "!! " + "; ".join(bad)
                res_s = "%.2e" % res if res is not None else "absent (legacy)"
                if bad or n in (2.0,):
                    print("  %-5d %-5d %-4g %-5.2f %-7.4f %-16s %-8s %s"
                          % (N, Re, n, al, ma, res_s, rec.get("status"), verdict))
                report["conformity"].setdefault(str(N), {}).setdefault(
                    str(int(Re)), []).append(
                    {"n": n, "alpha": al, "Ma": ma, "residual": res,
                     "status": rec.get("status"), "ok": not bad})

    # ------------------------------------------------------------------ B
    section("B. Les quatre termes du contraste et leurs signes")
    keys = (("D_star", "D*"), ("P_star", "P*"))
    for Re in res_list:
        for N in (N_LO, N_HI):
            tab = tabs[(N, Re)]
            if not tab:
                continue
            for na, nb in PAIRS:
                if (na, 0.0) not in tab or (nb, 0.10) not in tab:
                    continue
                for key, lab in keys:
                    c = contrast(tab, key, na, nb)
                    print("  Re=%-5d N=%-4d %s paire(%g,%g)  "
                          "Q(n_a,0)=%+.5f Q(n_a,a)=%+.5f  "
                          "Q(n_b,0)=%+.5f Q(n_b,a)=%+.5f  "
                          "d_a=%+.5f d_b=%+.5f  I=%+.6f  [%s]"
                          % (Re, N, lab, na, nb,
                             c["Q_na_0"], c["Q_na_a"], c["Q_nb_0"], c["Q_nb_a"],
                             c["delta_a"], c["delta_b"], c["I"],
                             "neg" if c["I"] < 0 else "POSIF"))

    # ------------------------------------------------------------------ C
    section("C. Normalisation : grandeur extensive vs aire fluide")
    print("D* et K* sont integres sur le masque fluide (diagnostics.py:255-256) ;")
    print("P* est une integrale lineaire sur la lid, dont le domaine est inchange.")
    print("L'aire fluide varie entre alpha=0 et alpha=0.10 -> testons les deux\n"
          "conventions de contraste.\n")
    print("  Re     N    paire   aire(a=0)  aire(a=.1)  dA/A      "
          "I(D*)      I(D*/A)    ecart")
    for Re in res_list:
        for N in (N_LO, N_HI):
            tab = tabs[(N, Re)]
            if not tab:
                continue
            for na, nb in PAIRS:
                if (na, 0.0) not in tab or (nb, 0.10) not in tab:
                    continue
                A0 = tab[(na, 0.0)]["fluid_area"]
                Aa = tab[(na, 0.10)]["fluid_area"]
                c_ext = contrast(tab, "D_star", na, nb)
                c_are = contrast(tab, "D_star", na, nb, renormalize_area=True)
                dA = (Aa - A0) / A0 * 100.0
                ecart = 100.0 * (abs(c_are["I"]) - abs(c_ext["I"])) / abs(c_ext["I"])
                print("  %-6d %-4d (%g,%g) %9.5f %10.5f %+7.3f%% "
                      "%+.6f  %+.6f  %+.2f%%"
                      % (Re, N, na, nb, A0, Aa, dA, c_ext["I"], c_are["I"], ecart))
                report["normalization"].setdefault(str(int(Re)), {}).setdefault(
                    str(N), []).append(
                    {"paire": [na, nb], "area_0": A0, "area_a": Aa,
                     "dA_pct": dA, "I_extensive": c_ext["I"],
                     "I_area": c_are["I"], "ecart_pct": ecart})

    # ------------------------------------------------------- C2 (modele nul)
    section("C2. Modele nul : que produirait un simple changement d'aire ?")
    print("Hypothese nulle : la densite de dissipation reste inchangee par")
    print("l'ondulation, seul le volume fluide varie -> D* = rho_D * A.")
    print("  dD*_aire = D*(n,0) * (A(a) - A(0))/A(0)")
    print("  Re     N    paire   dD*_aire(n_a) dD*_aire(n_b)  I(modele nul)"
          "   I(mesure)    rapport")
    nulls = {}
    for Re in res_list:
        for N in (N_LO, N_HI):
            tab = tabs[(N, Re)]
            if not tab:
                continue
            for na, nb in PAIRS:
                if (na, 0.0) not in tab or (nb, 0.10) not in tab:
                    continue
                pred = {}
                for n in (na, nb):
                    A0, Aa = tab[(n, 0.0)]["fluid_area"], tab[(n, 0.10)]["fluid_area"]
                    pred[n] = tab[(n, 0.0)]["D_star"] * (Aa - A0) / A0
                I_null = pred[nb] - pred[na]
                I_meas = contrast(tab, "D_star", na, nb)["I"]
                ratio = I_null / I_meas if abs(I_meas) > 1e-15 else float("inf")
                print("  %-6d %-4d (%g,%g) %+.5f      %+.5f      %+.6f    %+.6f   %+.1fx"
                      % (Re, N, na, nb, pred[na], pred[nb], I_null, I_meas, ratio))
                nulls[(Re, N, na, nb)] = (I_null, I_meas)
                report.setdefault("area_null", {}).setdefault(
                    str(int(Re)), {})[
                    "%d_%g%g" % (N, na, nb)] = {
                    "I_null": I_null, "I_measured": I_meas, "ratio": ratio}

    print("\n  Le modele nul n'explique donc pas le signe ni l'ordre de")
    print("  grandeur du contraste mesure : la reponse n'est pas un artefact")
    print("  de volume.")

    # ------------------------------------------------------------ C3 (robustesse)
    section("C3. Robustesse du signe et de la structure aux deux conventions")
    print("  Re     N    paire   I(D*)       I(D*/A)     meme signe ?")
    same_sign = 0
    total = 0
    for Re in res_list:
        for N in (N_LO, N_HI):
            tab = tabs[(N, Re)]
            if not tab:
                continue
            for na, nb in PAIRS:
                if (na, 0.0) not in tab or (nb, 0.10) not in tab:
                    continue
                ie = contrast(tab, "D_star", na, nb)["I"]
                ia = contrast(tab, "D_star", na, nb, renormalize_area=True)["I"]
                ok = (ie < 0) == (ia < 0)
                same_sign += ok
                total += 1
                print("  %-6d %-4d (%g,%g) %+.6f  %+.6f   %s"
                      % (Re, N, na, nb, ie, ia, "oui" if ok else "NON"))
    print("\n  signe conserve : %d/%d contrastes" % (same_sign, total))
    report["normalization"]["sign_preserved"] = "%d/%d" % (same_sign, total)

    # Part relative de chaque intervalle dans le changement total n=2 -> 8,
    # pour les deux conventions (le relecteur a demande de verifier que les
    # fractions restent additives -- elles le sont ici car tous les termes
    # ont le meme signe).
    print("\n  Repartition du changement total n=2->8 (part signee additive) :")
    print("  Re     N    convention  2->4      4->6      6->8      somme")
    for Re in res_list:
        for N in (N_LO, N_HI):
            tab = tabs[(N, Re)]
            if not tab:
                continue
            for lab, renorm in (("D* extensive", False), ("D*/aire", True)):
                if any((n, 0.10) not in tab for n in NS):
                    continue
                d = {n: contrast(tab, "D_star", 2.0, n, renorm)["delta_b"]
                     if n != 2.0 else contrast(tab, "D_star", 2.0, 4.0,
                                               renorm)["delta_a"]
                     for n in (2.0, 4.0, 6.0, 8.0)}
                # delta(n) = D*(n,a) - D*(n,0), recalcule proprement
                for n in NS:
                    q0 = tab[(n, 0.0)]["D_star"] / (tab[(n, 0.0)]["fluid_area"]
                                                    if renorm else 1.0)
                    qa = tab[(n, 0.10)]["D_star"] / (tab[(n, 0.10)]["fluid_area"]
                                                     if renorm else 1.0)
                    d[n] = qa - q0
                steps = [(d[4.0] - d[2.0]), (d[6.0] - d[4.0]), (d[8.0] - d[6.0])]
                tot = d[8.0] - d[2.0]
                parts = [s / tot * 100.0 for s in steps] if abs(tot) > 1e-15 else [0, 0, 0]
                print("  %-6d %-4d %-11s %+6.1f%%   %+6.1f%%   %+6.1f%%   %+.1f%%"
                      % (Re, N, lab, parts[0], parts[1], parts[2], sum(parts)))
                report.setdefault("shares", {}).setdefault(
                    str(int(Re)), {}).setdefault(lab, {})[str(N)] = {
                    "steps": steps, "total": tot, "pct": parts}

    # ------------------------------------------------------------------ D
    section("D. Bilan energetique (P-D)/D")
    print("  N     Re    n    a      (P-D)/D")
    bal = {}
    for N in (N_LO, N_HI):
        for Re in res_list:
            tab = tabs[(N, Re)]
            if not tab:
                continue
            for (n, al), rec in sorted(tab.items()):
                b = rec.get("balance")
                bal[(N, Re, n, al)] = b
                print("  %-5d %-5d %-4g %-5.2f  %+.4f %%"
                      % (N, Re, n, al, b * 100.0 if b is not None else float("nan")))
                report["balance"].setdefault(str(N), {}).setdefault(
                    str(int(Re)), {})["%g_%g" % (n, al)] = b

    print("\n  (a) independance vis-a-vis de la geometrie a N fixe :")
    for N in (N_LO, N_HI):
        for Re in res_list:
            vals = [v for (Nk, Rek, _n, _al), v in bal.items()
                    if Nk == N and Rek == Re and v is not None]
            if len(vals) < 2:
                continue
            spread = (max(vals) - min(vals)) * 100.0
            mean = sum(vals) / len(vals) * 100.0
            print("      N=%-4d Re=%-5d  (P-D)/D = %+.3f .. %+.3f %%"
                  "   moy=%+.3f pt%%   dispersion=%.3f pt%%   "
                  "disp/moy=%.2f%%"
                  % (N, Re, min(vals) * 100, max(vals) * 100, mean,
                     spread, 100.0 * spread / abs(mean) if mean else float("nan")))
            report.setdefault("balance_spread", {}).setdefault(
                str(N), {})[str(int(Re))] = {
                "min": min(vals), "max": max(vals), "mean": sum(vals) / len(vals),
                "spread_pp": spread}

    print("\n  (a-bis) dependance systematique du defaut a l'ondulation :")
    print("          (moyenne sur n des 4 formes, alpha=0 vs alpha=0.10)")
    print("  N     Re       <b>(a=0)    <b>(a=.1)   delta")
    for N in (N_LO, N_HI):
        for Re in res_list:
            g0 = [v for (Nk, Rek, _n, al), v in bal.items()
                  if Nk == N and Rek == Re and al == 0.0 and v is not None]
            g1 = [v for (Nk, Rek, _n, al), v in bal.items()
                  if Nk == N and Rek == Re and al == 0.10 and v is not None]
            if not g0 or not g1:
                continue
            m0, m1 = sum(g0) / len(g0), sum(g1) / len(g1)
            print("  %-5d %-6d %+.4f%%  %+.4f%%  %+.4f pt%%  "
                  "(%.1f%% du defaut)"
                  % (N, Re, m0 * 100, m1 * 100, (m1 - m0) * 100,
                     100.0 * (m1 - m0) / m0 if m0 else float("nan")))
            report.setdefault("balance_alpha", {}).setdefault(
                str(N), {})[str(int(Re))] = {"mean_a0": m0, "mean_a010": m1,
                                             "delta": m1 - m0}

    print("\n  (b) ordre de convergence apparent du defaut D-P :")
    print("      (defaut mesure : ratio N=128/N=256 = 1 -> ordre 0, 2 -> ordre 1)")
    for Re in res_list:
        b1 = bal.get((N_LO, Re, 2.0, 0.0))
        b2 = bal.get((N_HI, Re, 2.0, 0.0))
        if b1 is None or b2 is None or not b2:
            continue
        ratio = b1 / b2
        order = math.log(abs(ratio), 2) if abs(ratio) > 1 else float("nan")
        print("      Re=%-5d  %+.4f -> %+.4f   ratio=%.2f   ordre apparent~%.2f"
              % (Re, b1 * 100, b2 * 100, ratio, order))
        report["balance"].setdefault("convergence", {})[str(int(Re))] = {
            "bal_128": b1, "bal_256": b2, "ratio": ratio, "order": order}

    # ------------------------------------------------------------------ E
    section("E. Incertitude differentielle N=128 -> N=256")
    print("  ecart = (|I_128| - |I_256|) / |I_256|   (ref = grille fine)")
    print("  Re     paire   grandeur  I(N=128)     I(N=256)     ecart")
    for Re in res_list:
        for na, nb in PAIRS:
            for key, lab in (("D_star", "D*"), ("P_star", "P*")):
                t1 = tabs.get((N_LO, Re))
                t2 = tabs.get((N_HI, Re))
                if not t1 or not t2:
                    continue
                needed = ((na, 0.0), (na, 0.10), (nb, 0.0), (nb, 0.10))
                if any(k not in t1 for k in needed) or any(k not in t2 for k in needed):
                    continue
                i1 = contrast(t1, key, na, nb)["I"]
                i2 = contrast(t2, key, na, nb)["I"]
                e = (abs(i2) - abs(i1)) / abs(i2) * 100.0
                flag = "  <10%" if abs(e) < 10 else "  >10% !!"
                print("  %-6d (%g,%g) %-9s %+.6f   %+.6f   %+.2f%%%s"
                      % (Re, na, nb, lab, i1, i2, e, flag))
                report["uncertainty"].setdefault(str(int(Re)), {})[
                    "%s_%g%g" % (key, na, nb)] = {
                    "I_128": i1, "I_256": i2, "ecart_pct": e}

    # --------------------------------------------------------------- E2
    section("E2. Controle croise P* vs D* : barre d'erreur independante")
    print("  I_P et I_D sont calcules sur des domaines distincts (lid vs")
    print("  masque fluide) : leur ecart mesure l'incertitude reelle, qui")
    print("  couvre aussi les erreurs de modele, pas seulement la tolerance.")
    print("  Re     paire   I_D(N=256)   I_P(N=256)   I_P/I_D - 1")
    for Re in res_list:
        t2 = tabs.get((N_HI, Re))
        if not t2:
            continue
        for na, nb in PAIRS:
            try:
                d = contrast(t2, "D_star", na, nb)["I"]
                p = contrast(t2, "P_star", na, nb)["I"]
            except KeyError:
                continue
            if abs(d) < 1e-15:
                continue
            r = p / d - 1.0
            print("  %-6d (%g,%g) %+.6f   %+.6f   %+.2f%%"
                  % (Re, na, nb, d, p, r * 100.0))
            report.setdefault("crosscheck", {}).setdefault(
                str(int(Re)), {})["%g%g" % (na, nb)] = {
                "I_D": d, "I_P": p, "rel_diff_pct": r * 100.0}

    # ------------------------------------------------------------------ F
    section("F. Cohérence avec les fichiers interaction_*.json stockes")
    for N in (N_LO, N_HI):
        for Re in res_list:
            path = None
            for suf in ("_tol1e-08", ""):
                cand = os.path.join(SCALARS_DIR,
                                    "interaction_Re%g_N%d%s.json" % (Re, N, suf))
                if os.path.exists(cand):
                    path = cand
                    break
            if path is None:
                continue
            with open(path, encoding="utf-8") as fh:
                stored = json.load(fh)
            tab = tabs.get((N, Re))
            if not tab:
                continue
            worst = 0.0
            for pname, pv in stored.get("pairs", {}).items():
                na = float(pname.split("_")[0][1:])
                nb = float(pname.split("_")[1][1:])
                if (na, 0.0) not in tab:
                    continue
                fresh = contrast(tab, "D_star", na, nb)["I"]
                old = pv["D_star"]["I"]
                if abs(old) > 1e-12:
                    worst = max(worst, abs(fresh - old) / abs(old) * 100)
            print("  N=%-4d Re=%-5d  %d paires, ecart max I_D vs recalcule = %.4f %%"
                  % (N, Re, len(stored.get("pairs", {})), worst))
            report["consistency"]["%d_%d" % (N, Re)] = worst

    os.makedirs(AUDIT_DIR, exist_ok=True)
    out = os.path.join(AUDIT_DIR, "audit_contrasts.json")
    with open(out, "w", encoding="utf-8") as fh:
        json.dump(report, fh, indent=2, ensure_ascii=False)
    print("\n[ok] rapport = %s" % out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
