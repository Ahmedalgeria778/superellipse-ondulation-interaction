# -*- coding: utf-8 -*-
"""
solver/fit_mach_resolution.py
=============================
Diagnostic de sensibilite eps = f(h, Ma) sur les 8 points de Balayage Ma :

    Re=400, n=2, alpha=0 :  4 valeurs de Ma a N=128 (h=1/128)
                            4 valeurs de Ma a N=256 (h=1/256)

Le dessein est factoriel complet (2 niveaux de h x 4 niveaux de Ma), donc les
effets principaux de h et de Ma sont estimables sans correlation par
construction ; le terme h*Ma (interaction) est lui aussi identifiable.

MODELES COMPARES
    M0h  eps = b0 + bh*h                      (grille seule)
    M0M  eps = b0 + bM*Ma                     (compressibilite seule)
    M1   eps = b0 + bh*h + bM*Ma              (ADDITIF)
    M2   eps = b0 + bh*h + bM*Ma + bhM*h*Ma   (AVEC INTERACTION)
    M3   eps = b0 + bh*h + bM*Ma^2            (additif, Ma^2)
    M4   eps = b0 + bh*h + bM*Ma + bM2*Ma^2   (additif, Ma + Ma^2)

RAPPEL DU COUPLAGE (a ecrire explicitement dans le texte)
    A N,Re fixes :  tau = 0.5 + sqrt(3)*N*Ma/Re = 0.5 + sqrt(3)*Ma/(Re*h).
    Ma, tau et N sont donc lies dans cette base : ce balayage NE LES SEPARERA
    PAS.  Seul le croisement des deux resolutions rend h et Ma estimables
    l'un sans l'autre (design factoriel), pas tau.

Lecture : CECI EST UN DIAGNOSTIC DE SENSIBILITE, PAS UNE LOI GENERALE.
Un seul Reynolds (400), une seule forme (n=2, alpha=0), 8 points, 2 niveaux
de h : aucun ajustement multi-parametre ne peut etre presente comme une
asymptotique.

Sortie console + results/ma_sweep/fit_h_ma.json
"""

from __future__ import annotations

import json
import math
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
MA_DIR = os.path.join(_ROOT, "results", "ma_sweep")

MODELS = [
    ("M0h", ["h"], "grille seule"),
    ("M0M", ["Ma"], "compressibilite seule"),
    ("M1", ["h", "Ma"], "ADDITIF"),
    ("M2", ["h", "Ma", "hMa"], "avec INTERACTION h*Ma"),
    ("M3", ["h", "Ma2"], "additif, Ma^2"),
    ("M4", ["h", "Ma", "Ma2"], "additif, Ma + Ma^2"),
]


def load_points():
    pts = []
    for N in (128, 256):
        path = os.path.join(MA_DIR, "masweep_Re400_N%d.json" % N)
        if not os.path.exists(path):
            print("MANQUANT : %s" % path, file=sys.stderr)
            return None
        with open(path, "r", encoding="utf-8") as fh:
            doc = json.load(fh)
        for r in doc["cases"]:
            pts.append({
                "N": int(r["N"]), "h": 1.0 / float(r["N"]),
                "Ma": float(r["Ma"]), "tau": float(r["tau"]),
                "U_ref": float(r["U_ref"]),
                "D_star": float(r["D_star"]), "P_star": float(r["P_star"]),
                "balance_pct": 100.0 * float(r["balance"]),
                "iters": r.get("iters"), "residual": r.get("residual"),
                "status": r.get("status"),
            })
    pts.sort(key=lambda p: (p["N"], p["Ma"]))
    return pts


def design(pts, terms):
    rows = []
    for p in pts:
        row = [1.0]
        for t in terms:
            if t == "h":
                row.append(p["h"])
            elif t == "Ma":
                row.append(p["Ma"])
            elif t == "Ma2":
                row.append(p["Ma"] ** 2)
            elif t == "hMa":
                row.append(p["h"] * p["Ma"])
            else:
                raise ValueError(t)
        rows.append(row)
    return np.array(rows, dtype=float)


def fit(pts, terms):
    y = np.array([p["balance_pct"] for p in pts], dtype=float)
    X = design(pts, terms)
    n, k = X.shape
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    ss_res = float(resid @ resid)
    ss_tot = float(((y - y.mean()) ** 2).sum())
    dof = n - k
    s2 = ss_res / dof if dof else float("nan")
    try:
        cov = s2 * np.linalg.inv(X.T @ X)
        se = np.sqrt(np.diag(cov))
    except np.linalg.LinAlgError:
        se = np.full(k, float("nan"))
    r2 = 1.0 - ss_res / ss_tot if ss_tot else float("nan")
    aicc = (n * math.log(ss_res / n) + 2 * k
            + 2.0 * k * (k + 1) / (n - k - 1)) if n - k - 1 > 0 else float("nan")
    return {"terms": terms, "beta": beta, "se": se, "resid": resid,
            "rmse": math.sqrt(ss_res / n), "max_abs_resid": float(np.abs(resid).max()),
            "r2": r2, "aicc": aicc, "k": k, "n": n, "dof": dof}


def main():
    pts = load_points()
    if not pts:
        return 1

    print("=" * 78)
    print("Diagnostic de sensibilite  eps = f(h, Ma)   |   8 points (Re=400, n=2, a=0)")
    print("=" * 78)
    print("  N     h          tau      Ma       D*        P*        eps       it/status")
    for p in pts:
        print("  %-5d %-10.6f %-8.4f %-8.4f %-9.5f %-9.5f %+8.4f%%  %s/%s"
              % (p["N"], p["h"], p["tau"], p["Ma"], p["D_star"], p["P_star"],
                 p["balance_pct"], p["iters"], p["status"]))

    print()
    print("  Couplage : tau = 0.5 + sqrt(3)*Ma/(Re*h)  (Re=400 fixe)  ->  Ma, tau et N")
    print("  ne sont PAS separables ; seul le croisement 2 x 4 rend h et Ma estimables.")
    print()

    results = {}
    print("  %-5s %-3s %-9s %-9s %-9s %-9s %-9s %-8s %-7s %-8s %s"
          % ("mod", "k", "b0[%]", "bh", "bMa", "bMa2", "bhMa",
             "RMSE", "R2", "AICc", "max|r|"))
    print("  " + "-" * 104)
    for name, terms, label in MODELS:
        r = fit(pts, terms)
        results[name] = {
            "label": label, "terms": terms, "k": r["k"], "n": r["n"],
            "beta": r["beta"].tolist(), "se": r["se"].tolist(),
            "rmse_pct": r["rmse"], "max_abs_resid_pct": r["max_abs_resid"],
            "r2": r["r2"], "aicc": r["aicc"],
            "resid_pct": r["resid"].tolist(),
        }
        vals = {t: float("nan") for t in ("h", "Ma", "Ma2", "hMa")}
        for i, t in enumerate(terms):
            vals[t] = r["beta"][i + 1]
        print("  %-5s %-3d %+9.4f %+9.4f %+9.4f %+9.4f %+9.4f %-8.4f %-7.4f %+8.1f %+.4f"
              % (name, r["k"], r["beta"][0], vals["h"], vals["Ma"], vals["Ma2"],
                 vals["hMa"], r["rmse"], r["r2"], r["aicc"], r["max_abs_resid"]))
    print("  unites : bh en %/h, bMa et bMa2 en %/Ma(^p), bhMa en %/(h*Ma) ;"
          " le design 2x4 (Ma identiques aux deux h) rend h orthogonal a tout"
          " polynome en Ma -> bh est le meme dans tous les modeles")

    b0 = {m: results[m]["beta"][0] for m in results}
    lo, hi = min(b0.values()), max(b0.values())
    print()
    print("  Intercept b0 (defaut extrapole a h->0, Ma->0) selon le modele :")
    for m in MODELS:
        print("    %-4s %+8.4f %%   (ecart au point Ma=0,0902 : %+7.4f pt%%)"
              % (m[0], b0[m[0]], b0[m[0]] - min(
                  p["balance_pct"] for p in pts if p["N"] == 256)))
    print("    -> amplitude du choix de modele : %.4f pt%%  (%s)"
          % (hi - lo, "sensibilite au modele, pas loi"))

    ref = [p["balance_pct"] for p in pts if p["N"] == 128
           and abs(p["Ma"] - math.sqrt(3.0) * (40.0 / 3.0) / 128.0) < 1e-9]
    ref128 = ref[0] if ref else float("nan")
    print()
    print("  b0 rapporte au defaut de campagne N=128 (Ma=0,1804) :")
    for m in MODELS:
        print("    %-4s %6.1f %%" % (m[0], b0[m[0]] / ref128 * 100.0))

    print()
    print("  Residus par point (pt%%) :")
    print("    %-5s %-10s %s" % ("mod", "N", " ".join("%7.4f" % p["Ma"] for p in pts)))
    for name in results:
        r = results[name]
        print("    %-5s %-10s %s" % (name, "tous",
                                     " ".join("%+7.4f" % v for v in r["resid_pct"])))

    out = {
        "Re": 400, "n": 2, "alpha": 0.0,
        "design": "factoriel complet 2 (h) x 4 (Ma)",
        "coupling": "tau = 0.5 + sqrt(3)*Ma/(Re*h) : Ma, tau et N lies",
        "points": pts, "models": results,
        "b0_by_model": b0, "b0_span_pct_pt": hi - lo,
        "warning": ("diagnostic de sensibilite sur 8 points, un seul Re, une seule "
                    "forme, 2 niveaux de h ; pas une loi generale"),
    }
    path = os.path.join(MA_DIR, "fit_h_ma.json")
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=2, ensure_ascii=False)
    print("\n[ok] %s" % path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
