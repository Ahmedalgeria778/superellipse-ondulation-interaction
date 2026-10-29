# -*- coding: utf-8 -*-
"""
solver/plot_resweep.py
======================
Cartographie de la dependance en Reynolds du contraste d'interaction.

Charge les tableaux `results/scalars/pilot_Re*_N{N}.csv`, recalcule le
contraste apparié (n=2, alpha=0) vs (n=4, alpha=0.10) pour D*, P*, K*, puis
compare trois modeles de la dependance en Re :

  M1  loi de puissance      |I| = A * Re^gamma
  M2  saturation            |I| = A * Re^gamma / (1 + Re/Re_s)
  M3  lineaire              |I| = A * (Re - Re_0)

avec validation croisee leave-one-Re-out (LOGO) et selection par RMSE.

Usage
-----
  python -m solver.plot_resweep --N 128
  python -m solver.plot_resweep --N 128 --N-ref 256
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import re
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
for _p in (_HERE, _ROOT):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import numpy as np                                     # noqa: E402
import matplotlib                                      # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt                        # noqa: E402

SCALARS_DIR = os.path.join(_ROOT, "results", "scalars")
FIG_DIR = os.path.join(_ROOT, "figures", "analysis")
OUT_DIR = os.path.join(_ROOT, "results", "resweep")

# Stabilite relative N=128 -> N=256 mesuree sur I_D a Re = 100/400/1000
# (voir docs/06_tests_tolerance_resolution.md) : utilisee comme barre
# d'erreur de reference pour le graphe.
STAB_PCT = {"100.0": 3.53, "400.0": 1.96, "1000.0": 2.83}
STAB_DEFAULT_PCT = 3.53

PAIR = (2.0, 4.0)


def load_table(N, Re):
    """Retourne ({alpha: {n: record}}, chemin) pour un cas donne, ou None.

    Priorite aux JSON individuels : ils portent le `Ma` reellement simule, ce
    qui evite qu'un CSV herite (calcule a une autre vitesse de maquette)
    ecrase une campagne plus recente. Le CSV ne sert que de repli.
    """
    import csv

    for suf in ("_tol1e-08", ""):
        out, n_found = {}, 0
        for n in (2.0, 4.0, 6.0, 8.0):
            for al in (0.0, 0.10):
                cid = "se_n%g_w%.2f_Nw1_Re%d_N%d%s" % (n, al, int(Re), N, suf)
                path = os.path.join(SCALARS_DIR, cid + ".json")
                if not os.path.exists(path):
                    continue
                try:
                    with open(path, encoding="utf-8") as fh:
                        rec = json.load(fh)
                except Exception:
                    continue
                out.setdefault(al, {})[n] = rec
                n_found += 1
        if n_found >= 4:
            return out, "(JSON %s, %d cas)" % (suf or "legacy", n_found)

    for pat in ("pilot_Re%g_N%d_tol1e-08.csv" % (Re, N),
                "pilot_Re%g_N%d.csv" % (Re, N)):
        path = os.path.join(SCALARS_DIR, pat)
        if not os.path.exists(path):
            continue
        out = {}
        with open(path, encoding="utf-8") as fh:
            for row in csv.DictReader(fh):
                out.setdefault(float(row["alpha"]), {})[float(row["n"])] = row
        return out, os.path.basename(path)
    return None, None


def contrast(tbl, key):
    """Contraste double difference appariée.

    a = (n=2, alpha=0)   b = (n=4, alpha=0.10)
    I = [D(b,0.10) - D(b,0)] - [D(a,0.10) - D(a,0)]
    soit exactement la grandeur `pairs.n2_n4_a0.100` des JSON d'interaction.
    """

    def val(n, al):
        r = tbl.get(al, {}).get(n)
        return float(r[key]) if r is not None else None

    a0, a1 = val(PAIR[0], 0.0), val(PAIR[0], 0.10)
    b0, b1 = val(PAIR[1], 0.0), val(PAIR[1], 0.10)
    if None in (a0, a1, b0, b1):
        return None
    da = a1 - a0
    db = b1 - b0
    I = db - da
    return {"I": I, "delta_a": da, "delta_b": db,
            "val_a": a0, "val_b": b0,
            "I_rel_a_pct": 100.0 * I / a0}


def discover(N):
    """Liste des Re disponibles pour une grille donnee (CSV ou JSON)."""
    Res = set()
    for path in glob.glob(os.path.join(SCALARS_DIR, "pilot_Re*_N%d*.csv" % N)):
        m = re.search(r"pilot_Re([0-9.]+)_N%d" % N, os.path.basename(path))
        if m:
            Res.add(float(m.group(1)))
    for path in glob.glob(os.path.join(SCALARS_DIR, "se_*_N%d*.json" % N)):
        m = re.search(r"_Re([0-9]+)_N%d" % N, os.path.basename(path))
        if m:
            Res.add(float(m.group(1)))
    return sorted(Res)


# ---------------------------------------------------------------- modeles
def _fit_power(Re, y):
    b, c = np.polyfit(np.log(Re), np.log(y), 1)
    return {"gamma": float(b), "A": float(np.exp(c)), "n": 2}


def _fit_saturation(Re, y):
    """|I| = A * Re^gamma / (1 + Re/Re_s), Re_s par grille log-spacée."""
    best = None
    for Re_s in np.logspace(np.log10(20.0), np.log10(4000.0), 400):
        x = Re / (1.0 + Re / Re_s)
        try:
            b, c = np.polyfit(np.log(x), np.log(y), 1)
        except Exception:
            continue
        pred = np.exp(c) * x ** b
        rms = float(np.sqrt(np.mean((np.log(pred) - np.log(y)) ** 2)))
        if best is None or rms < best[0]:
            best = (rms, {"gamma": float(b), "A": float(np.exp(c)),
                          "Re_s": float(Re_s), "n": 3})
    return best[1] if best else None


def _fit_linear(Re, y):
    b, c = np.polyfit(Re, y, 1)
    Re0 = -c / b if b != 0 else float("nan")
    return {"A": float(b), "Re_0": float(Re0), "n": 2}


MODELS = {
    "power": (_fit_power, lambda p, R: p["A"] * R ** p["gamma"]),
    "saturation": (_fit_saturation,
                   lambda p, R: p["A"] * R ** p["gamma"] / (1.0 + R / p["Re_s"])),
    "linear": (_fit_linear, lambda p, R: p["A"] * (R - p["Re_0"])),
}


def logo_cv(Re, y, fit_fn, pred_fn, n_refit=1):
    """Leave-one-Re-out ; RMSE sur log|I| (echelle relative)."""
    if len(Re) < 4:
        return None
    errs = []
    for i in range(len(Re)):
        m = np.ones(len(Re), bool)
        m[i] = False
        if m.sum() < n_refit + 1:
            return None
        try:
            p = fit_fn(Re[m], y[m])
        except Exception:
            return None
        if p is None:
            return None
        pred = pred_fn(p, Re[i])
        if pred <= 0 or not np.isfinite(pred):
            return None
        errs.append(np.log(pred) - np.log(y[i]))
    e = np.array(errs)
    return {"rmse_log": float(np.sqrt(np.mean(e ** 2))),
            "bias_log": float(np.mean(e)), "n": int(len(e)),
            "rmse_pct": float(100.0 * (np.exp(np.sqrt(np.mean(e ** 2))) - 1.0))}


def analyse_series(Re, y, label):
    out = {"label": label, "Re": Re.tolist(), "I": y.tolist(), "models": {}}
    for name, (fit_fn, pred_fn) in MODELS.items():
        try:
            p = fit_fn(np.array(Re), y)
        except Exception:
            p = None
        if p is None:
            continue
        pred = pred_fn(p, np.array(Re))
        p = dict(p)
        p["pred"] = np.asarray(pred, float).tolist()
        if np.all(pred > 0) and np.all(y > 0):
            resid = np.log(pred) - np.log(y)
            p["r2_log"] = float(1.0 - np.sum(resid ** 2)
                                / np.sum((np.log(y) - np.log(y).mean()) ** 2))
            p["cv"] = logo_cv(np.array(Re), y, fit_fn, pred_fn)
        else:
            p["r2_log"] = float("nan")
            p["cv"] = None
        out["models"][name] = p
    if out["models"]:
        scored = [(m["cv"]["rmse_log"], n)
                  for n, m in out["models"].items() if m.get("cv")]
        if scored:
            out["best"] = min(scored)[1]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int, default=128)
    ap.add_argument("--N-ref", type=int, default=256)
    ap.add_argument("--key", type=str, default="D_star")
    args = ap.parse_args()

    os.makedirs(FIG_DIR, exist_ok=True)
    os.makedirs(OUT_DIR, exist_ok=True)

    Res = discover(args.N)
    print("[resweep] N=%d : %d Reynolds disponibles : %s"
          % (args.N, len(Res), ", ".join("%g" % r for r in Res)), flush=True)

    rows = []
    for Re in Res:
        tbl, _src = load_table(args.N, Re)
        if tbl is None:
            continue
        c_abs = contrast(tbl, args.key)
        if c_abs is None:
            continue
        it = None
        try:
            it = max(int(float(r["iters"]))
                     for r in tbl.values() for r in r.values()
                     if r.get("iters") not in (None, "", "None"))
        except Exception:
            it = None
        rows.append({"Re": Re, "I": c_abs["I"],
                     "I_rel_a_pct": c_abs["I_rel_a_pct"],
                     "iters": it, "src": _src})

    if len(rows) < 3:
        print("[resweep] moins de 3 Re exploitables, arret", flush=True)
        return 1

    Re_a = np.array([r["Re"] for r in rows])
    y_signed = np.array([r["I"] for r in rows])
    y_abs = np.abs(y_signed)
    y_rel = np.abs(np.array([r["I_rel_a_pct"] for r in rows]))

    res = {"N": args.N, "key": args.key,
           "pair": {"a": {"n": PAIR[0], "alpha": 0.0},
                    "b": {"n": PAIR[1], "alpha": 0.10}},
           "points": rows,
           "abs": analyse_series(Re_a, y_abs, "|I|"),
           "rel": analyse_series(Re_a, y_rel, "|I_rel_a|")}

    # barres d'erreur : stabilite N->N-ref mesuree sur les 3 Re de reference
    stab = []
    for Re in Re_a:
        stab.append(STAB_PCT.get(str(Re), STAB_DEFAULT_PCT))
    res["stab_pct"] = stab

    # comparaison N vs N-ref si dispo
    if args.N_ref:
        ref = []
        for Re in Re_a:
            tbl, _ = load_table(args.N_ref, Re)
            if tbl is None:
                continue
            c = contrast(tbl, args.key)
            if c is None:
                continue
            ref.append({"Re": Re, "I_ref": c["I"],
                        "ecart_pct": 100.0 * (abs(c["I"]) - y_abs[list(Re_a).index(Re)])
                        / abs(c["I"])})
        res["vs_ref"] = ref

    out_json = os.path.join(OUT_DIR, "resweep_N%d.json" % args.N)
    with open(out_json, "w", encoding="utf-8") as fh:
        json.dump(res, fh, indent=2, ensure_ascii=False)

    # ---------------------------------------------------------------- rapport
    print("\n[resweep] contraste %s  paire (n=%g,a=0) vs (n=%g,a=0.10)"
          % (args.key, PAIR[0], PAIR[1]))
    print("  %-8s %-14s %-12s %s" % ("Re", "I", "I_rel_a_%", "iters"))
    for r in rows:
        print("  %-8g %-14.6g %-12.4f %s"
              % (r["Re"], r["I"], r["I_rel_a_pct"], r["iters"] or "?"))

    for tag, ser in (("absolu", res["abs"]), ("relatif", res["rel"])):
        print("\n  --- %s ---" % tag)
        if "best" in ser:
            print("  meilleur modele (CV LOGO) : %s" % ser["best"])
        for name, m in ser["models"].items():
            cv = m.get("cv") or {}
            print("    %-11s n=%d R2_log=%+.4f  %s"
                  % (name, m["n"], m["r2_log"],
                     ("CV RMSE=%.2f%% bias=%+.2f%% (n=%d)"
                      % (cv.get("rmse_pct", float("nan")),
                         100.0 * cv.get("bias_log", 0.0), cv.get("n", 0))
                      if cv else "CV indisponible")))
            if name == "power":
                print("                gamma = %+.4f  A = %.5g"
                      % (m["gamma"], m["A"]))
            if name == "saturation":
                print("                gamma = %+.4f  Re_s = %.4g  A = %.5g"
                      % (m["gamma"], m["Re_s"], m["A"]))
            if name == "linear":
                print("                A = %.5g  Re_0 = %.4g" % (m["A"], m["Re_0"]))

    if res.get("vs_ref"):
        print("\n  --- confrontation N=%d vs N=%d ---" % (args.N, args.N_ref))
        for r in res["vs_ref"]:
            print("    Re=%-6g I(N=%d)=%+.6g  I(N=%d)=%+.6g  ecart=%+.2f%%"
                  % (r["Re"], args.N_ref, r["I_ref"], args.N,
                     y_signed[list(Re_a).index(r["Re"])],
                     r["ecart_pct"]))

    # ---------------------------------------------------------------- figures
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.6))
    ax = axes[0]
    ax.errorbar(Re_a, y_rel, yerr=np.asarray(stab) * y_rel / 100.0,
                fmt="o", color="C3", ms=7, lw=1.4, capsize=3.5,
                label=r"$|I_{rel}|$ (N=%d)" % args.N)
    if args.N_ref:
        rr = res.get("vs_ref") or []
        if rr:
            ax.plot([r["Re"] for r in rr],
                    [abs(y_rel[list(Re_a).index(r["Re"])]) for r in rr],
                    "s", mfc="none", mec="C0", ms=7, mew=1.4,
                    label=r"$|I_{rel}|$ (N=%d)" % args.N_ref)
    best = res["rel"].get("best")
    if best:
        m = res["rel"]["models"][best]
        Rg = np.linspace(Re_a.min() * 0.9, Re_a.max() * 1.1, 300)
        _, pred_fn = MODELS[best]
        ax.plot(Rg, pred_fn(m, Rg), "-", color="0.4", lw=1.4,
                label="%s (CV %.1f%%)" % (best, m["cv"]["rmse_pct"]))
    ax.set_xscale("log")
    ax.set_xlabel("Re")
    ax.set_ylabel(r"$|I_{rel}|$  (%)")
    ax.set_title(r"Contraste d'interaction $(2,0)\to(4,0.10)$")
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=8)

    ax = axes[1]
    ax.errorbar(np.log(Re_a), np.log(y_rel),
                yerr=[np.log(1.0 + p / 100.0) for p in stab],
                fmt="o", color="C3", ms=7, lw=1.4, capsize=3.5)
    if best:
        Rg = np.linspace(Re_a.min() * 0.9, Re_a.max() * 1.1, 300)
        _, pred_fn = MODELS[best]
        ax.plot(np.log(Rg), np.log(pred_fn(m, Rg)), "-", color="0.4", lw=1.4)
        if best == "power":
            lab = r"$\gamma=%.3f$" % m["gamma"]
        elif best == "saturation":
            lab = r"$\gamma=%.3f,\ Re_s=%.0f$" % (m["gamma"], m["Re_s"])
        else:
            lab = r"lineaire, $Re_0=%.0f$" % m["Re_0"]
        ax.text(0.03, 0.95, lab,
                transform=ax.transAxes, va="top", fontsize=9,
                bbox=dict(fc="white", alpha=0.8, ec="0.7"))
    ax.set_xlabel("ln Re")
    ax.set_ylabel(r"ln $|I_{rel}|$")
    ax.set_title("log-log : loi de puissance ?")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    out_png = os.path.join(FIG_DIR, "resweep_ID_vs_Re_N%d.png" % args.N)
    fig.savefig(out_png, dpi=160)
    plt.close(fig)

    print("\n[resweep] ecrit : %s" % os.path.relpath(out_json, _ROOT))
    print("[resweep] ecrit : %s" % os.path.relpath(out_png, _ROOT))
    return 0


if __name__ == "__main__":
    sys.exit(main())
