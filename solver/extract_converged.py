# -*- coding: utf-8 -*-
"""
solver/extract_converged.py
===========================
Extraction des resultats converges et regeneration des lignes de tableaux.

Ce script est l'unique source des chiffres de la table SM.X : il relit les JSON
bruts, recalcule les contrastes et reimprime les lignes. Il ne modifie jamais
le manuscrit : les propositions sont ecrites dans results/extract/.

Definitions (figees dans le manuscrit) :
  I_D(qi,qj) = 100 * (D*_j - D*_i) / D*_i            [%)
  I_P(qi,qj) = 100 * (P*_j - P*_i) / P*_i            [%)
  delta      = 100 * |I_P/I_D - 1|                   [%]   ( "--" si |I_D|<1e-4 )
  diff       = I_P - I_D                             [pp]
  I_area     = 100 * (Af_j - Af_i) / Af_i            [%]   (prediction extensive)
  span       = max(I_D) - min(I_D) sur les N disponibles   [pp]

Usage :
  python -m solver.extract_converged            # extraction + validation
  python -m solver.extract_converged --report   # couverture des campagnes
"""

from __future__ import annotations

import csv
import glob
import io
import json
import os
import re
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
SCAL = os.path.join(_ROOT, "results", "scalars")
LOGS = os.path.join(_ROOT, "results", "logs")
OUT = os.path.join(_ROOT, "results", "extract")

TABLE_TEX = os.path.join(_ROOT, "paper", "supplementary",
                         "contrasts_table.tex")

EXPO = (2, 4, 6, 8)
TRANS = ((2, 4), (4, 6), (6, 8))
RES = (100, 400, 1000)
ALPHAS = (0.0, 0.1)
N_LEVELS = (128, 256, 384, 512)

# cibles des campagnes encore ouvertes
TARGETS = {
    "P2.1 (N=512, alpha=0)": [
        dict(n=n, alpha=0.0, Nw=1, Re=Re, N=512)
        for Re in (400, 1000) for n in (2, 4)
    ] + [dict(n=n, alpha=0.0, Nw=1, Re=1000, N=512) for n in (6, 8)],
    "P3.1 (Nw=2)": [dict(n=n, alpha=0.0, Nw=2, Re=Re, N=384)
                    for Re in (400, 1000) for n in (2, 4)],
    "P3.2 (alpha=0.05/0.20)": [dict(n=n, alpha=a, Nw=1, Re=1000, N=384)
                              for a in (0.05, 0.20) for n in (2, 4)],
}


def case_id(n, alpha, Nw, Re, N):
    return "se_n%d_w%.2f_Nw%d_Re%d_N%d_tol1e-08" % (n, alpha, Nw, Re, N)


def ok(rec):
    return rec.get("status") in ("converged", "reanalyzed")


def load_converged():
    """index (n, alpha, Nw, Re, N) -> dict de diagnostics, statuts utilisables."""
    recs = {}
    for p in sorted(glob.glob(os.path.join(SCAL, "se_*.json"))):
        try:
            d = json.load(io.open(p, encoding="utf-8"))
        except ValueError:
            continue
        if d.get("status") not in ("converged", "reanalyzed"):
            continue
        for k in ("D_star", "P_star", "fluid_area"):
            if d.get(k) is None:
                break
        else:
            key = (int(round(d["n"])), round(float(d["alpha"]), 4),
                   int(d.get("Nw", 1)), int(round(d["Re"])), int(d["N"]))
            recs[key] = dict(d, path=p)
    return recs


def contrast(a, b):
    """a = reference (petit n), b = compare (grand n)."""
    i_d = 100.0 * (b["D_star"] - a["D_star"]) / a["D_star"]
    i_p = 100.0 * (b["P_star"] - a["P_star"]) / a["P_star"]
    i_area = 100.0 * (b["fluid_area"] - a["fluid_area"]) / a["fluid_area"]
    per_area = 100.0 * ((b["D_star"] / b["fluid_area"])
                        / (a["D_star"] / a["fluid_area"]) - 1.0)
    delta = 100.0 * abs(i_p / i_d - 1.0) if abs(i_d) > 1e-12 else float("nan")
    return dict(I_D=i_d, I_P=i_p, delta=delta, diff=i_p - i_d,
                I_area=i_area, I_per_area=per_area,
                ratio=(i_area / i_d if abs(i_d) > 1e-12 else float("nan")))


def series(recs, Re, tr, alpha):
    """I_D, I_P, ... pour chaque N disponible (ordre croissant)."""
    n_i, n_j = tr
    rows = []
    for N in N_LEVELS:
        a = recs.get((n_i, alpha, 1, Re, N))
        b = recs.get((n_j, alpha, 1, Re, N))
        if not a or not b or not ok(a) or not ok(b):
            continue
        c = contrast(a, b)
        c.update(N=N, Ma_a=a.get("Ma"), Ma_b=b.get("Ma"),
                 iters_a=a.get("iters"), iters_b=b.get("iters"))
        rows.append(c)
    return rows


def fmt(v, w=8, prec=5):
    return "$%+.*f$" % (prec, v) if v == v else "--"


def tex_rows(recs, with_span_at=None):
    """Lignes au format SM.X ; span calcule sur les N disponibles."""
    out = []
    for Re in RES:
        for tr in TRANS:
            block = []
            for alpha in ALPHAS:
                rows = series(recs, Re, tr, alpha)
                if not rows:
                    continue
                span = max(r["I_D"] for r in rows) - min(r["I_D"] for r in rows)
                last = rows[-1]["N"] if with_span_at is None else with_span_at
                for r in rows:
                    d = ("--" if abs(r["I_D"]) < 1e-4
                         else "%5.2f" % r["delta"])
                    sp = ("$%+.5f$" % span) if r["N"] == last else ""
                    block.append("%-4d & %-3d & $(%d,%d)$ & %.2f & %s & %s & %s & %s & %s \\\\"
                                 % (Re, r["N"], tr[0], tr[1], alpha,
                                    fmt(r["I_D"]), fmt(r["I_P"]), d,
                                    fmt(r["diff"]), sp))
                block.append("\\hline")
            out.extend(block)
    return out


def parse_existing(path):
    """(Re, N, transition, alpha) -> (I_D, I_P, delta, diff) lus du .tex gelé."""
    if not os.path.exists(path):
        return {}
    txt = io.open(path, encoding="utf-8").read()
    got = {}
    pat = re.compile(
        r"^(\d+)\s*&\s*(\d+)\s*&\s*\$\((\d+),(\d+)\)\$\s*&\s*([\d.]+)\s*&\s*"
        r"\$([+\-][\d.]+)\$\s*&\s*\$([+\-][\d.]+)\$\s*&\s*([\d.\-]+)\s*&\s*"
        r"\$([+\-][\d.]+)\$", re.M)
    for m in pat.finditer(txt):
        Re, N, ni, nj, al = (int(m.group(1)), int(m.group(2)), int(m.group(3)),
                             int(m.group(4)), float(m.group(5)))
        got[(Re, N, (ni, nj), al)] = (
            float(m.group(6)), float(m.group(7)),
            None if m.group(8).strip() == "--" else float(m.group(8)),
            float(m.group(9)))
    return got


def validate(recs):
    """Confronte les lignes recalculees a la table SM.X gelée."""
    frozen = parse_existing(TABLE_TEX)
    if not frozen:
        return ["table gelee introuvable ou non analysee : %s" % TABLE_TEX]
    mine = {}
    for Re in RES:
        for tr in TRANS:
            for alpha in ALPHAS:
                for r in series(recs, Re, tr, alpha):
                    mine[(Re, r["N"], tr, alpha)] = (r["I_D"], r["I_P"],
                                                     r["delta"], r["diff"])
    bad = []
    for k, (id_, ip, dl, df) in sorted(frozen.items()):
        m = mine.get(k)
        if m is None:
            bad.append("absent du recalcul : Re=%d N=%d %s alpha=%.2f"
                       % (k[0], k[1], k[2], k[3]))
            continue
        for nm, x, y in (("I_D", id_, m[0]), ("I_P", ip, m[1]),
                         ("diff", df, m[3])):
            if abs(x - y) > 1.5e-5:
                bad.append("%s Re=%d N=%d %s a=%.2f : gel=%.5f recalc=%.5f"
                           % (nm, k[0], k[1], k[2], k[3], x, y))
        if dl is None and abs(m[0]) >= 1e-4:
            bad.append("delta gele '--' mais |I_D|=%.5f recalc=%.2f  %s"
                       % (m[0], m[2], k))
        elif dl is not None and m[2] != m[2]:
            bad.append("delta gele=%.2f mais |I_D|=%.5f < 1e-4 (sans sens)  %s"
                       % (dl, m[0], k))
        elif dl is not None and abs(dl - m[2]) > 0.02:
            bad.append("delta Re=%d N=%d %s a=%.2f : gel=%.2f recalc=%.2f"
                       % (k[0], k[1], k[2], k[3], dl, m[2]))
    extra = [k for k in mine if k not in frozen]
    return bad + ["ligne recalculee hors table : Re=%d N=%d %s a=%.2f" % k
                  for k in sorted(extra)]


def coverage():
    """Etat des cibles de campagnes + derniere progression vue dans les logs."""
    recs = load_converged()
    lines, todo = [], []
    for camp, cases in TARGETS.items():
        lines.append("== %s ==" % camp)
        for c in cases:
            cid = case_id(**c)
            key = (c["n"], round(c["alpha"], 4), c["Nw"], c["Re"], c["N"])
            if key in recs and ok(recs[key]):
                r = recs[key]
                lines.append("  [OK]      %-42s D*=%.5f P*=%.5f it=%s res=%.1e"
                             % (cid, r["D_star"], r["P_star"], r.get("iters"),
                                r.get("residual") or 0))
                continue
            log = os.path.join(LOGS, cid + ".log")
            prog = "aucun log"
            if os.path.exists(log):
                t = io.open(log, encoding="utf-8", errors="replace").read()
                m = re.findall(r"it=\s*(\d+)\s+dU=([\d.eE+-]+)", t)
                prog = ("it=%s dU=%s" % m[-1]) if m else "log sans progression"
            lines.append("  [EN ATT]  %-42s %s" % (cid, prog))
            todo.append(cid)
        lines.append("")
    return lines, todo, recs


def write_csv(recs):
    with open(os.path.join(OUT, "converged_runs.csv"), "w", newline="",
              encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["n", "alpha", "Nw", "Re", "N", "D_star", "P_star", "K_star",
                    "Z", "fluid_area", "Ma", "iters", "residual", "balance",
                    "n_vortices", "case_id"])
        for key in sorted(recs):
            r = recs[key]
            w.writerow([r["n"], r["alpha"], r.get("Nw"), r["Re"], r["N"],
                        r["D_star"], r["P_star"], r.get("K_star"), r.get("Z"),
                        r["fluid_area"], r.get("Ma"), r.get("iters"),
                        r.get("residual"), r.get("balance"),
                        r.get("n_vortices"), os.path.basename(r["path"])[:-5]])


def write_contrasts(recs):
    with open(os.path.join(OUT, "contrasts.csv"), "w", newline="",
              encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["Re", "transition", "alpha", "N", "I_D", "I_P", "delta",
                    "IP_minus_ID_pp", "I_area", "I_per_area",
                    "area_over_measured", "span_pp", "Ma_ref", "Ma_cmp"])
        for Re in RES:
            for tr in TRANS:
                for alpha in ALPHAS + (0.05, 0.20):
                    rows = series(recs, Re, tr, alpha)
                    if not rows:
                        continue
                    span = max(r["I_D"] for r in rows) - min(r["I_D"] for r in rows)
                    for r in rows:
                        w.writerow([Re, "(%d,%d)" % tr, alpha, r["N"],
                                    "%.6f" % r["I_D"], "%.6f" % r["I_P"],
                                    "%.3f" % r["delta"], "%.6f" % r["diff"],
                                    "%.6f" % r["I_area"], "%.6f" % r["I_per_area"],
                                    "%.1f" % r["ratio"], "%.6f" % span,
                                    r["Ma_a"], r["Ma_b"]])


def main():
    os.makedirs(OUT, exist_ok=True)
    report_only = "--report" in sys.argv
    lines, todo, recs = coverage()
    print("\n".join(lines))
    print("cas converges exploitables : %d" % sum(1 for r in recs.values() if ok(r)))
    if report_only:
        return 0

    bad = validate(recs)
    print("=== VALIDATION vs table SM.X gelée ===")
    if bad:
        for b in bad:
            print("  ECART  " + b)
    else:
        print("  OK : toutes les lignes gelees sont reproduites a 1.5e-5 pres")

    rows = tex_rows(recs)
    p = os.path.join(OUT, "SMX_rows_auto.tex")
    with io.open(p, "w", encoding="utf-8") as fh:
        fh.write("%% genere par solver/extract_converged.py -- ne pas editer\n")
        fh.write("%% lignes recalculees a partir des JSON converges (N jusqu'a %d)\n"
                 % max(N_LEVELS))
        fh.write("$Re$ & $N$ & Transition & $\\alpha$ & $I_D$ [\\%] & $I_P$ [\\%]"
                 " & $\\delta_{IP/ID}$ [\\%] & $I_P-I_D$ [pp] & span [pp] \\\\\n")
        fh.write("\n".join(rows) + "\n")
    print("\n=== sorties ===")
    write_csv(recs)
    write_contrasts(recs)
    for f in ("converged_runs.csv", "contrasts.csv", "SMX_rows_auto.tex"):
        print("  results/extract/%s (%d octets)"
              % (f, os.path.getsize(os.path.join(OUT, f))))
    if todo:
        print("\n%d cible(s) de campagne encore attendues ; relancer ce script"
              " une fois convergées." % len(todo))
        print("  " + "\n  ".join(todo))
    return 0


if __name__ == "__main__":
    sys.exit(main())
