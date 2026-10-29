# -*- coding: utf-8 -*-
"""
solver/run_res_test.py
======================
Test decisif : stabilité de l'interaction I_Q entre deux resolutions N.

Contexte
--------
A N=128, l'interaction I_D vaut -1.42590e-01 a Re=1000 (paire 2,4), soit
-0.470 % de D(n=2, alpha=0). Entre N=128 et N=256, D* bouge de +2.99 % et
P* de +6.69 % : le rapport signal/resolution est donc de 5 a 6.

Le contraste I_Q est une difference de differences : les erreurs de
resolution communes aux quatre cas devraient s'annuler, mais cette
annulation n'est pas garantie. Ce script la mesure.

Protocole
---------
  1) 4 cas (n in {2,4}, alpha in {0, 0.10}) a Re, N_lo fixes ;
  2) les 4 memes cas a N_hi ;
  3) I_Q (D_star, P_star, K_star) aux deux resolutions ;
  4) ecart relatif |I(N_hi) - I(N_lo)| / |I(N_lo)| et verdicts :
        < 10 %   -> EXCELLENTE stabilite (papier publiable en l'etat)
        10-30 %  -> stabilite acceptable (discuter l'incertitude)
        > 30 %   -> probleme, comprendre avant de continuer
  5) rapport signal/resolution : |I| / |dD*|_resolution moyen.

Usage
-----
  python solver/run_res_test.py                    # compare (cas deja presents)
  python solver/run_res_test.py --run              # lance les cas manquants
  python solver/run_res_test.py --Re 400           # un autre Reynolds
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
_PROV_A = os.path.join(_ROOT, "provenance", "Lid_driven_ondulation")
for _p in (_HERE, _ROOT, _PROV_A):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from solver.run_case import case_id, SCALARS_DIR          # noqa: E402
from solver.run_pilot import interaction                  # noqa: E402

OUT_DIR = os.path.join(_ROOT, "results", "resolution")
LOG_DIR = os.path.join(_ROOT, "results", "logs")
Q_KEYS = ("D_star", "P_star", "K_star")

THR_EXCELLENT_PCT = 10.0
THR_ACCEPTABLE_PCT = 30.0


def _pct(v_ref, v_new):
    if v_ref is None or v_new is None or v_ref == 0.0:
        return None
    return 100.0 * abs(v_new - v_ref) / abs(v_ref)


def _fmt_pct(p):
    return "--" if p is None else "%.3f %%" % p


def _verdict(rel):
    if rel is None:
        return "n/a"
    if rel < THR_EXCELLENT_PCT:
        return "EXCELLENT (<10%)"
    if rel < THR_ACCEPTABLE_PCT:
        return "ACCEPTABLE (10-30%)"
    return "PROBLEME (>30%)"


def load_block(args, N):
    """Charge (ou lance si --run) les 4 cas de la paire a une resolution."""
    out = {}
    for n in args.n_list:
        for al in args.alpha_list:
            cid = case_id(n, al, args.Nw, args.Re, N, args.tol)
            path = os.path.join(SCALARS_DIR, cid + ".json")
            if not os.path.exists(path):
                if not args.run:
                    sys.exit("[absent] %s  (relancer avec --run)" % path)
                from solver.run_pilot import get_or_run
                os.makedirs(LOG_DIR, exist_ok=True)
                out[(n, al)] = get_or_run(
                    n, al, args.Nw, args.Re, N, args.a, args.b, args.H,
                    args.max_iter, reuse=not args.no_reuse, tol=args.tol,
                    log_path=os.path.join(LOG_DIR, cid + ".log"))
                continue
            with open(path, encoding="utf-8") as fh:
                rec = json.load(fh)
            if rec.get("status") != "converged":
                sys.exit("[non converge] %s status=%s"
                         % (cid, rec.get("status")))
            if float(rec.get("tol", args.tol)) != float(args.tol):
                sys.exit("[tol differente] %s tol=%s attendu %g"
                         % (cid, rec.get("tol"), args.tol))
            print("[reuse] %s" % cid, flush=True)
            out[(n, al)] = rec
    return out


def build_report(args, res, iqs):
    L = []
    def say(s=""):
        L.append(s)

    n_lo, n_hi = args.N_list
    na, nb = args.n_list
    alpha = [a for a in args.alpha_list if a != 0.0][0]

    say("# Test de stabilite de I_Q en resolution (N=%d -> N=%d)"
        % (n_lo, n_hi))
    say("")
    say("Genere par `solver/run_res_test.py` le %s."
        % time.strftime("%Y-%m-%d %H:%M:%S"))
    say("")
    say("Configuration : Re = %g, Nw = %d, paire (n_a=%g, n_b=%g), "
        "alpha = %g, tol = %g." % (args.Re, args.Nw, na, nb, alpha,
                                   args.tol))
    say("")
    say("Question : l'erreur de resolution, qui deplace D* de %+.2f %% "
        "entre %d et %d pour le cas (n=2, alpha=0), s'annule-t-elle dans "
        "le contraste I_Q ?" % (
            100.0 * (res[n_hi][(na, 0.0)]["D_star"]
                     - res[n_lo][(na, 0.0)]["D_star"])
            / res[n_lo][(na, 0.0)]["D_star"], n_lo, n_hi))
    say("")

    say("## 1. Cas calcules")
    say("")
    say("| N | cas | status | iters | residual | D* | P* | K* | bilan |")
    say("|---|---|---|---|---|---|---|---|---|")
    for N in (n_lo, n_hi):
        for key in sorted(res[N]):
            s = res[N][key]
            say("| %d | n=%g alpha=%.2f | %s | %d | %.2e | %.6f | %.6f | "
                "%.6f | %+.4f%% |"
                % (N, key[0], key[1], s.get("status"), s.get("iters"),
                   s.get("residual"), s["D_star"], s["P_star"], s["K_star"],
                   100.0 * s.get("balance", 0.0)))
    say("")

    say("## 2. Sensibilite absolue a la resolution")
    say("")
    say("| cas | dD*/D* | dP*/P* | dK*/K* |")
    say("|---|---|---|---|")
    for key in sorted(res[n_lo]):
        row = "| n=%g alpha=%.2f |" % (key[0], key[1])
        for q in Q_KEYS:
            r = _pct(res[n_lo][key][q], res[n_hi][key][q])
            row += " %s |" % _fmt_pct(r)
        say(row)
    say("")

    say("## 3. Interactions I_Q aux deux resolutions")
    say("")
    say("| grandeur | I (N=%d) | I (N=%d) | ecart absolu | ecart relatif | "
        "verdict |" % (n_lo, n_hi))
    say("|---|---|---|---|---|---|")
    verdicts = {}
    for q in Q_KEYS:
        a = iqs[n_lo][q]["I"]
        b = iqs[n_hi][q]["I"]
        rel = _pct(a, b)
        verdicts[q] = rel
        say("| `%s` | %+.6f | %+.6f | %+.3e | %s | %s |"
            % (q, a, b, b - a, _fmt_pct(rel), _verdict(rel)))
    say("")

    say("## 4. Composantes intermediaires")
    say("")
    say("| grandeur | delta(n=%g) @N=%d | delta(n=%g) @N=%d | "
        "delta(n=%g) @N=%d | delta(n=%g) @N=%d |"
        % (na, n_lo, nb, n_lo, na, n_hi, nb, n_hi))
    say("|---|---|---|---|---|")
    for q in Q_KEYS:
        d = [iqs[n_lo][q]["delta_a"], iqs[n_lo][q]["delta_b"],
             iqs[n_hi][q]["delta_a"], iqs[n_hi][q]["delta_b"]]
        say("| `%s` | %+.6f | %+.6f | %+.6f | %+.6f |" % (q, d[0], d[1],
                                                          d[2], d[3]))
    say("")

    say("## 5. Rapport signal / resolution")
    say("")
    say("D* bouge de %s entre les deux resolutions (cas n=%g, alpha=0). "
        "Le signal d'interaction vaut %s a N=%d."
        % (_fmt_pct(_pct(res[n_lo][(na, 0.0)]["D_star"],
                         res[n_hi][(na, 0.0)]["D_star"])),
           na, "%+.6f" % iqs[n_hi]["D_star"]["I"], n_hi))
    say("")
    say("| grandeur | |I(N=%d)| / |dQ*_res| | lecture |"
        % n_hi)
    say("|---|---|---|")
    for q in Q_KEYS:
        d_res = abs(res[n_hi][(na, 0.0)][q] - res[n_lo][(na, 0.0)][q])
        i_abs = abs(iqs[n_hi][q]["I"])
        ratio = i_abs / d_res if d_res else None
        if ratio is None:
            lect = "--"
        elif ratio > 0.2:
            lect = "signal clairement discernable"
        elif ratio > 0.05:
            lect = "zone limite : annulation partielle attendue"
        else:
            lect = "signal petit devant le bruit de resolution"
        say("| `%s` | %s | %s |" % (q,
                                    "--" if ratio is None else "%.2f" % ratio,
                                    lect))
    say("")

    say("## 6. Verdict global")
    say("")
    rel_D = verdicts["D_star"]
    if rel_D is not None and rel_D < THR_EXCELLENT_PCT:
        decision = "STABLE"
        msg = ("I_D est stable a %s entre N=%d et N=%d : le signal "
               "d'interaction est robuste face a la resolution. Les "
               "campagnes N=256 puis N=512 peuvent etre lancees telles "
               "qu'elles sont planifiees."
               % (_fmt_pct(rel_D), n_lo, n_hi))
    elif rel_D is not None and rel_D < THR_ACCEPTABLE_PCT:
        decision = "STABLE_AVEC_INCERTITUDE"
        msg = ("I_D varie de %s entre les deux resolutions : stabilite "
               "acceptable. Discuter honnetement de l'incertitude dans le "
               "papier et preferer N=512 pour le cas de reference."
               % _fmt_pct(rel_D))
    else:
        decision = "INSTABLE"
        msg = ("I_D varie de %s entre les deux resolutions : au-dela du "
               "seuil de 30 %%. Comprendre l'origine de l'instabilite "
               "(bilan energetique, comptage de cellules fluides, effet "
               "de paroi) avant de poursuivre." % _fmt_pct(rel_D))
    say("| critere | mesure | seuil | verdict |")
    say("|---|---|---|---|")
    say("| `I_D` : ecart relatif N=%d -> N=%d | %s | 10 %% / 30 %% | %s |"
        % (n_lo, n_hi, _fmt_pct(rel_D), _verdict(rel_D)))
    say("| `I_P` : ecart relatif | %s | 10 %% / 30 %% | %s |"
        % (_fmt_pct(verdicts["P_star"]), _verdict(verdicts["P_star"])))
    say("| `I_K` : ecart relatif | %s | 10 %% / 30 %% | %s |"
        % (_fmt_pct(verdicts["K_star"]), _verdict(verdicts["K_star"])))
    say("")
    say("**Decision : %s.** %s" % (decision, msg))
    say("")

    say("## 7. Paragraphe de redaction (papier 3)")
    say("")
    say(r"""> The interaction contrast must be shown to be free of resolution
> contamination. We therefore recomputed the full \((2,4)\) pair at
> \(Re = %g\) for \(N = %d\) and \(N = %d\). The dissipation of the reference
> configuration \((n=2,\alpha=0)\) shifts by %.2f \%% between the two grids,
> yet the contrast changes only from \(I_D = %+.5f\) at \(N=%d\) to
> \(I_D = %+.5f\) at \(N=%d\), i.e. %.2f \%% of relative change (%s).
> The matched-resolution cancellation is therefore effective to within that
> figure, which is used as the systematic uncertainty of the reported
> interaction.""" % (
        args.Re, n_lo, n_hi,
        100.0 * (res[n_hi][(na, 0.0)]["D_star"]
                 - res[n_lo][(na, 0.0)]["D_star"])
        / res[n_lo][(na, 0.0)]["D_star"],
        iqs[n_lo]["D_star"]["I"], n_lo,
        iqs[n_hi]["D_star"]["I"], n_hi,
        rel_D if rel_D is not None else float("nan"),
        _verdict(rel_D).lower()))
    return "\n".join(L), decision, verdicts


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--Re", type=float, default=1000)
    ap.add_argument("--N-list", type=str, default="128,256")
    ap.add_argument("--Nw", type=int, default=1)
    ap.add_argument("--n-list", type=str, default="2,4")
    ap.add_argument("--alpha-list", type=str, default="0,0.10")
    ap.add_argument("--tol", type=float, default=1e-9)
    ap.add_argument("--max_iter", type=int, default=1200000)
    ap.add_argument("--a", type=float, default=0.5)
    ap.add_argument("--b", type=float, default=0.5)
    ap.add_argument("--H", type=float, default=1.0)
    ap.add_argument("--no-reuse", action="store_true")
    ap.add_argument("--run", action="store_true",
                    help="lance les cas absents (checkpoint actif)")
    args = ap.parse_args()

    args.N_list = [int(x) for x in args.N_list.split(",")]
    args.n_list = [float(x) for x in args.n_list.split(",")]
    args.alpha_list = [float(x) for x in args.alpha_list.split(",")]
    assert len(args.N_list) == 2, "--N-list attend 2 valeurs"
    assert len(args.n_list) == 2, "--n-list attend 2 valeurs"

    res = {N: load_block(args, N) for N in args.N_list}
    na, nb = args.n_list
    alpha = [a for a in args.alpha_list if a != 0.0][0]
    iqs = {N: interaction(res[N], na, nb, alpha) for N in args.N_list}

    os.makedirs(OUT_DIR, exist_ok=True)
    stem = "res_test_Re%g_N%d_N%d_n%g_n%g" % (args.Re, args.N_list[0],
                                              args.N_list[1], na, nb)
    md, decision, verdicts = build_report(args, res, iqs)
    p_md = os.path.join(OUT_DIR, stem + ".md")
    with open(p_md, "w", encoding="utf-8") as fh:
        fh.write(md + "\n")
    p_json = os.path.join(OUT_DIR, stem + ".json")
    with open(p_json, "w", encoding="utf-8") as fh:
        json.dump({"Re": args.Re, "N_list": args.N_list, "n_a": na,
                   "n_b": nb, "alpha": alpha, "tol": args.tol,
                   "criteria": {"excellent_pct": THR_EXCELLENT_PCT,
                                "acceptable_pct": THR_ACCEPTABLE_PCT},
                   "iq": {str(N): iqs[N] for N in args.N_list},
                   "scalars": {str(N): {("%g,%g" % k): v
                                        for k, v in res[N].items()}
                               for N in args.N_list},
                   "verdicts": {k: v for k, v in verdicts.items()},
                   "decision": decision}, fh, indent=2, ensure_ascii=False)

    print("\n=== Resume ===")
    for q in Q_KEYS:
        print("  I_%s : %s  (%s)" % (q[0], _fmt_pct(verdicts[q]),
                                     _verdict(verdicts[q])))
    print("  DECISION : %s" % decision)
    print("[ok] %s" % p_md)
    print("[ok] %s" % p_json)
    return 0


if __name__ == "__main__":
    sys.exit(main())
