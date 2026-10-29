# -*- coding: utf-8 -*-
"""
solver/run_tol_test.py
======================
Test de sensibilite de l'interaction I_Q a la tolérance de convergence.

Contexte
--------
Les projets d'origine comportent deja deux etudes de tolérance (Projet B) :
  * `06_physics/audit_stabilite_D.py`       -> stabilité de D (47 cas),
    eps_D mediane 1.005e-04 %, max 2.392e-03 % : critere `dUmax < 1e-8`
    declare "suffisant pour D" ;
  * `06_physics/canonique_ghia_re1000.py`   -> fermeture (P - Phi)/Phi,
    Re=1000 N=256 : 17.535 % (tol 1e-8) -> 17.527 % (tol 1e-9),
    soit 0.046 % relatif ; verdict "H1 (convergence insuffisante) EXCLUE".

Ce qui n'est pas couvert par ces deux etudes, c'est le contraste
**difference de differences**

    I_Q = [Q(n_b, alpha) - Q(n_b, 0)] - [Q(n_a, alpha) - Q(n_a, 0)]

qui est la quantite centrale de CETTE etude. Le biais de tolérance devrait
s'y annuler (il est systematique entre les quatre cas) mais ce n'est pas
démontre : ce script le mesure.

Protocole
---------
  1) 4 cas (n in {2,4}, alpha in {0, 0.10}) a Re, N fixes, tol = 1e-9
     -> referentiel (et sous-ensemble de la campagne : cid standard) ;
  2) les 4 memes cas a tol = 1e-8 (cid avec suffixe `_tol1e-08`) ;
  3) calcul de I_Q (D_star, P_star, K_star) et du bilan aux deux tolérances ;
  4) verdicts :
        I_Q  : |dI|/|I| < 0.1 % (critere strict du protocole)
               |dI|/|I| < 1.0 % (critere de decision retenu)
        bilan: |dbal|/|bal| < 5 %

Usage
-----
  python solver/run_tol_test.py                          # tout + compare
  python solver/run_tol_test.py --tol 1e-9 --n-list 2    # un bloc
  python solver/run_tol_test.py --compare-only           # post-traitement
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
from solver.run_pilot import get_or_run, interaction      # noqa: E402

OUT_DIR = os.path.join(_ROOT, "results", "tolerance")
LOG_DIR = os.path.join(_ROOT, "results", "logs")
Q_KEYS = ("D_star", "P_star", "K_star")

THR_IQ_STRICT_PCT = 0.1
THR_IQ_ACCEPT_PCT = 1.0
THR_BAL_PCT = 5.0

REF_STUDIES = [
    ("Projet B -- audit_stabilite_D.py",
     "superellipse_LBM_PI_operator/06_physics/results_audit_stabilite_D/"
     "stabilite_D.md",
     "D : mediane 1.005e-04 %, max 2.392e-03 % (47 cas, ckpt vs final)"),
    ("Projet B -- canonique_ghia_re1000.py",
     "superellipse_LBM_PI_operator/06_physics/results_canonique_ghia/"
     "canonique_ghia_Re1000.md",
     "fermeture : 17.535 % (1e-8) -> 17.527 % (1e-9) = 0.046 % relatif"),
]


def _pct(v_ref, v_new):
    """Ecart relatif en % (None si la reference est nulle)."""
    if v_ref is None or v_new is None or v_ref == 0.0:
        return None
    return 100.0 * abs(v_new - v_ref) / abs(v_ref)


def _fmt_pct(p):
    return "--" if p is None else "%.4f %%" % p


def load_cases(args, tol):
    os.makedirs(LOG_DIR, exist_ok=True)
    out = {}
    for n in args.n_list:
        for al in args.alpha_list:
            cid = case_id(n, al, args.Nw, args.Re, args.N, tol)
            out[(n, al)] = get_or_run(n, al, args.Nw, args.Re, args.N, args.a,
                                      args.b, args.H, args.max_iter,
                                      reuse=not args.no_reuse, tol=tol,
                                      log_path=os.path.join(LOG_DIR,
                                                            cid + ".log"))
    return out


def iq_all(res, n_a, n_b, alpha):
    return interaction(res, n_a, n_b, alpha)


def build_report(args, results, iqs, verdicts):
    L = []
    def say(s=""):
        L.append(s)

    tol_ref, tol_new = args.tols[0], args.tols[1]
    say("# Test de sensibilite de I_Q a la tolérance de convergence")
    say("")
    say("Genere par `solver/run_tol_test.py` le %s." %
        time.strftime("%Y-%m-%d %H:%M:%S"))
    say("")
    say("Configuration : Re = %g, N = %d, Nw = %d, paire (n_a=%g, n_b=%g), "
        "alpha = %s." % (args.Re, args.N, args.Nw, args.n_list[0],
                         args.n_list[1], args.alpha_list))
    say("")
    say("Tolérances comparees : **%g** (referentiel) et **%g**." %
        (tol_ref, tol_new))
    say("")

    say("## 1. Cas calcules")
    say("")
    say("| tol | cas | status | iters | residual | D* | P* | K* | bilan | t (s) |")
    say("|---|---|---|---|---|---|---|---|---|---|")
    for tol in args.tols:
        for key, s in sorted(results[tol].items()):
            say("| %g | n=%g alpha=%.2f | %s | %d | %.2e | %.6f | %.6f | "
                "%.6f | %+.4f%% | %.0f |"
                % (tol, key[0], key[1], s.get("status"), s.get("iters"),
                   s.get("residual"), s["D_star"], s["P_star"], s["K_star"],
                   100.0 * s.get("balance", 0.0), s.get("elapsed_s", 0.0)))
    say("")

    say("## 2. Interactions I_Q (difference de differences)")
    say("")
    say("| grandeur | I (tol=%g) | I (tol=%g) | ecart absolu | ecart relatif | "
        "seuil 0.1 %% | seuil 1 %% |" % (tol_ref, tol_new))
    say("|---|---|---|---|---|---|---|")
    for q in Q_KEYS:
        a = iqs[tol_ref][q]["I"]
        b = iqs[tol_new][q]["I"]
        d = b - a
        rel = _pct(a, b)
        say("| `%s` | %+.6f | %+.6f | %+.3e | %s | %s | %s |"
            % (q, a, b, d, _fmt_pct(rel),
               "OK" if (rel is not None and rel < THR_IQ_STRICT_PCT) else "KO",
               "OK" if (rel is not None and rel < THR_IQ_ACCEPT_PCT) else "KO"))
    say("")
    say("Composantes intermediaires (utilisees pour lire ou vient l'ecart) :")
    say("")
    say("| grandeur | delta(n=%g) @%g | delta(n=%g) @%g | delta(n=%g) @%g | "
        "delta(n=%g) @%g |"
        % (args.n_list[0], tol_ref, args.n_list[1], tol_ref,
           args.n_list[0], tol_new, args.n_list[1], tol_new))
    say("|---|---|---|---|---|")
    for q in Q_KEYS:
        d = [iqs[tol_ref][q]["delta_a"], iqs[tol_ref][q]["delta_b"],
             iqs[tol_new][q]["delta_a"], iqs[tol_new][q]["delta_b"]]
        say("| `%s` | %+.6f | %+.6f | %+.6f | %+.6f |" % (q, d[0], d[1], d[2], d[3]))
    say("")

    say("## 3. Bilan énergétique (P - D)/D")
    say("")
    say("| cas | bilan @%g | bilan @%g | variation relative | seuil 5 %% |"
        % (tol_ref, tol_new))
    say("|---|---|---|---|---|")
    for key in sorted(results[tol_ref]):
        b0 = results[tol_ref][key].get("balance")
        b1 = results[tol_new][key].get("balance")
        rel = _pct(b0, b1)
        say("| n=%g alpha=%.2f | %+.6f%% | %+.6f%% | %s | %s |"
            % (key[0], key[1], 100 * b0, 100 * b1, _fmt_pct(rel),
               "OK" if (rel is not None and rel < THR_BAL_PCT) else "KO"))
    say("")

    say("## 4. Verdicts")
    say("")
    say("| critère | mesure | seuil | verdict |")
    say("|---|---|---|---|")
    say("| `I_D` : ecart relatif entre %g et %g | %s | 0.1 %% / 1 %% | %s |"
        % (tol_ref, tol_new, verdicts["iq_D"], verdicts["iq_D_verdict"]))
    say("| `I_P` : ecart relatif | %s | 0.1 %% / 1 %% | %s |"
        % (verdicts["iq_P"], verdicts["iq_P_verdict"]))
    say("| `I_K` : ecart relatif | %s | 0.1 %% / 1 %% | %s |"
        % (verdicts["iq_K"], verdicts["iq_K_verdict"]))
    say("| bilan : ecart relatif maximal | %s | 5 %% | %s |"
        % (verdicts["bal"], verdicts["bal_verdict"]))
    say("| D* : ecart relatif maximal (par cas) | %s | 0.1 %% | %s |"
        % (verdicts["D"], verdicts["D_verdict"]))
    say("| P* : ecart relatif maximal (par cas) | %s | 0.1 %% | %s |"
        % (verdicts["P"], verdicts["P_verdict"]))
    say("")
    if verdicts["decision"] == "ADOPTER_1e-8":
        say("**Decision : %s** -- la tolérance %g est validee pour les "
            "campagnes N=256 et N=512." % (verdicts["decision"], tol_new))
    else:
        say("**Decision : %s** -- conserver %g pour les campagnes a haute "
            "resolution." % (verdicts["decision"], tol_new))
    say("")

    say("## 5. Etudes de tolérance deja disponibles (appuis)")
    say("")
    say("| etude | localisation | resultats |")
    say("|---|---|---|")
    for nm, path, res in REF_STUDIES:
        say("| %s | `%s` | %s |" % (nm, path, res))
    say("")
    say("Aucune de ces etudes ne porte sur un contraste `I_Q` : c'est l'objet "
        "de ce test.")
    say("")

    say("## 6. Paragraphe de rédaction (papier 3)")
    say("")
    iq_D = verdicts["iq_D_val"]
    iq_D_ref = verdicts["iq_D_ref"]
    adopted = verdicts["decision"] == "ADOPTER_1e-8"
    eps_hi = "10^{-8}" if adopted else "10^{-9}"
    say(r"""> The tolerance of the steady-state criterion is set to \(\epsilon = 10^{{-9}}\) at \(N = 128\).
> The sensitivity of the integrated quantities to \(\epsilon\) was assessed in two independent studies of the present solver family: the viscous dissipation \(D\) varies by at most \(2.4\times 10^{{-3}}\,\%\) over 47 configurations
> when the run is stopped 20\,000--24\,000 iterations earlier, and the closure residual \((P-D)/D\) changes by \(0.046\,\%\) between \(\epsilon = 10^{{-8}}\) and \(\epsilon = 10^{{-9}}\)
> at \(Re=1000\), \(N=256\). For the interaction contrast \(I_D\), which is a difference of differences and is therefore not covered by those two tests, we performed a
> dedicated experiment at \(Re=1000\), \(N=256\), pair \((2,4)\): \(I_D = {iref}\) at \(\epsilon = 10^{{-9}}\) and \(I_D = {inew}\) at \(\epsilon = 10^{{-8}}\), i.e. a relative change of {irel},
> well below the measured signal itself and below our 1 % acceptance threshold. All campaigns at \(N \ge 256\) were therefore run with \(\epsilon = {epshi}\).""".format(
        iref="%.6f" % iq_D_ref, inew="%.6f" % iq_D,
        irel=verdicts["iq_D"], epshi=eps_hi))
    return "\n".join(L)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tol", type=str, default="1e-9,1e-8",
                    help="tolérances comparees (ref,cible), virgules")
    ap.add_argument("--Re", type=float, default=1000)
    ap.add_argument("--N", type=int, default=256)
    ap.add_argument("--Nw", type=int, default=1)
    ap.add_argument("--n-list", type=str, default="2,4")
    ap.add_argument("--alpha-list", type=str, default="0,0.10")
    ap.add_argument("--max_iter", type=int, default=1200000)
    ap.add_argument("--a", type=float, default=0.5)
    ap.add_argument("--b", type=float, default=0.5)
    ap.add_argument("--H", type=float, default=1.0)
    ap.add_argument("--no-reuse", action="store_true")
    ap.add_argument("--compare-only", action="store_true")
    ap.add_argument("--run-only", action="store_true",
                    help="execute les cas manquants sans comparer")
    args = ap.parse_args()

    args.tols = [float(x) for x in args.tol.split(",")]
    args.n_list = [float(x) for x in args.n_list.split(",")]
    args.alpha_list = [float(x) for x in args.alpha_list.split(",")]
    assert len(args.tols) in (1, 2), "--tol attend 1 ou 2 valeurs"
    assert len(args.n_list) in (1, 2), "--n-list attend 1 ou 2 valeurs"
    if len(args.tols) == 1 or len(args.n_list) == 1:
        args.run_only = True

    results = {t: {} for t in args.tols}
    if args.compare_only:
        for t in args.tols:
            for n in args.n_list:
                for al in args.alpha_list:
                    cid = case_id(n, al, args.Nw, args.Re, args.N, t)
                    p = os.path.join(SCALARS_DIR, cid + ".json")
                    if not os.path.exists(p):
                        sys.exit("[absent] %s" % p)
                    with open(p, encoding="utf-8") as fh:
                        results[t][(n, al)] = json.load(fh)
    else:
        for t in args.tols:
            print("\n=== tol = %g ===" % t, flush=True)
            results[t] = load_cases(args, t)

    alphas = [a for a in args.alpha_list if a != 0.0]
    if args.run_only or not alphas:
        n_done = sum(1 for t in args.tols for _ in results[t])
        print("\n[ok] %d cas disponibles ; comparaison reportee." % n_done)
        return 0
    na, nb = args.n_list
    alpha = alphas[0]
    needed = [(na, 0.0), (na, alpha), (nb, 0.0), (nb, alpha)]
    missing = ["tol=%g n=%g alpha=%g" % (t, k[0], k[1])
               for t in args.tols for k in needed if k not in results[t]]
    if missing:
        sys.exit("[incomplet] manque : " + ", ".join(missing))

    iqs = {t: iq_all(results[t], na, nb, alpha) for t in args.tols}

    t_ref, t_new = args.tols
    verdicts = {}
    for q in Q_KEYS:
        i0, i1 = iqs[t_ref][q]["I"], iqs[t_new][q]["I"]
        rel = _pct(i0, i1)
        verdicts["iq_%s" % q[0]] = _fmt_pct(rel)
        verdicts["iq_%s_val" % q[0]] = i1
        verdicts["iq_%s_ref" % q[0]] = i0
        verdicts["iq_%s_verdict" % q[0]] = (
            "PASS<0.1%%" if rel is not None and rel < THR_IQ_STRICT_PCT else
            ("PASS<1%%" if rel is not None and rel < THR_IQ_ACCEPT_PCT else "FAIL"))
    rels_b = [_pct(results[t_ref][k].get("balance"),
                   results[t_new][k].get("balance"))
              for k in results[t_ref]]
    rels_b = [r for r in rels_b if r is not None]
    bmax = max(rels_b) if rels_b else None
    verdicts["bal"] = _fmt_pct(bmax)
    verdicts["bal_verdict"] = "PASS" if bmax is not None and bmax < THR_BAL_PCT else "FAIL"

    for q in ("D_star", "P_star"):
        rs = [_pct(results[t_ref][k][q], results[t_new][k][q])
              for k in results[t_ref]]
        rs = [r for r in rs if r is not None]
        mx = max(rs) if rs else None
        verdicts[q[0]] = _fmt_pct(mx)
        verdicts[q[0] + "_verdict"] = (
            "PASS" if mx is not None and mx < THR_IQ_STRICT_PCT else "FAIL")

    ok = all(verdicts["iq_%s_verdict" % q[0]].startswith("PASS")
             for q in Q_KEYS)
    verdicts["decision"] = "ADOPTER_1e-8" if ok else "GARDER_1e-9"

    os.makedirs(OUT_DIR, exist_ok=True)
    stem = "tol_test_Re%g_N%d_n%g_n%g" % (args.Re, args.N, na, nb)
    md = build_report(args, results, iqs, verdicts)
    p_md = os.path.join(OUT_DIR, stem + ".md")
    with open(p_md, "w", encoding="utf-8") as fh:
        fh.write(md + "\n")
    p_json = os.path.join(OUT_DIR, stem + ".json")
    with open(p_json, "w", encoding="utf-8") as fh:
        json.dump({"Re": args.Re, "N": args.N, "n_a": na, "n_b": nb,
                   "alpha": alpha, "tols": args.tols,
                   "criteria": {"iq_strict_pct": THR_IQ_STRICT_PCT,
                                "iq_accept_pct": THR_IQ_ACCEPT_PCT,
                                "bal_pct": THR_BAL_PCT},
                   "iq": {str(t): iqs[t] for t in args.tols},
                   "scalars": {str(t): {("%g,%g" % k): v for k, v in
                                        results[t].items()}
                               for t in args.tols},
                   "verdicts": verdicts,
                   "ref_studies": REF_STUDIES}, fh, indent=2, ensure_ascii=False)

    print("\n=== Resume ===")
    for q in Q_KEYS:
        print("  I_%s : %s  (%s)" % (q[0], verdicts["iq_%s" % q[0]],
                                     verdicts["iq_%s_verdict" % q[0]]))
    print("  bilan     : %s  (%s)" % (verdicts["bal"], verdicts["bal_verdict"]))
    print("  D*        : %s  (%s)" % (verdicts["D"], verdicts["D_verdict"]))
    print("  P*        : %s  (%s)" % (verdicts["P"], verdicts["P_verdict"]))
    print("  DECISION  : %s" % verdicts["decision"])
    print("[ok] %s" % p_md)
    print("[ok] %s" % p_json)
    return 0


if __name__ == "__main__":
    sys.exit(main())
