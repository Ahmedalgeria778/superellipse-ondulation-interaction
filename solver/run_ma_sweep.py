# -*- coding: utf-8 -*-
"""
solver/run_ma_sweep.py
======================
Priorite 2 (bilan energetique) : le defaut (P*-D*)/D* est-il une erreur de
discretisation (en h) ou un effet de compressibilite (en Ma^2) ?

L'affinement N=128 -> N=256 modifie a la fois h ET Ma (Ma = U*sqrt(3) avec
U = (40/3)/N : 0.1804 -> 0.0902). L'ordre apparent ~1 mesure donc les deux
a la fois. Ici on fixe N et Re et on fait varier U (donc nu et tau) :
le ecoulement adimensionne est le meme (Re fixe), seul Ma varie.

    nu = U*N/Re ,  tau = 0.5 + 3*nu ,  Ma = sqrt(3)*U

Le defaut est ajuste par  eps(Ma) = c0 + c1*Ma^2 :
    - c0 ~ 0  -> le defaut est une pure compressibilite ;
    - c0 = eps a Ma -> 0 visible -> le defaut est de discretisation.

Les cas sont ecrits dans results/ma_sweep/ (JAMAIS dans results/scalars,
pour ne ecraser aucun cas deja calcule).

Usage :
  python -m solver.run_ma_sweep --Re 400 --N 128 --U-list 0.026,0.052,0.104,0.208
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
_PROV_A = os.path.join(_ROOT, "provenance", "Lid_driven_ondulation")
for _p in (_HERE, _PROV_A):
    if _p not in sys.path:
        sys.path.insert(0, _p)

OUT_DIR = os.path.join(_ROOT, "results", "ma_sweep")


def one(args):
    n, alpha, Nw, Re, N, U, max_iter, tol = args
    from solver.run_case import run_one
    t0 = time.time()
    summ = run_one(n, alpha, Nw, Re, N, max_iter=max_iter, tol=tol,
                   verbose=False, save=False, ckpt=False, U_ref=U)
    nu = U * N / Re
    tau = 0.5 + 3.0 * nu
    return {
        "n": n, "alpha": alpha, "Nw": Nw, "Re": Re, "N": N,
        "U_ref": U, "nu": nu, "tau": tau, "Ma": math.sqrt(3.0) * U,
        "D_star": summ["D_star"], "P_star": summ["P_star"],
        "balance": summ["balance"], "K_star": summ["K_star"],
        "iters": summ.get("iters"), "residual": summ.get("residual"),
        "status": summ.get("status"), "elapsed_s": time.time() - t0,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--Re", type=float, default=400)
    ap.add_argument("--N", type=int, default=128)
    ap.add_argument("--n", type=float, default=2.0)
    ap.add_argument("--alpha", type=float, default=0.0)
    ap.add_argument("--Nw", type=int, default=1)
    ap.add_argument("--U-list", type=str, default="0.0260416667,0.0520833333,0.1041666667,0.2083333333")
    ap.add_argument("--max-iter", type=int, default=600000)
    ap.add_argument("--tol", type=float, default=1e-8)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--refit", action="store_true",
                    help="ne pas relancer : reprendre le JSON existant et "
                         "recalculer seul les ajustements")
    args = ap.parse_args()

    _tag = "" if (args.n == 2.0 and args.alpha == 0.0) else "_n%g_a%.2f" % (args.n, args.alpha)
    path = os.path.join(OUT_DIR, "masweep_Re%g_N%d%s.json" % (args.Re, args.N, _tag))
    if args.refit:
        if not os.path.exists(path):
            print("aucun JSON a refaire : %s" % path)
            return 1
        with open(path, "r", encoding="utf-8") as fh:
            old = json.load(fh)
        out = sorted(old["cases"], key=lambda r: r["Ma"])
        _report(out, args.N)
        old["fit"] = _fit(out)
        with open(path, "w", encoding="utf-8") as fh:
            json.dump(old, fh, indent=2, ensure_ascii=False)
        print("\n[ok] %s (ajustements seulement)" % path)
        return 0

    Us = [float(x) for x in args.U_list.split(",")]
    jobs = [(args.n, args.alpha, args.Nw, args.Re, args.N, U,
             args.max_iter, args.tol) for U in Us]

    os.makedirs(OUT_DIR, exist_ok=True)
    print("[ma-sweep] Re=%g N=%d n=%g alpha=%.2f | %d valeurs de U"
          % (args.Re, args.N, args.n, args.alpha, len(jobs)), flush=True)
    print("  U          nu         tau       Ma       bal(P-D)/D", flush=True)

    out = []
    with ProcessPoolExecutor(max_workers=min(args.workers, len(jobs))) as ex:
        futs = {ex.submit(one, j): j for j in jobs}
        for fu in as_completed(futs):
            try:
                r = fu.result()
            except Exception as exc:                     # noqa: BLE001
                print("  ECHEC %s: %s" % (futs[fu][5], exc), flush=True)
                continue
            out.append(r)
            print("  %-10.6f %-10.6f %-9.4f %-8.4f %+.4f %%   %s (%.0fs)"
                  % (r["U_ref"], r["nu"], r["tau"], r["Ma"],
                     r["balance"] * 100.0, r["status"], r["elapsed_s"]),
                  flush=True)

    out.sort(key=lambda r: r["Ma"])
    if not out:
        return 1

    _report(out, args.N)

    payload = {"Re": args.Re, "N": args.N, "n": args.n,
               "alpha": args.alpha, "cases": out, "fit": None}
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2, ensure_ascii=False)

    payload["fit"] = _fit(out)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2, ensure_ascii=False)
    print("\n[ok] %s" % path)
    return 0


def _linfit(x, y):
    """Regression lineaire y = c0 + c1*x, retourne (c0, c1, R^2)."""
    n = len(x)
    if n < 2:
        return float("nan"), float("nan"), float("nan")
    mx, my = sum(x) / n, sum(y) / n
    sxx = sum((a - mx) ** 2 for a in x)
    if sxx <= 0:
        return float("nan"), float("nan"), float("nan")
    c1 = sum((a - mx) * (b - my) for a, b in zip(x, y)) / sxx
    c0 = my - c1 * mx
    ss_tot = sum((b - my) ** 2 for b in y)
    ss_res = sum((b - (c0 + c1 * a)) ** 2 for a, b in zip(x, y))
    r2 = 1.0 - ss_res / ss_tot if ss_tot else float("nan")
    return c0, c1, r2


def _fit(out):
    """Ajuste le defaut eps(Ma) a h fixe : en Ma^2 et en Ma."""
    import math
    ma = [r["Ma"] for r in out]
    eps = [r["balance"] for r in out]
    if len(out) < 2:
        return None
    q = _linfit([a * a for a in ma], eps)
    l = _linfit(ma, eps)
    # point de reference : la loi de campagne U = (40/3)/N  ->  Ma fixe
    N = out[0]["N"]
    ma_camp = math.sqrt(3.0) * (40.0 / 3.0) / N
    ref = min(out, key=lambda r: abs(r["Ma"] - ma_camp))
    eps_ref = ref["balance"]
    spread = max(eps) - min(eps)
    print()
    print("=" * 76)
    print("Ajustement du defaut, h fixe (h = 1/%d)" % N)
    print("=" * 76)
    print("  reference : Ma_camp = %.4f  ->  eps = %+.4f %%" % (ma_camp, eps_ref * 100.0))
    print("  eps = c0 + c1*Ma^2 : c0 = %+7.4f %%  c1 = %+8.4f %%/Ma^2  R2 = %.4f"
          % (q[0] * 100.0, q[1] * 100.0, q[2]))
    print("  eps = c0 + c1*Ma   : c0 = %+7.4f %%  c1 = %+8.4f %%/Ma    R2 = %.4f"
          % (l[0] * 100.0, l[1] * 100.0, l[2]))
    print()
    print("  extrapolation Ma -> 0 / defaut de reference :")
    print("      en Ma^2 : %5.1f %%" % (q[0] / eps_ref * 100.0 if eps_ref else float("nan")))
    print("      en Ma   : %5.1f %%" % (l[0] / eps_ref * 100.0 if eps_ref else float("nan")))
    print("  variation de eps sur la plage balayee  : %+.3f pt%%  (%.0f %% du defaut)"
          % (spread * 100.0, spread / abs(eps_ref) * 100.0 if eps_ref else float("nan")))
    print()
    print("  c0 ~ eps -> le defaut est de discretisation (h) ;")
    print("  c0 ~ 0   -> il est de compressibilite (Ma).")
    return {"ma2": {"c0": q[0], "c1": q[1], "r2": q[2],
                    "frac_at_zero": q[0] / eps_ref if eps_ref else float("nan")},
            "ma": {"c0": l[0], "c1": l[1], "r2": l[2],
                   "frac_at_zero": l[0] / eps_ref if eps_ref else float("nan")},
            "ma_campaign": ma_camp, "eps_campaign": eps_ref,
            "spread_pt_pct": spread * 100.0,
            "eps_at_each_Ma": eps}


def _report(out, N):
    print()
    print("  U          nu         tau       Ma       bal(P-D)/D")
    for r in out:
        print("  %-10.6f %-10.6f %-9.4f %-8.4f %+.4f %%   %s (%.0fs)"
              % (r["U_ref"], r["nu"], r["tau"], r["Ma"],
                 r["balance"] * 100.0, r["status"], r["elapsed_s"]))
    print()
    print("  Decomposition P* et D* :")
    print("     Ma      D*         P*        eps")
    d_ref = out[0]["D_star"]
    p_ref = out[0]["P_star"]
    for r in out:
        print("   %.4f  %.5f   %.5f  %+.4f %%   dD %+.2f%%  dP %+.2f%%"
              % (r["Ma"], r["D_star"], r["P_star"], r["balance"] * 100.0,
                 (r["D_star"] - d_ref) / d_ref * 100.0,
                 (r["P_star"] - p_ref) / p_ref * 100.0))


if __name__ == "__main__":
    raise SystemExit(main())
