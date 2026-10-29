# -*- coding: utf-8 -*-
"""Bilan global des campagnes P1-P4 (run unique, lecture des JSON)."""
import json, math, os, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SC = os.path.join(ROOT, "results", "scalars")
FEM = os.path.join(ROOT, "results", "fem")
MA = os.path.join(ROOT, "results", "ma_sweep")

def load(N, Re, n, al, tol="_tol1e-08"):
    for s in (tol, ""):
        p = os.path.join(SC, "se_n%g_w%.2f_Nw1_Re%d_N%d%s.json" % (n, al, Re, N, s))
        if os.path.exists(p):
            try:
                with open(p, encoding="utf-8") as fh:
                    return json.load(fh)
            except Exception:
                return None
    return None

def ok(rec):
    return rec is not None and rec.get("status") in ("converged", "reanalyzed")

# ---------- 1) P1.1 : sequence complete n=2,4,6,8 ----------
print("=" * 78)
print("1) P1.1 — sequence n=2,4,6,8 a N=384 : I_D(na,nb)=(D*(nb,a)-D*(na,a))/D*(na,a)*100")
print("=" * 78)
have = {}
for Re in (100, 400, 1000):
    rows = []
    for n in (2, 4, 6, 8):
        for al in (0.0, 0.10):
            rec = load(384, Re, n, al)
            have[(Re, n, al)] = ok(rec)
            rows.append("%s" % ("O" if ok(rec) else "-"))
    print("Re=%-5g dispo n=2,2a,4,4a,6,6a,8,8a : %s" % (Re, " ".join(rows)))

def contrast(tab_key, Re, na, nb, key="D_star"):
    """Formule manuscrit (Table SM.X) : contraste relatif entre exposants.
    I_D(na,nb) = (Q(nb,a) - Q(na,a)) / Q(na,a) * 100, a l'alpha donne.
    (NE PAS confondre avec l'interaction appariee en alpha d'audit_contrasts.)
    """
    qn = {}
    for n in (na, nb):
        r = load(384, Re, n, tab_key)
        if not ok(r):
            return None
        qn[n] = r[key]
    return qn[na], qn[nb], 100.0 * (qn[nb] - qn[na]) / qn[na]

for Re in (100, 400, 1000):
    for al in (0.0, 0.10):
        for (na, nb) in ((2,4), (4,6), (6,8)):
            c = contrast(al, Re, na, nb)
            if c:
                qa, qb, I = c
                print("  Re=%-5g a=%.2f I_D(%d,%d) : D*(n%d)=%.4f D*(n%d)=%.4f  I=%+.4f%%"
                      % (Re, al, na, nb, na, qa, nb, qb, I))
            else:
                print("  Re=%-5g a=%.2f I_D(%d,%d) : incomplet" % (Re, al, na, nb))

# valeurs brutes des nouveaux cas
print("\n  scalars des cas P1.1 (D*, P*, K*, Z, it):")
for Re in (100, 1000):
    for n in (6, 8):
        for al in (0.0, 0.10):
            r = load(384, Re, n, al)
            if ok(r):
                print("   Re=%-5g n=%d a=%.2f : D*=%.4f P*=%.4f K*=%.5f Z=%.5f it=%s"
                      % (Re, n, al, r["D_star"], r["P_star"], r["K_star"], r["Z"],
                         r.get("iters")))

# ---------- 2) P1.2 : FEM vs LBM ----------
print("\n" + "=" * 78)
print("2) P1.2 — FEM (Taylor-Hood) vs LBM N=128, alpha=0")
print("=" * 78)
for Re in (100, 1000):
    for n in (2, 4):
        for nx in (64, 128, 256):
            p = os.path.join(FEM, "FEMX-n%d-a0.00-R%d-nx%d.json" % (n, Re, nx))
            if not os.path.exists(p):
                print("   Re=%-5g n=%d nx=%-3d : —" % (Re, n, nx))
                continue
            with open(p, encoding="utf-8") as fh:
                f = json.load(fh)
            lbm = load(128, Re, n, 0.0)
            line = "   Re=%-5g n=%d nx=%-3d : D*nat=%.4f P*nat=%.4f K*nat=%.5f it=%s" % (
                Re, n, nx, f.get("D_star_native", 0), f.get("P_star_native", 0),
                f.get("K_star_native", 0), f.get("iters"))
            if lbm and ok(lbm):
                ddn = 100.0 * abs(f["D_star_native"] - lbm["D_star"]) / lbm["D_star"]
                dkn = 100.0 * abs(f["K_star_native"] - lbm["K_star"]) / lbm["K_star"]
                line += "  | LBM: dD*=%.2f%% dK*=%.2f%%" % (ddn, dkn)
            print(line)

# ---------- 3) P2.1/P3.x : restants ----------
print("\n" + "=" * 78)
print("3) P2.1 (N=512, a=0) et P3.x — etat")
print("=" * 78)
for Re in (100, 400, 1000):
    for n in (2, 4):
        r = load(512, Re, n, 0.0)
        print("   N512 Re=%-5g n=%d : %s" % (Re, n, "CONVERGE" if ok(r) else "a lancer"))
for Re in (400, 1000):
    for n in (2, 4):
        p = os.path.join(SC, "se_n%g_w0.00_Nw2_Re%d_N384_tol1e-08.json" % (n, Re))
        st = "-"
        if os.path.exists(p):
            with open(p, encoding="utf-8") as fh:
                st = "CONVERGE" if json.load(fh).get("status") in ("converged","reanalyzed") else "?"
        print("   Nw2 Re=%-5g n=%d : %s" % (Re, n, st))
for n in (2, 4):
    for al in (0.05, 0.20):
        r = load(384, 1000, n, al)
        print("   a=%.2f Re=1000 n=%d : %s" % (al, n, "CONVERGE" if ok(r) else "a relancer"))

# ---------- 4) P4.1 : sweep Re a N=256 ----------
print("\n" + "=" * 78)
print("4) P4.1 — sweep Re a N=256, a=0 : Z(n=2), Z(n=4), exposant local")
print("=" * 78)
res = {}
for Re in (100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350, 375, 400):
    z2 = load(256, Re, 2, 0.0); z4 = load(256, Re, 4, 0.0)
    if ok(z2) and ok(z4):
        g = math.log(z2["Z"]/z4["Z"]) / math.log(4.0/2.0) if z2["Z"]>0 and z4["Z"]>0 else float("nan")
        res[Re] = (z2["Z"], z4["Z"], g)
        print("   Re=%-4g Z2=%.6f Z4=%.6f  gamma_24=%.3f" % (Re, z2["Z"], z4["Z"], g))
    elif ok(z2):
        print("   Re=%-4g Z2=%.6f Z4=—" % (Re, z2["Z"]))
    elif ok(z4):
        print("   Re=%-4g Z2=— Z4=%.6f" % (Re, z4["Z"]))
    else:
        print("   Re=%-4g : non converge" % Re)

# ---------- 5) P2.2 : factorial Ma ----------
print("\n" + "=" * 78)
print("5) P2.2 — factorial (tau,U) : eps=balance et variation sur la plage Ma")
print("=" * 78)
for Re in (100, 1000):
    for N in (128, 256):
        for n in (2, 4):
            for tag in ("_n%g_a0.00" % n if n != 2 else "",):
                p = os.path.join(MA, "masweep_Re%g_N%d%s.json" % (Re, N, ("_n%g_a0.00" % n) if (n!=2 or True) else ""))
            # nommage reel du script: "" si n=2 et alpha=0 sinon _n%g_a%.2f
            tag = "" if (n == 2.0 and 0.0 == 0.0) else "_n%g_a%.2f" % (n, 0.0)
            # le script n'ecrit "" que si n==2 ET alpha==0
            p = os.path.join(MA, "masweep_Re%d_N%d%s.json" % (Re, N, tag))
            if not os.path.exists(p):
                print("   Re=%-5g N=%d n=%d : pas encore de JSON" % (Re, N, n))
                continue
            with open(p, encoding="utf-8") as fh:
                d = json.load(fh)
            eps = [c["balance"] for c in d["cases"]]
            mas = [c["Ma"] for c in d["cases"]]
            if eps:
                print("   Re=%-5g N=%d n=%d : %d cas, Ma %.3f-%.3f, eps %.3f%% a %.3f%% (spread %.3f pt)"
                      % (Re, N, n, len(eps), min(mas), max(mas), min(eps)*100, max(eps)*100,
                         (max(eps)-min(eps))*100))
print("\nFIN BILAN")
