"""Couche de donnees unique pour le manuscrit et le supplement NEUFS a 120 cas.

Regle du projet : aucun nombre du TeX n'est tape a la main.  Tout ce que les
deux composeurs ecrivent dans leurs .tex passe par ce module, qui lit :

  * 06_physics/rapports/chiffres_120.json          grille, observables, argmin
  * 06_physics/figures_120/tables/*_120.csv|.json   analyse de paroi, monotonie
  * 06_physics/results_120/vortex_catalog.csv      catalogue des tourbillons
  * 06_physics/tables_120/*.tex                    fragments de tableaux deja
                                                   generes (on les inclut tels
                                                   quels, ils sont deja
                                                   produits a partir des
                                                   donnees)
  * 08_operator_results/A_B0_B1_summary.json      operateur, gele
  * 08_operator_results/PH13_final/*              operateur, tableaux
                                                   consolides

Le module se verifie lui-meme : `verifier()` recalcule chaque affirmation que
les documents将与 porteront et echoue bruyamment si l'une d'elles est fausse.
Rien n'est ecrit avant que `verifier()` soit passe.

Convention de dissipation, celle du papier :
    D = 2 nu integral(S:S dA),  nu = 2/Re,  D/A normalise par area_exact.
La colonne `dissipation` du CSV de production vaut (1/Re) integral(2 S:S) :
elle est donc multipliee par 2 pour retrouver D, et divisee par l'aire exacte
pour obtenir D/A.  C'est exactement ce que fait `T(n, Re)["D_papier"]`.
"""
from __future__ import annotations

import csv
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
PROJ = HERE.parent

GRILLE = 120   # taille de la grille ; pilote de tous les chemins
CHIFFRES = HERE / "rapports" / ("chiffres_%d.json" % GRILLE)
NEARWALL_CSV = (HERE / ("figures_%d" % GRILLE) / "tables"
                / ("PH13_nearwall_curvature_analysis_%d.csv" % GRILLE))
NEARWALL_JSON = (HERE / ("figures_%d" % GRILLE) / "tables"
                / ("PH13_nearwall_curvature_analysis_%d.json" % GRILLE))
MONO_JSON = (HERE / ("figures_%d" % GRILLE) / "tables"
             / ("PH13_nearwall_monotonicity_%d.json" % GRILLE))
CATALOGUE = HERE / ("results_%d" % GRILLE) / "vortex_catalog.csv"
TABLES = HERE / ("tables_%d" % GRILLE)

OPER_SUMMARY = PROJ / "08_operator_results" / "A_B0_B1_summary.json"
OPER_FIELD_R2 = PROJ / "08_operator_results" / "PH13_final" / "PH13_final_field_R2.csv"
OPER_PHYSICS_CSV = PROJ / "08_operator_results" / "PH13_final" / "PH13_final_physics.csv"
OPER_COST = PROJ / "08_operator_results" / "PH13_final" / "PH13_cost_curve.json"
OPER_B1_TR98 = (PROJ / "08_operator_results" / "B1_PI_train98_pinit"
                / "B1_results.json")

# Grille figee, telle que publiee : elle ne fait pas partie des 120 cas.
GRILLE_FIGEE = (100.0, 400.0, 1000.0)
N_CAS_FIGES = 30


# ---------------------------------------------------------------------------
# formatage -- aucune valeur n'est ecrite a la main dans les .tex
# ---------------------------------------------------------------------------
def g(x: float) -> str:
    """%g, comme dans le TeX publie : 100, 2.5, 1.25."""
    return "%g" % x


def dec(x: float, n: int = 4) -> str:
    return "%.*f" % (n, x)


def sci(x: float, n: int = 4) -> str:
    """$9.9990\\times10^{-9}$ ; mantisse et exposant sont calcules."""
    e = int(math.floor(math.log10(abs(x))))
    m = x / 10.0 ** e
    if round(m, n) >= 10.0:
        m, e = m / 10.0, e + 1
    return "$%.*f\\times10^{-%d}$" % (n, m, -e)


def pc(x: float, n: int = 1) -> str:
    """$1.5\\,\\%$ a partir d'une fraction."""
    return "$%.*f\\,\\%%$" % (n, 100.0 * x)


def pct(x: float, n: int = 1) -> str:
    """Comme pc mais sans dollar : pour un fragment deja en mode math."""
    return "%.*f\\,\\%%" % (n, 100.0 * x)


# ---------------------------------------------------------------------------
# chargement
# ---------------------------------------------------------------------------
_chiffres: dict | None = None
_nw: dict[tuple[float, float], dict] | None = None
_mono: dict | None = None
_cat: list[dict] | None = None
_oper: dict | None = None
_b1tr98: dict | None = None
_r2: dict[tuple[str, str], float] | None = None
_phys: dict[str, dict] | None = None


def chiffres() -> dict:
    global _chiffres
    if _chiffres is None:
        _chiffres = json.loads(CHIFFRES.read_text(encoding="utf-8"))
    return _chiffres


def NS() -> list[float]:
    return chiffres()["ns"]


def RES() -> list[float]:
    return chiffres()["res"]


def NCAS() -> int:
    return chiffres()["n_cases"]


def T(n: float, re: float) -> dict:
    """Une case de la grille, cles en %g."""
    return chiffres()["table"][g(n)][g(re)]


def argmin(re: float) -> dict:
    return chiffres()["argmin"][g(re)]


def marge(re: float) -> float:
    """Ecart relatif en % entre les deux meilleurs D/A a ce Reynolds."""
    return chiffres()["marge_argmin_DA_pct"][g(re)]


def qualite() -> dict:
    return chiffres()["qualite"]


def bandes() -> dict:
    return chiffres()["bandes"]


def W(n: float, re: float) -> dict:
    """Ligne d'analyse de paroi, colonnes converties en flottants."""
    global _nw
    if _nw is None:
        _nw = {}
        with NEARWALL_CSV.open(encoding="utf-8", newline="") as fh:
            for row in csv.DictReader(fh):
                cle = (float(row["n"]), float(row["Re"]))
                _nw[cle] = {k: (v if k == "case_id" else float(v))
                            for k, v in row.items()}
    return _nw[(n, re)]


def nw_rows(re: float) -> list[dict]:
    return [W(n, re) for n in sorted(NS())]


def mono() -> dict:
    global _mono
    if _mono is None:
        _mono = json.loads(MONO_JSON.read_text(encoding="utf-8"))["diagnostic"]
    return _mono


def catalogue() -> list[dict]:
    global _cat
    if _cat is None:
        with CATALOGUE.open(encoding="utf-8", newline="") as fh:
            _cat = [dict(r, n=float(r["n"]), Re=float(r["Re"]),
                         threshold_pct=float(r["threshold_pct"]),
                         rank=int(r["rank"]), psi=float(r["psi"]),
                         abs_psi=float(r["abs_psi"]))
                    for r in csv.DictReader(fh)]
    return _cat


def oper() -> dict:
    global _oper
    if _oper is None:
        _oper = json.loads(OPER_SUMMARY.read_text(encoding="utf-8"))
    return _oper


def r2(model: str, split: str) -> float:
    global _r2
    if _r2 is None:
        _r2 = {}
        with OPER_FIELD_R2.open(encoding="utf-8", newline="") as fh:
            for r in csv.DictReader(fh):
                _r2[(r["model"], r["split"])] = float(r["R2"])
    return _r2[(model, split)]


def phys_op(model: str) -> dict:
    """Ligne de PH13_final_physics.csv, indexee par modele."""
    global _phys
    if _phys is None:
        _phys = {}
        with OPER_PHYSICS_CSV.open(encoding="utf-8", newline="") as fh:
            for r in csv.DictReader(fh):
                _phys[r["model"]] = r
    return _phys[model]


# ---------------------------------------------------------------------------
# grandeurs derivees, chacune verifiee plus bas
# ---------------------------------------------------------------------------
def tuples_argmin() -> tuple[str, str]:
    """n^*_D et n^*_{D/A} ecrits sur les 12 Reynolds, dans l'ordre de la grille.

    Le papier publie ecrivait deux triplets (trois Reynolds).  Sur la grille
    complete ce sont deux suites de %d valeurs ; les ecrire suches serait faux.
    """ % NCAS()
    td = ",\\,".join(g(argmin(r)["n_star_D_papier"]) for r in RES())
    ta = ",\\,".join(g(argmin(r)["n_star_D_A_papier"]) for r in RES())
    return td, ta


def re_basculement_D() -> float:
    """Plus petit Re tel que l'argmin de D quitte 2.5."""
    for r in RES():
        if argmin(r)["n_star_D_papier"] != argmin(RES()[0])["n_star_D_papier"]:
            return r
    raise AssertionError("l'argmin de D ne bascule sur aucun Reynolds")


def re_basculement_DA() -> float:
    for r in RES():
        if argmin(r)["n_star_D_A_papier"] != argmin(RES()[0])["n_star_D_A_papier"]:
            return r
    raise AssertionError("l'argmin de D/A ne bascule sur aucun Reynolds")


def re_marge_minimale() -> tuple[float, float]:
    """Reynolds ou les deux meilleurs D/A sont les plus proches, et l'ecart."""
    r = min(RES(), key=marge)
    return r, marge(r)


def montee_presquarree() -> dict[float, float]:
    """D/A de n=8 relatif a n=4, en %, par Reynolds.

    Base FIXE n=4, et non l'argmin echantillonne.  Les deux ne coincident qu'a
    partir de Re=300 : a Re=100 l'argmin est n=3, et l'ecart vaut 20.49 %
    au lieu de 18.03 %.  Le papier compare n=4 a n=8, donc c'est la base fixe
    qu'il faut ; la verifier() verrouille cette distinction, parce que c'est
    precisement le genre d'ecart qui ferait porter au papier un chiffre faux.
    """
    out = {}
    for r in RES():
        base = T(4.0, r)["D_A_papier"]
        out[r] = 100.0 * (T(8.0, r)["D_A_papier"] - base) / base
    return out


def montee_depuis_argmin() -> dict[float, float]:
    """D/A de n=8 relatif a l'argmin echantillonne, en %.

    Cite tel quel seulement avec sa base explicite, car il differe de
    `montee_presquarree` aux Reynolds bas.
    """
    out = {}
    for r in RES():
        base = T(argmin(r)["n_star_D_A_papier"], r)["D_A_papier"]
        out[r] = 100.0 * (T(8.0, r)["D_A_papier"] - base) / base
    return out


def meilleur_re() -> float:
    return max(RES())


def pire_re() -> float:
    return min(RES())


def n_avec_tourbillon_secondaire(seuil: float = 1.0) -> list[tuple[float, float, int]]:
    """Cas ou plus de `seuil`% du tourbillon primaire est primaire."""
    out = []
    for n in NS():
        for r in RES():
            c = T(n, r)
            if c["n_vortices_1pct"] > 1:
                out.append((n, r, c["n_vortices_1pct"]))
    return out


def familles() -> dict[tuple[float, float], dict]:
    """Familles de violation de f_0.05, datees.

    Le JSON de monotonie donne, pour chaque Reynolds, les paires (n_from, n_to)
    dont le gain est positif.  Regroupees par paire, elles doivent former un
    intervalle de Reynolds : c'est ce que la datation verifie ici, parce que le
    papier affirmera que chaque famille apparaît à partir d'un Reynolds
    identifié et non à un seul endroit de la grille.
    """
    out: dict[tuple[float, float], dict] = {}
    for rtxt, items in mono()["0.05"]["violations"].items():
        r = float(rtxt)
        for it in items:
            cle = (it["n_from"], it["n_to"])
            e = out.setdefault(cle, {"res": [], "max_pct": 0.0, "max_re": None})
            e["res"].append(r)
            if it["rise_pct"] > e["max_pct"]:
                e["max_pct"] = it["rise_pct"]
                e["max_re"] = r
    for e in out.values():
        e["res"] = sorted(e["res"])
        assert e["res"] == [r for r in RES() if r in set(e["res"])], (
            "les Reynolds d'une famille ne suivent pas l'ordre de la grille")
        precedents = [r for r in RES() if r < e["res"][0]]
        e["re_avant"] = max(precedents) if precedents else None
        e["re_debut"] = e["res"][0]
        e["re_fin"] = e["res"][-1]
    return out


def famille_presente(cle: tuple[float, float], re: float) -> bool:
    return cle in familles() and re in familles()[cle]["res"]


def tableau_correlations(re: float) -> list[tuple[str, float, float, str, str]]:
    """(predicteur, r de Pearson, rho de Spearman, libelle, unite) a un Re.

    Les deux statistiques sont recalculees ici, et non relues : l'audit du
    papier publie compare des valeurs, pas la methode.
    """
    lignes = nw_rows(re)
    y = [T(w["n"], re)["D_A_papier"] for w in lignes]
    preds = (("kappa_w_rms", r"$\kappa_{w,\rm rms}$", ""),
             ("kappa_D_weighted", r"$\langle\kappa\rangle_D$", ""),
             ("f_delta_025", r"$f_{0.025}$", ""),
             ("f_delta_05", r"$f_{0.05}$", ""),
             ("f_sharp_total", r"$f_{\rm sharp}$", ""),
             ("A", r"$A$", ""))
    out = []
    for cle, lab, unit in preds:
        x = [w[cle] for w in lignes]
        out.append((cle, pearson(x, y), spearman(x, y), lab, unit))
    return out


def pearson(x: list[float], y: list[float]) -> float:
    n = len(x)
    mx, my = sum(x) / n, sum(y) / n
    sxy = sum((a - mx) * (b - my) for a, b in zip(x, y))
    sxx = sum((a - mx) ** 2 for a in x)
    syy = sum((b - my) ** 2 for b in y)
    return sxy / math.sqrt(sxx * syy) if sxx > 0 and syy > 0 else float("nan")


def spearman(x: list[float], y: list[float]) -> float:
    return pearson(ranks(x), ranks(y))


def ranks(v: list[float]) -> list[float]:
    ordre = sorted(range(len(v)), key=lambda i: v[i])
    r = [0.0] * len(v)
    i = 0
    while i < len(ordre):
        j = i
        while j + 1 < len(ordre) and v[ordre[j + 1]] == v[ordre[i]]:
            j += 1
        moy = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            r[ordre[k]] = moy
        i = j + 1
    return r


def corrs(predicteur: str) -> dict[float, tuple[float, float]]:
    """{Re: (r de Pearson, rho de Spearman)} pour un predicteur donne."""
    return {r: [t for t in tableau_correlations(r) if t[0] == predicteur][0][1:3]
            for r in RES()}


def pearson_extremes(predicteur: str) -> tuple[float, float, float, float]:
    """(r min, Re du r min, r max, Re du r max)."""
    c = corrs(predicteur)
    lo = min(c, key=lambda r: c[r][0])
    hi = max(c, key=lambda r: c[r][0])
    return c[lo][0], lo, c[hi][0], hi


def rho_min_global(predicteur: str) -> tuple[float, float]:
    """(rho minimal sur la grille, Reynolds ou il est atteint)."""
    c = corrs(predicteur)
    r = min(c, key=lambda r: c[r][1])
    return c[r][1], r


def bandes_KE_A() -> tuple[float, float, list[float]]:
    b = bandes()["KE_A_exact"]
    return b["min"], b["max"], b["band_a_citer"]


def aire(n: float) -> float:
    """Aire exacte du domaine,/demo-axes 1, donc facteur d'echelle 1.

    A(n) = 4 Gamma(1+1/n)^2 / Gamma(1+2/n).  La formule reproduit les
    120 aires de la grille a l'identite machine, y compris pour n<2 ou la
    forme est concave : verifier() la reverrouille sur les dix exposants.
    """
    return 4.0 * math.gamma(1.0 + 1.0 / n) ** 2 / math.gamma(1.0 + 2.0 / n)


def d_a_min(re: float) -> float:
    """D/A au minimum echantillonne de ce Reynolds, convention papier.

    C'est la serie des trois valeurs protegees du papier publie, une par
    Reynolds de la grille figee ; verifier() la reverrouille.
    """
    return T(argmin(re)["n_star_D_A_papier"], re)["D_A_papier"]


def r2_mesure(split: str, cle: str = "r2_field_mean") -> tuple[float, float]:
    """(moyenne, ecart-type) sur les trois graines, pour B1-train98.

    Le resume consolide ne donne que la moyenne ; l'ecart-type que publie le
    papier n'est que dans le B1_results.json de l'experimentation.  C'est donc
    la qu'il est lu, et non recalcule a partir d'un seul point.
    """
    global _b1tr98
    if _b1tr98 is None:
        _b1tr98 = json.loads(OPER_B1_TR98.read_text(encoding="utf-8"))
    s = _b1tr98["summary"][split][cle]
    return s["mean"], s["std"]


def r2_chaine(split: str, ndec: int = 4) -> tuple[str, str]:
    """(moyenne formee, ecart-type forme), arrondis au plus proche.

    A ndec=5 on obtient 0.97542 et 0.00145, dont la forme protegee du papier
    a quatre decimales est 0.9754 et 0.0015.  Le script ne fabrique pas la
    chaine protegee : il verifie seulement que les deux arrondis sont
    compatibles, pour que l'arrondi protege ne derive jamais en silence des
    donnees.
    """
    m, s = r2_mesure(split)

    def arrondi(x: float) -> str:
        return "%.*f" % (ndec, math.floor(x * 10 ** ndec + 0.5) / 10 ** ndec)

    return arrondi(m), arrondi(s)


def d_a_min_figee() -> dict[float, float]:
    """Les minima de D/A aux trois Reynolds de la grille figee."""
    return {r: d_a_min(r) for r in GRILLE_FIGEE}


def fragment(nom: str) -> str:
    """Un fragment de tableau deja genere, inclus tel quel."""
    return (TABLES / nom).read_text(encoding="utf-8").rstrip("\n")


# ---------------------------------------------------------------------------
# campagne par cas : ce tableau d'observables est derivable, le tableau
# publie ne l'etait pas : il portait le champ moyen sur les trois graines,
# que rien dans les artefacts ne permet de regenerer sans rejouer l'operateur).
# ---------------------------------------------------------------------------
CAMPAGNE = PROJ / "08_operator_results" / "PH13_final" / "PH13_case_corr_campaign.csv"
CAMPAGNE_NOM = {"B0": "B0", "B1-PI": "B1_PI_pinit_seed40",
                "B1-train98": "B1_PI_train98_pinit_seed40",
                "B1-rand98": "B1_PI_rand98_pinit_seed40"}
#: Lignes de tab:obs derivees de la campagne, et leur identifiant de source.
CAMPAGNE_METRIQUES = (("r2_field", r"\Rtwo", False),
                      ("ke_rel_err", r"$KE$", True),
                      ("D_rel_err", r"$D$", True),
                      ("E_rel_err", r"$\mathcal{E}$", True))
_camp: dict[str, list[dict]] | None = None

#: L'ancien CSV des 30 cas figes, verrouille, avec la correspondance de
#: colonnes.  `dissipation` y est en convention brute, d'ou le mapping vers
#: D_csv et non D_papier ; confondre les deux ferait croire a un ecart de 100 %.
ANCIEN_CSV = HERE / "results" / "physics_observables.csv"
ANCIEN_MAP = {"kinetic_energy": "KE", "dissipation": "D_csv",
              "umax": "umax", "vorticity_max_abs": "vorticity_max_abs",
              "area_exact": "A_exact", "area_mask": "A_mask"}


def verifier_grille_figee() -> tuple[int, int, float, int]:
    """Compare les 30 cas figes a l'ancien CSV verrouille.

    Rend (valeurs comparees, valeurs identiques au bit, ecart relatif
    maximal, nombre de cas).  C'est ce controle qui autorise le papier a
    ecrire que les 90 simulations nouvelles n'ont pas perturbe la grille
    figee sur laquelle l'operateur a ete construit.
    """
    ancien: dict[tuple[float, float], dict] = {}
    with ANCIEN_CSV.open(encoding="utf-8", newline="") as fh:
        for r in csv.DictReader(fh):
            ancien[(float(r["n"]), float(r["Re"]))] = r
    assert len(ancien) == N_CAS_FIGES, (
        "l'ancien CSV devrait contenir %d cas, il en contient %d"
        % (N_CAS_FIGES, len(ancien)))
    total = ident = 0
    pire = 0.0
    for (n, re), row in ancien.items():
        assert re in GRILLE_FIGEE, "le cas fige n=%g Re=%g sort de la grille" % (n, re)
        case = T(n, re)
        for colonne, cible in ANCIEN_MAP.items():
            if colonne not in row:
                continue
            a, b = float(row[colonne]), case[cible]
            total += 1
            if a == b:
                ident += 1
            else:
                pire = max(pire, abs(a - b) / max(abs(a), 1e-300))
    return total, ident, pire, len(ancien)


def campagne(modele: str) -> list[dict]:
    global _camp
    if _camp is None:
        _camp = {}
        with CAMPAGNE.open(encoding="utf-8", newline="") as fh:
            for r in csv.DictReader(fh):
                if r["split"] == "ood":
                    _camp.setdefault(r["model"], []).append(r)
    return _camp[CAMPAGNE_NOM[modele]]


def campagne_stat(modele: str, metrique: str, relatif: bool) -> tuple[float, float]:
    """(moyenne, pire) sur les configurations OOD, en pourcentage.

    `relatif` indique si la colonne est une erreur relative deja en fraction :
    elle est alors convertie en pourcentage, comme les lignes d'energie,
    dissipation et enstrophie du tableau des observables.
    """
    v = [float(r[metrique]) for r in campagne(modele) if r[metrique] not in ("", "nan")]
    assert v, "aucune valeur pour %s / %s" % (modele, metrique)
    k = 100.0 if relatif else 1.0
    return k * (sum(v) / len(v)), k * max(v)


# ---------------------------------------------------------------------------
# operateur : moyennes sur les trois realisations, calculees
# ---------------------------------------------------------------------------
def div_ood_moyen(modeles: tuple[str, ...]) -> float:
    """Moyenne, sur les trois graines, du Ldiv OOD de `modeles`."""
    vals = []
    for m in modeles:
        d = oper()["physics_diagnostic"]["models"].get(m)
        if d is not None:
            vals.append(d["ood"]["L_div_mean"])
    if not vals:
        raise AssertionError("aucun modele de %r dans le resume operateur" % (modeles,))
    return sum(vals) / len(vals)


def div_split(modele: str, split: str) -> float:
    return oper()["physics_diagnostic"]["models"][modele][split]["L_div_mean"]


def n_graines() -> int:
    return oper()["ExpB1_PI_field"]["n_seeds"]


#: Les trois graines, et le patron de cle de chaque modele dans le resume.
GRAINES = (40, 41, 42)
MOTIFS = {
    "B0": "B0",
    "B1": "B1_PI_seed%d",
    "B1-pinit": "B1_PI_pinit_seed%d",
    "B1-train98": "B1_PI_train98_pinit_seed%d",
}
#: Decoupage lisible -> decoupage des cles du resume operateur.
SPLITS = {"train": "train", "ID": "id_test", "OOD": "ood"}


def oper_blocs(modele: str) -> list[dict]:
    """Les blocs d'un modele, un par graine (un seul bloc pour B0)."""
    motif = MOTIFS[modele]
    if "%d" not in motif:
        return [oper()["physics_diagnostic"]["models"][motif]]
    return [oper()["physics_diagnostic"]["models"][motif % s] for s in GRAINES]


def oper_ldiv(modele: str, split: str, metrique: str = "L_div_mean") -> tuple[float, float, int]:
    """(moyenne, ecart-type, n cas) d'une metrique de divergence, par graine."""
    vals = [b[SPLITS[split]][metrique] for b in oper_blocs(modele)]
    n = oper_blocs(modele)[0][SPLITS[split]]["n_cases"]
    mu, sd = moyenne_ecart(vals)
    return mu, sd, n


def n_cases(split: str) -> int:
    """Nombre de cas du decoupage, tel que compte par le solveur de reference."""
    return oper()["physics_diagnostic"]["lbm_fd_div"][SPLITS[split]]["n"]


def moyenne_ecart(v: list[float]) -> tuple[float, float]:
    m = sum(v) / len(v)
    if len(v) < 2:
        return m, 0.0
    return m, math.sqrt(sum((x - m) ** 2 for x in v) / (len(v) - 1))


def ecart_type(v: list[float]) -> float:
    return moyenne_ecart(v)[1]


# ---------------------------------------------------------------------------
# verification : chaque affirmation est recalculee avant toute ecriture
# ---------------------------------------------------------------------------
def verifier() -> list[str]:
    notes: list[str] = []
    ch = chiffres()

    # --- coherence de la grille -------------------------------------------
    assert ch["n_cases"] == len(NS()) * len(RES()), "grille incoherente"
    assert len(NS()) == 10 and len(RES()) * 10 == GRILLE, (
        "grille 10x%d attendue, trouvee 10x%d" % (GRILLE // 10, len(RES())))
    notes.append("grille : %d cas = %d exponents x %d Reynolds"
                 % (ch["n_cases"], len(NS()), len(RES())))

    # --- la grille figee est bien un sous-ensemble -------------------------
    for r in GRILLE_FIGEE:
        assert r in RES(), "le Reynolds fige %g manque a la grille" % r
    notes.append("grille figee : %d Reynolds x %d exponents = %d cas inclus"
                 % (len(GRILLE_FIGEE), len(NS()), N_CAS_FIGES))

    # --- bascule de l'argmin de D/A ---------------------------------------
    r_da = re_basculement_DA()
    td, ta = tuples_argmin()
    notes.append("argmin D/A : %s (bascule a Re=%s)" % (ta, g(r_da)))
    assert ta.split(",\\,")[0] != ta.split(",\\,")[-1], "argmin D/A constant"
    notes.append("argmin D   : %s (bascule a Re=%s)" % (td, g(re_basculement_D())))
    assert td.split(",\\,")[0] != td.split(",\\,")[-1], "argmin D constant"

    # --- la marge minimale est bien la degenerescence annoncee -------------
    r_min, m_min = re_marge_minimale()
    assert m_min < 0.05, "la marge minimale %.4f %% n'est pas degéneree" % m_min
    des = ch["desaccords_normalisation_aire"]
    assert len(des) == 1 and des[0]["Re"] == r_min, (
        "le desaccord de normalisation devrait etre unique et a Re=%g" % r_min)
    notes.append("degenerescence : Re=%s, ecart %.4f %%, et le seul desaccord "
                 "de normalisation de la grille" % (g(r_min), m_min))

    # --- monotonie ---------------------------------------------------------
    d025 = mono()["0.025"]
    d05 = mono()["0.05"]
    assert len(d025["monotone_at_Re"]) == len(RES()) and not d025["violations"], (
        "f_0.025 devrait etre strictement decroissant a tous les Reynolds")
    notes.append("f_0.025 : decroissant aux %d Reynolds" % len(RES()))
    assert len(d05["monotone_at_Re"]) == 1, "f_0.05 : monotonie inattendue"
    notes.append("f_0.05 : monotone a %d Reynolds sur %d, %d familles datees"
                 % (len(d05["monotone_at_Re"]), len(RES()), len(familles())))
    for (a, b), e in sorted(familles().items()):
        assert e["re_avant"] is not None, "famille %g->%g presente des le debut" % (a, b)
        notes.append("   %g->%g : absente a Re=%s, presente de Re=%s a Re=%s, "
                     "max %+.4f %% a Re=%s"
                     % (a, b, g(e["re_avant"]), g(e["re_debut"]),
                        g(e["re_fin"]), e["max_pct"], g(e["max_re"])))

    # --- tourbillon secondaire --------------------------------------------
    sec = n_avec_tourbillon_secondaire()
    assert len(sec) == 1, "attendu un seul cas secondaire, obtenu %d" % len(sec)
    n_s, r_s, _ = sec[0]
    cat = [c for c in catalogue()
           if c["n"] == n_s and c["Re"] == r_s and c["threshold_pct"] == 1.0
           and c["rank"] == 2]
    assert len(cat) == 1, "le catalogue ne contient pas le secondaire attendu"
    part = 100.0 * cat[0]["abs_psi"] / abs(
        [c for c in catalogue()
         if c["n"] == n_s and c["Re"] == r_s and c["threshold_pct"] == 1.0
         and c["rank"] == 1][0]["psi"])
    notes.append("tourbillon secondaire : 1 cas sur %d (n=%g, Re=%g), "
                 "%.4f %% du primaire" % (ch["n_cases"], n_s, r_s, part))
    # absent a Re=900 : l'apparition est donc localisee dans ]900,1000]
    assert T(n_s, 900.0)["n_vortices_1pct"] == 1, "le secondaire ne devrait "\
        "pas etre present a Re=900"
    notes.append("apparition localisee dans ]%g,%g]" % (900.0, meilleur_re()))

    # --- qualite -----------------------------------------------------------
    q = ch["qualite"]
    assert q["cas_echoues"] == 0, "des cas ont echoue"
    assert q["area_rel_error_max"] < q["area_rel_error_seuil"], "aire hors seuil"
    assert q["residual_lbm_max"] < q["residual_lbm_seuil"], "residu hors seuil"
    assert q["rho_min"] > 0.0, "densite negative"
    notes.append("qualite : 0 echec, aire max %.4e, residuel max %.4e, rho min %.4f"
                 % (q["area_rel_error_max"], q["residual_lbm_max"], q["rho_min"]))

    # --- bande d'energie ---------------------------------------------------
    kmin, kmax, cite = bandes_KE_A()
    assert cite[0] <= kmin and kmax <= cite[1], "la bande citee ne contient pas KE/A"
    notes.append("KE/A dans [%.4f,%.4f], contenu dans la bande citee [%g,%g]"
                 % (kmin, kmax, cite[0], cite[1]))

    # --- correlations : ce que le papier en dit ----------------------------
    r_min, re_min, r_max, re_max = pearson_extremes("kappa_D_weighted")
    rho_min, re_rho = rho_min_global("kappa_D_weighted")
    assert r_min > 0 and rho_min < 0, (
        "Pearson devrait rester positif et Spearman changer de signe")
    assert re_min != re_rho, (
        "Pearson et Spearman localisent ici le meme Reynolds : le papier "
        "affirme le contraire")
    notes.append("corr. <k>_D : Pearson de %+.3f (Re=%g) a %+.3f (Re=%g) ; "
                 "Spearman de %+.3f a %+.3f (minimum a Re=%g ; Pearson et "
                 "Spearman n'y aboutissent pas au meme Reynolds)"
                 % (r_max, re_max, r_min, re_min,
                    corrs("kappa_D_weighted")[re_max][1], rho_min, re_rho))
    for r in RES():
        rangs_p = [t[0] for t in sorted(tableau_correlations(r), key=lambda t: -t[1])]
        rangs_s = [t[0] for t in sorted(tableau_correlations(r), key=lambda t: -t[2])]
        assert rangs_p != rangs_s, "Pearson et Spearman classeraient pareil a Re=%g" % r
    notes.append("Pearson et Spearman classent les predicteurs differemment "
                 "aux %d Reynolds" % len(RES()))

    # --- montee vers le carre : base fixe n=4, valeur protegee 18.0 % -----
    mp = montee_presquarree()
    r_bas = pire_re()
    assert abs(round(mp[r_bas], 1) - 18.0) < 1e-9, (
        "la montee vers le carre a Re=%g devrait valoir 18.0 %%, vaut %.4f %%"
        % (r_bas, mp[r_bas]))
    assert mp[meilleur_re()] < mp[r_bas], "la penalite presquarree ne s'affaiblit pas"
    notes.append("montee n=4 -> n=8 : %+.1f %% a Re=%g, %+.1f %% a Re=%g "
                 "(base fixe n=4 ; depuis l'argmin ce serait %+.1f %% a Re=%g)"
                 % (mp[r_bas], r_bas, mp[meilleur_re()], meilleur_re(),
                    montee_depuis_argmin()[r_bas], r_bas))

    # --- geometrie : l'aire exacte, formulee et mesuree --------------------
    for n in NS():
        for r in RES():
            assert aire(n) == T(n, r)["A_exact"], (
                "la formule d'aire ne reproduit pas A_exact a n=%g Re=%g" % (n, r))
    a_min, a_max = min(aire(n) for n in NS()), max(aire(n) for n in NS())
    notes.append("aire exacte : 4 Gamma(1+1/n)^2 / Gamma(1+2/n) reproduit les "
                 "%d valeurs de la grille a l'identite machine, de %.6f "
                 "(n=%g) a %.6f (n=%g)"
                 % (len(NS()) * len(RES()), a_min, min(NS(), key=aire),
                    a_max, max(NS(), key=aire)))

    # --- grille figee se reproduit dans la grille complete ------------------
    total, ident, pire, ncas = verifier_grille_figee()
    assert pire < 1e-14, (
        "les %d cas figes ne se reproduisent pas : ecart relatif %.3e"
        % (ncas, pire))
    notes.append("grille figee : %d valeurs comparees sur %d cas, %d "
                 "identiques au bit, ecart relatif maximal %s (erreur machine)"
                 % (total, ncas, ident, sci(pire, 3)))

    # --- valeurs protegees du papier publie, reverrouillees sur la grille ---
    protegees = {100.0: 3.228e-2, 400.0: 1.291e-2, 1000.0: 6.901e-3}
    for r, attendu in protegees.items():
        v = d_a_min(r)
        assert abs(v - attendu) / attendu < 2e-4, (
            "le minimum de D/A a Re=%g vaut %.4e, le papier publie %.4e"
            % (r, v, attendu))
    notes.append("valeurs de D/A au minimum, aux trois Reynolds figes : " +
                 ", ".join("%s -> %s" % (g(r), sci(d_a_min(r), 4))
                           for r in GRILLE_FIGEE))

    m5, s5 = r2_chaine("ood", 5)
    m4, s4 = r2_chaine("ood", 4)
    assert (m5, s5) == ("0.97542", "0.00145"), (m5, s5)
    # Le papier publie ecrit \Rtwo=0.9754\pm0.0015.  L'ecart-type reel est
    # 0.0014478 : arrondi directement a quatre decimales il donne 0.0014.  Le
    # 0.0015 du publie vient d'un double arrondi, cinq decimales puis quatre.
    # Le manuscrit neuf porte la valeur correcte ; l'ecart est documente ici
    # et verrouille, pour que la correction soit visible et non silencieuse.
    assert (m4, s4) == ("0.9754", "0.0014"), (m4, s4)
    brut_m, brut_s = r2_mesure("ood")
    notes.append("R2 OOD B1-train98 sur trois graines : %s +/- %s a cinq "
                 " decimales, %s +/- %s a quatre (ecart-type reel %.7f ; le "
                 "papier publie 0.0015 par double arrondi, le manuscrit neuf "
                 "porte 0.0014)"
                 % (m5, s5, m4, s4, brut_s))

    # --- operateur : les deux nombres proteges -----------------------------
    d_pinit = div_ood_moyen(("B1_PI_pinit_seed40", "B1_PI_pinit_seed41",
                             "B1_PI_pinit_seed42"))
    d_tr98 = div_ood_moyen(("B1_PI_train98_pinit_seed40",
                            "B1_PI_train98_pinit_seed41",
                            "B1_PI_train98_pinit_seed42"))
    assert 0.049 < d_pinit < 0.052, "Ldiv OOD B1-pinit inattendu : %g" % d_pinit
    assert 0.0060 < d_tr98 < 0.0064, "Ldiv OOD B1-train98 inattendu : %g" % d_tr98
    assert d_tr98 < d_pinit / 8.0, "le gain d'un ordre de grandeur ne se lit plus"
    notes.append("operateur : Ldiv OOD moyen %s (B1-pinit) -> %s (B1-train98), "
                 "facteur %.1f" % (sci(d_pinit, 2), sci(d_tr98, 3), d_pinit / d_tr98))
    assert n_graines() == 3, "trois realisations attendues"
    notes.append("operateur : %d realisations, R2 OOD B1-train98 = %.5f"
                 % (n_graines(), r2("B1-train98", "OOD")))

    return notes


if __name__ == "__main__":
    for ligne in verifier():
        print("  OK ", ligne)
    print("\n%d verifications passees" % len(verifier()))
