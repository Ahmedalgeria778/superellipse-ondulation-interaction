# -*- coding: utf-8 -*-
"""
audit/build_inventory.py
========================
Etape 1 du plan : audit READ-ONLY des deux bases existantes et production d'une
table unifiee des cas reellement disponibles.

Il ne modifie AUCUN fichier des projets d'origine. Il lit :

  Projet A  Lid_driven_ondulation
      results/scalar_table.json            (cas principaux convergés)
      results_fneq/scalar_table.json       (sous-ensemble avec populations f)
      results/convergence/resolution.csv   (etudes de resolution)

  Projet B  superellipse_LBM_PI_operator
      05_dataset/**/manifest*.csv          (toutes les campagnes)
      06_physics/results_190/physics_observables.csv  (analyse physique 190 cas)

Sorties (dans audit/) :
      inventaire_cas.csv        table unifiee ligne par ligne
      synthese_inventaire.json  compteurs et croisements

Usage :
    python audit/build_inventory.py
"""

from __future__ import annotations

import csv
import glob
import json
import os
import re
import sys
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
PROJECT = os.path.dirname(HERE)

PROJ_A = os.environ.get("SEGI_PROJ_A", "")  # external base project (not shipped); see README "Provenance"
PROJ_B = os.environ.get("SEGI_PROJ_B", "")  # external base project (not shipped); see README "Provenance"

# Colonnes de la table unifiee.
COLUMNS = [
    "project", "family", "case_id", "n", "alpha", "Nw", "Re", "N",
    "lid", "status", "converged", "fields_npz",
    "area_fluid", "area_exact", "K", "P", "D", "balance_pct",
    "umax", "source", "notes",
]


def _fnum(v, default=None):
    """Conversion float tolerante."""
    if v is None:
        return default
    try:
        s = str(v).strip()
        if s == "" or s.lower() in ("nan", "none", "null"):
            return default
        return float(s)
    except (TypeError, ValueError):
        return default


def _b(v):
    if isinstance(v, bool):
        return v
    return str(v).strip().lower() in ("true", "1", "yes", "converged")


def load_json(path):
    try:
        with open(path, "r", encoding="utf-8") as fh:
            return json.load(fh)
    except (OSError, ValueError):
        return None


def read_csv_rows(path):
    try:
        with open(path, "r", encoding="utf-8-sig", newline="") as fh:
            return list(csv.DictReader(fh))
    except OSError:
        return []


_CASE_RE = re.compile(r"se_n(?P<n>[0-9p]+)_Re(?P<Re>[0-9]+)_N(?P<N>[0-9]+)")


def superellipse_key(case_id):
    """Cle (n, Re, N) normalisee.

    Les deux bases du projet B utilisent deux ecritures differentes du meme cas :
        manifest            : se_n4p0_Re400_N128
        physics_observables : se_n4_Re400_N128
    Le nommage est donc normalise par (n, Re, N) pour pouvoir apparier.
    """
    if not case_id:
        return None
    m = _CASE_RE.search(str(case_id))
    if not m:
        return None
    n = float(m.group("n").replace("p", "."))
    return (n, float(m.group("Re")), int(m.group("N")))


# ----------------------------------------------------------------------
#  Projet A
# ----------------------------------------------------------------------
def collect_project_A(rows):
    """Retourne (rows, notes) pour le projet A."""
    base = os.path.join(PROJ_A, "results")
    notes = []

    # 1) table principale + fneq
    for tag, tpath in (("A:scalar_table",
                        os.path.join(base, "scalar_table.json")),
                       ("A:scalar_table_fneq",
                        os.path.join(PROJ_A, "results_fneq", "scalar_table.json"))):
        table = load_json(tpath)
        if not table:
            notes.append(f"absent/illisible: {tpath}")
            continue
        for r in table:
            cid = r.get("case")
            N = _fnum(r.get("N"))
            npz = os.path.join(base, cid or "", f"{cid}_N{int(N)}.npz") \
                if cid and N else ""
            rows.append({
                "project": "A", "family": "wavy_square",
                "case_id": cid, "n": None,
                "alpha": _fnum(r.get("h_over_H")), "Nw": _fnum(r.get("Nw")),
                "Re": _fnum(r.get("Re")), "N": N, "lid": r.get("lid"),
                "status": r.get("status"),
                "converged": _b(r.get("status") == "converged"),
                "fields_npz": os.path.exists(npz),
                "area_fluid": _fnum(r.get("fluid_area")),
                "area_exact": None,
                "K": _fnum(r.get("K_star")), "P": _fnum(r.get("P_star")),
                "D": _fnum(r.get("D_star")),
                "balance_pct": (100.0 * _fnum(r.get("balance"))
                                if _fnum(r.get("balance")) is not None else None),
                "umax": None,
                "source": tag,
                "notes": ("flat_height=%s" % r.get("flat_height")
                          if r.get("flat_height") not in (None, 0.0) else ""),
            })

    # 2) etude de resolution
    rpath = os.path.join(base, "convergence", "resolution.csv")
    for r in read_csv_rows(rpath):
        cid = r.get("case")
        N = _fnum(r.get("N"))
        variant = r.get("variant")
        rows.append({
            "project": "A", "family": "wavy_square_conv",
            "case_id": cid, "n": None,
            "alpha": _fnum(r.get("h_over_H")), "Nw": _fnum(r.get("Nw")),
            "Re": _fnum(r.get("Re")), "N": N,
            "lid": "regularized",
            "status": r.get("status"),
            "converged": _b(r.get("status") == "converged"),
            "fields_npz": True,
            "area_fluid": None, "area_exact": None,
            "K": _fnum(r.get("K_star")), "P": _fnum(r.get("P_star")),
            "D": _fnum(r.get("D_star")),
            "balance_pct": _fnum(r.get("balance_pct")),
            "umax": None,
            "source": "A:resolution.csv",
            "notes": f"variant={variant}",
        })
    if not os.path.exists(rpath):
        notes.append(f"absent: {rpath}")
    return notes


# ----------------------------------------------------------------------
#  Projet B
# ----------------------------------------------------------------------
def _npz_exists_B(case_id):
    for pat in ("data", "data_all190", "data_all100", "data_all120",
                "data_re_ext", "data_re_ext_b", "data_refine", "data_refine_seed",
                "data_conv", "data_conv256", "data_convHR", "data_corner"):
        p = os.path.join(PROJ_B, "05_dataset", pat, f"{case_id}.npz")
        if os.path.exists(p):
            return True
    return False


def collect_project_B(rows):
    notes = []

    # 1) tous les manifest*.csv de 05_dataset
    manifests = glob.glob(os.path.join(PROJ_B, "05_dataset", "**",
                                       "manifest*.csv"), recursive=True)
    seen = set()
    for m in manifests:
        for r in read_csv_rows(m):
            cid = r.get("case_id") or r.get("case") or r.get("id")
            if not cid or cid in seen:
                continue
            seen.add(cid)
            N = _fnum(r.get("N"))
            rows.append({
                "project": "B", "family": "superellipse",
                "case_id": cid, "n": _fnum(r.get("n")), "alpha": 0.0, "Nw": 0.0,
                "Re": _fnum(r.get("Re")), "N": N,
                "lid": "tangent_arc",
                "status": "converged",
                "converged": True,
                "fields_npz": _npz_exists_B(cid),
                "area_fluid": None, "area_exact": None,
                "K": None, "P": None, "D": None, "balance_pct": None,
                "umax": _fnum(r.get("umax_lbm")),
                "source": "B:" + os.path.relpath(m, PROJ_B).replace("\\", "/"),
                "notes": "",
            })
    if not manifests:
        notes.append("aucun manifest*.csv trouve dans 05_dataset")

    # 2) analyse physique des 190 cas (D, K, aire, coeur)
    obs = os.path.join(PROJ_B, "06_physics", "results_190",
                       "physics_observables.csv")
    obsmap = {}
    for r in read_csv_rows(obs):
        k = superellipse_key(r.get("case_id"))
        if k is not None:
            obsmap[k] = r
    for r in rows:
        if r["project"] != "B":
            continue
        o = obsmap.get(superellipse_key(r["case_id"]))
        if not o:
            continue
        r["area_fluid"] = _fnum(o.get("area_mask"))
        r["area_exact"] = _fnum(o.get("area_exact"))
        r["K"] = _fnum(o.get("kinetic_energy"))
        r["D"] = _fnum(o.get("dissipation"))
        # D_core conserve dans les notes (colonne dediee absente de la table unifiee)
        extras = []
        dcore = _fnum(o.get("dissipation_core"))
        if dcore is not None:
            extras.append("D_core=%.6g" % dcore)
        if str(o.get("pressure_available")).strip().lower() in ("false", "0", ""):
            extras.append("P_non_calculee")
        if extras:
            r["notes"] = (r["notes"] + " ; " + " ; ".join(extras)).strip(" ;")
        r["umax"] = _fnum(o.get("umax")) or r["umax"]
    if not os.path.exists(obs):
        notes.append(f"absent: {obs}")
    return notes


# ----------------------------------------------------------------------
def write_csv(path, rows):
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS, extrasaction="ignore")
        w.writeheader()
        for r in sorted(rows, key=lambda x: (str(x["project"]),
                                             str(x["family"]),
                                             str(x["case_id"]))):
            w.writerow(r)


def summarise(rows):
    by_proj = Counter(r["project"] for r in rows)
    by_family = Counter(r["family"] for r in rows)
    ns = sorted({r["n"] for r in rows if r["n"] is not None and r["project"] == "B"})
    reA = sorted({r["Re"] for r in rows if r["project"] == "A" and r["Re"] is not None})
    reB = sorted({r["Re"] for r in rows if r["project"] == "B" and r["Re"] is not None})
    alphas = sorted({r["alpha"] for r in rows if r["project"] == "A"
                     and r["alpha"] is not None})
    nws = sorted({r["Nw"] for r in rows if r["project"] == "A"
                  and r["Nw"] is not None})
    Ns = sorted({r["N"] for r in rows if r["N"] is not None})
    withfields = sum(1 for r in rows if r["fields_npz"])
    return {
        "n_rows": len(rows),
        "by_project": dict(by_proj),
        "by_family": dict(by_family),
        "rows_with_fields_npz": withfields,
        "project_B_exponents_n": ns,
        "project_A_Re": reA,
        "project_B_Re": reB,
        "project_A_alpha": alphas,
        "project_A_Nw": nws,
        "resolutions_N": Ns,
    }


def main():
    rows, notes = [], []
    if os.path.isdir(PROJ_A):
        notes += collect_project_A(rows)
    else:
        notes.append(f"projet A introuvable : {PROJ_A}")
    if os.path.isdir(PROJ_B):
        notes += collect_project_B(rows)
    else:
        notes.append(f"projet B introuvable : {PROJ_B}")

    out_csv = os.path.join(HERE, "inventaire_cas.csv")
    out_json = os.path.join(HERE, "synthese_inventaire.json")
    write_csv(out_csv, rows)
    synth = summarise(rows)
    synth["notes"] = notes
    with open(out_json, "w", encoding="utf-8") as fh:
        json.dump(synth, fh, indent=2, ensure_ascii=False)

    print(f"[inventaire] {len(rows)} lignes -> {out_csv}")
    print(f"[synthese]   -> {out_json}")
    print(json.dumps(synth, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
