# -*- coding: utf-8 -*-
"""
solver/plot_contrast_figure.py
===============================
Figure a deux panneaux demandee pour le relecteur :

  (a) |I_D(paire, Re)| pour les 3 paires de forme (2,4) (4,6) (6,8),
      aux 3 Reynolds et aux 2 resolutions (N=128 plein, N=256 ouvert),
      avec en fond la courbe continue des 19 Reynolds (paire 2,4, N=128) ;
  (b) decomposition ADDITIVE signee de l'effet en intervalles de forme
      2->4 / 4->6 / 6->8, en % du total, pour chaque Re et les 2 grilles.

Conventions de legende (demandees) :
  * couleur  = paire geometrique ;
  * trait/marqueur = grille (N=128 : plein + pastille ; N=256 : tiret + carre ouvert) ;
  * style propre (gris tirete) = resweep des 19 Reynolds ;
  * 2 legendes groupees d'UNE ligne chacune (pas de legende verticale a 7 entrees).

Entree  : results/audit/audit_contrasts.json  (audit A-F, priorite 1)
          results/resweep/resweep_N128.json    (19 Reynolds, optionnel)
Sortie  : figures/analysis/contraste_paire_Re.png (600 dpi) + .pdf (vecteur)
          results/audit/figure_contrast_data.json

Controles automatiques a l'export (voir report() ) :
  * chevauchement etiquette d'axe / legende ;
  * taille physique du PNG et corps de police effectif si le PDF est place
    a 190 mm de large (largeur de colonne d'un journal).

Usage :
  python -m solver.plot_contrast_figure
"""

from __future__ import annotations

import json
import os
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
from matplotlib.lines import Line2D                    # noqa: E402
from matplotlib.patches import Patch                   # noqa: E402
from matplotlib.ticker import FixedLocator              # noqa: E402

AUDIT_JSON = os.path.join(_ROOT, "results", "audit", "audit_contrasts.json")
RESWEEP_JSON = os.path.join(_ROOT, "results", "resweep", "resweep_N128.json")
FIG_DIR = os.path.join(_ROOT, "figures", "analysis")
AUDIT_DIR = os.path.join(_ROOT, "results", "audit")

RES = [100, 400, 1000]
PAIRS = [("24", r"$2\rightarrow4$", (2, 4)),
         ("46", r"$4\rightarrow6$", (4, 6)),
         ("68", r"$6\rightarrow8$", (6, 8))]
PAIR_COL = {"24": "C0", "46": "C1", "68": "C2"}

LS128 = "-"
LS256 = (0, (3.2, 1.8))
DASH_SWEEP = (0, (5.0, 2.2))

# tailles de police (points) -- voir report() pour le corps effectif imprime
FS_TITLE = 10.0
FS_LABEL = 10.0
FS_TICK = 9.0
FS_LEG = 9.0
FIGSIZE = (9.4, 4.5)
DPI_PNG = 600
JOURNAL_MM = 190.0          # largeur de placement visee (PDF)


def load():
    with open(AUDIT_JSON, encoding="utf-8") as fh:
        return json.load(fh)


def panel_a(ax, d):
    """|I_D| vs Re, 3 paires, 2 grilles, resweep en fond."""
    sweep = None
    if os.path.exists(RESWEEP_JSON):
        with open(RESWEEP_JSON, encoding="utf-8") as fh:
            sw = json.load(fh)
        pts = sw.get("points") or []
        if pts:
            sweep = (np.array([p["Re"] for p in pts], float),
                     np.array([abs(p["I"]) for p in pts], float))

    if sweep is not None:
        ax.plot(sweep[0], sweep[1], "-", color="0.86", lw=3.4, zorder=0)
        ax.plot(sweep[0], sweep[1], "--", color="0.45", lw=1.35,
                dash_capstyle="round", zorder=1)

    unc = d["uncertainty"]
    for key, lab, _pr in PAIRS:
        re128, y128, re256, y256 = [], [], [], []
        for Re in RES:
            rec = unc.get(str(Re), {}).get("D_star_" + key)
            if rec is None:
                continue
            re128.append(Re)
            y128.append(abs(rec["I_128"]))
            re256.append(Re)
            y256.append(abs(rec["I_256"]))
        if not y128:
            continue
        c = PAIR_COL[key]
        ax.plot(re128, y128, ls=LS128, color=c, marker="o", ms=5.4,
                mfc=c, mec=c, lw=1.5, zorder=3, label=lab)
        ax.plot(re256, y256, ls=LS256, color=c, marker="s", ms=5.6,
                mfc="none", mec=c, mew=1.5, lw=1.35, zorder=4, label=lab)

    ax.set_xscale("log")
    ax.set_yscale("log")
    # marge en haut : la zone [0.15, 1] ne contient aucune donnee (max 0.1466)
    # et sert de reserve a la legende interne (aucun recouvrement possible)
    ax.set_ylim(4.0e-5, 1.0)
    ax.set_xlabel("Re", fontsize=FS_LABEL)
    ax.set_ylabel(r"$|I_D|$  (contraste apparié)", fontsize=FS_LABEL)
    ax.set_title("(a) Contraste d'interaction vs Reynolds", fontsize=FS_TITLE)
    ax.grid(alpha=0.3, which="both")
    ax.xaxis.set_minor_locator(FixedLocator(
        [100, 150, 200, 300, 400, 600, 800, 1000]))
    ax.set_xticks([100, 400, 1000])
    ax.set_xticklabels(["100", "400", "1000"], fontsize=FS_TICK)
    ax.tick_params(axis="both", which="major", labelsize=FS_TICK)

    # legende A : les 3 paires, par couleur, UNE ligne (en haut a gauche,
    # dans la marge ouverte ci-dessus). Idiom a deux legends : on conserve
    # l'objet avant de creer la seconde, sinon il est ecrase.
    legA = [Line2D([], [], color=PAIR_COL[k], lw=2.4, label=lab)
            for k, lab, _pr in PAIRS]
    objA = ax.legend(handles=legA, fontsize=FS_LEG, ncol=3, loc="upper left",
                     framealpha=0.92, borderpad=0.45, handlelength=1.9,
                     columnspacing=1.0, handletextpad=0.5)
    ax.add_artist(objA)

    # legende B : style de grille + resweep, UNE ligne, sous l'axe
    legB = [
        Line2D([], [], color="0.35", ls=LS128, marker="o", ms=5.0,
               mfc="0.35", lw=1.5, label=r"$N=128$ (plein)"),
        Line2D([], [], color="0.35", ls=LS256, marker="s", ms=5.2,
               mfc="none", lw=1.5, label=r"$N=256$ (ouvert)"),
        Line2D([], [], color="0.45", ls=DASH_SWEEP, lw=1.6,
               label=r"19 Re, $2\rightarrow4$, $N=128$"),
    ]
    ax.legend(handles=legB, fontsize=FS_LEG, ncol=3, loc="upper center",
              bbox_to_anchor=(0.5, -0.30), framealpha=0.92,
              borderpad=0.5, columnspacing=1.3, handlelength=2.1,
              handletextpad=0.5)
    return objA


def panel_b(ax, d):
    """Decomposition additive signee en % du total, empilee."""
    shares = d["shares"]
    x = np.arange(len(RES))
    w = 0.38
    for k, N in enumerate((128, 256)):
        bottoms = np.zeros(len(RES))
        for iseg, (key, lab, _pr) in enumerate(PAIRS):
            vals = []
            for Re in RES:
                rec = (shares.get(str(Re), {})
                       .get("D* extensive", {})
                       .get(str(N)))
                vals.append(rec["pct"][iseg] if rec else 0.0)
            vals = np.asarray(vals, float)
            ax.bar(x + (k - 0.5) * w, vals, width=w, bottom=bottoms,
                   color=plt.cm.viridis(0.15 + 0.30 * iseg),
                   edgecolor="k" if N == 128 else "0.35",
                   linewidth=1.0 if N == 128 else 0.9,
                   hatch="" if N == 128 else "///")
            bottoms += vals

    ax.set_xticks(x)
    ax.set_xticklabels([str(r) for r in RES], fontsize=FS_TICK)
    # la mention des remplissages va dans l'etiquette d'axe, comme demande,
    # pour ne pas alourdir la legende (qui reste a 3 entrees / 1 ligne)
    ax.set_xlabel(r"Re  —  plein : $N=128$ ; hachuré : $N=256$",
                  fontsize=FS_LABEL)
    ax.set_ylim(0, 101)
    ax.set_yticks([0, 25, 50, 75, 100])
    ax.set_ylabel("Part additive de l'effet  (%)", fontsize=FS_LABEL)
    ax.set_title("(b) Décomposition 2→4 / 4→6 / 6→8", fontsize=FS_TITLE)
    ax.grid(alpha=0.3, axis="y")
    ax.set_axisbelow(True)
    ax.tick_params(axis="y", which="major", labelsize=FS_TICK)

    handles = [Patch(facecolor=plt.cm.viridis(0.15 + 0.30 * i),
                     edgecolor="k", label=PAIRS[i][1]) for i in range(3)]
    ax.legend(handles=handles, fontsize=FS_LEG, ncol=3, loc="upper center",
              bbox_to_anchor=(0.5, -0.30), framealpha=0.92,
              columnspacing=1.4, handlelength=1.7, borderpad=0.5,
              handletextpad=0.5)


def report(fig, axes, leg_a, png_path, pdf_ok, notes):
    """Controles d'export : chevauchements, taille physique, corps imprime."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    print("[chk] legendes / etiquette d'axe :")
    for i, ax in enumerate(axes):
        lg = ax.get_legend()
        lb = ax.xaxis.label
        bb_lg = lg.get_window_extent(r)
        bb_lb = lb.get_window_extent(r)
        overlap = bb_lg.overlaps(bb_lb)
        print("      panneau %s : legende-bas=%s, xlabel=%s, chevauchement=%s"
              % ("ab"[i],
                 tuple(int(v) for v in bb_lg.extents),
                 tuple(int(v) for v in bb_lb.extents),
                 "OUI ***" if overlap else "non"))
        if overlap:
            notes.append("xlabel/legende chevauches en panneau %s"
                         % "ab"[i])

    # legende interne (a) : verifie qu'elle ne recouvre AUCUNE courbe
    if leg_a is not None:
        bbA = leg_a.get_window_extent(r)
        hits = []
        for ln in axes[0].lines:
            p = ln.get_path()
            if len(p.vertices) == 0:
                continue
            if p.transformed(ln.get_transform()).intersects_bbox(bbA):
                hits.append(ln.get_label())
        print("[chk] legende interne (a) : %s"
              % ("recouvre %d courbe(s): %s" % (len(hits), hits)
                 if hits else "aucune courbe recouverte"))
        if hits:
            notes.append("legende interne recouvre : %s" % hits)

    w_in, h_in = fig.get_size_inches()
    for mm in (JOURNAL_MM,):
        f = (mm / 25.4) / w_in
        print("[chk] placement PDF a %.0f mm de large -> facteur %.3f"
              % (mm, f))
        print("      corps effectif imprime : titre %.1f pt, axe %.1f pt, "
              "graduation %.1f pt, legende %.1f pt"
              % (FS_TITLE * f, FS_LABEL * f, FS_TICK * f, FS_LEG * f))
        if FS_LEG * f < 6.0:
            notes.append("corps de legende imprime < 6 pt a %.0f mm" % mm)

    if os.path.exists(png_path):
        from matplotlib.image import imread
        img = imread(png_path)
        px_h, px_w = img.shape[0], img.shape[1]
        print("[chk] PNG %d x %d px  ->  %.1f x %.1f mm a %d dpi"
              % (px_w, px_h, px_w / DPI_PNG * 25.4, px_h / DPI_PNG * 25.4,
                 DPI_PNG))
    print("[chk] PDF : %s" % ("ecrit" if pdf_ok
                              else "NON ECRIT (fichier verrouille)"))
    if not pdf_ok:
        notes.append("pdf non ecrit : fichier verrouille")


def main():
    if not os.path.exists(AUDIT_JSON):
        print("audit introuvable : %s" % AUDIT_JSON)
        return 1
    d = load()

    os.makedirs(FIG_DIR, exist_ok=True)
    os.makedirs(AUDIT_DIR, exist_ok=True)

    fig, axes = plt.subplots(1, 2, figsize=FIGSIZE)
    leg_a = panel_a(axes[0], d)
    panel_b(axes[1], d)
    fig.tight_layout()
    notes = []

    out_png = os.path.join(FIG_DIR, "contraste_paire_Re.png")
    fig.savefig(out_png, dpi=DPI_PNG, bbox_inches="tight")
    out_pdf = os.path.join(FIG_DIR, "contraste_paire_Re.pdf")
    try:
        fig.savefig(out_pdf, bbox_inches="tight")
        pdf_ok = True
    except PermissionError:                 # ouvert dans une visionneuse
        pdf_ok = False
    report(fig, axes, leg_a, out_png, pdf_ok, notes)
    plt.close(fig)

    # --- donnees exposees (pour relecture / replot) -------------------------
    export = {"pairs": [p[1].replace("$", "").replace(r"\rightarrow", "->")
                        for p in PAIRS],
              "Re": RES, "data": {}}
    unc = d["uncertainty"]
    for key, _lab, _pr in PAIRS:
        export["data"][key] = {
            str(Re): {"I_128": unc[str(Re)]["D_star_" + key]["I_128"],
                      "I_256": unc[str(Re)]["D_star_" + key]["I_256"],
                      "ecart_pct": unc[str(Re)]["D_star_" + key]["ecart_pct"]}
            for Re in RES}
    out_json = os.path.join(AUDIT_DIR, "figure_contrast_data.json")
    with open(out_json, "w", encoding="utf-8") as fh:
        json.dump(export, fh, indent=2, ensure_ascii=False)

    print("[fig] %s" % os.path.relpath(out_png, _ROOT))
    if pdf_ok:
        print("[fig] %s" % os.path.relpath(out_pdf, _ROOT))
    else:
        print("[fig] %s NON ECRIT (fichier verrouille)" %
              os.path.relpath(out_pdf, _ROOT))
    print("[fig] %s" % os.path.relpath(out_json, _ROOT))
    for n in notes:
        print("[fig][WARN] %s" % n)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
