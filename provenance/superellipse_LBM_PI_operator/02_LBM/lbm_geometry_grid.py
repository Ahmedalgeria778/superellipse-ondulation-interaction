# -*- coding: utf-8 -*-
"""
02_LBM/lbm_geometry_grid.py

Interface geometrie <-> grille LBM pour la cavite superelliptique.

Produit, sur une grille N x N couvrant la boite englobante [-a,a] x [-b,b]
(de indices (i,j), i<->y, j<->x, comme le solveur de reference), :

  - le masque fluide/solide (interieur de Omega_n),
  - la liste des liens de frontiere (noeud fluide -> voisin solide) avec la
    fraction de paroi exacte q in (0,1) obtenue par bissection sur la courbe,
  - la vitesse de paroi associee (couvercle tangent pour n < infini).

Aucune modification du solveur de reference. Le schema de rebond interpole
(Bouzidi/Yu) est laisse au noyau LBM : ce module fournit seulement q, le point
de paroi et la vitesse de paroi par lien.

Conventions verrouillees :
    Omega_n = {(x,y): |x/a|^n + |y/b|^n < 1}, a=b=1, AR=1, boite fixe.
    couvercle n<infini : arc superieur y>0, u_lid = U_lid * t, t=(n_y,-n_x) -> +x.
    carre n=infini : geometrie distincte (SquareGeometry), couvercle horizontal.

Dependances : numpy ; importe 01_geometry/superellipse.py.
"""

from __future__ import annotations

import os
import sys
from typing import Optional

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_GEOM = os.path.join(os.path.dirname(_HERE), "01_geometry")
if _GEOM not in sys.path:
    sys.path.insert(0, _GEOM)

from superellipse import Superellipse, SquareGeometry  # noqa: E402

EX = np.array([0, 1, 0, -1, 0, 1, -1, -1, 1], dtype=np.float64)
EY = np.array([0, 0, 1, 0, -1, 1, 1, -1, -1], dtype=np.float64)
OPP = np.array([0, 3, 4, 1, 2, 7, 8, 5, 6], dtype=np.int32)

U_REF = 0.05


def build_grid(geom, N: int, half_cell: Optional[bool] = None):
    """Grille (Y,X) sur la boite englobante + masque fluide booleen.

    Y[i,j] = y_i, X[i,j] = x_j (i<->y, j<->x). Renvoie (X, Y, fluid).

    half_cell : decale la grille d'une demi-cellule pour que les parois planes
        tombent a mi-chemin entre noeuds (q = 0.5, rebond demi-pas exact).
        Auto-active pour SquareGeometry (parois uniquement planes), desactive
        sinon (les parois courbes ont deja un q fractionnaire exact). Utile pour
        reproduire le benchmark de Ghia avec le solveur courbe.
    """
    if N < 4:
        raise ValueError("N >= 4 requis.")
    if half_cell is None:
        half_cell = isinstance(geom, SquareGeometry)
    if half_cell:
        dx = 2.0 * geom.a / N
        dy = 2.0 * geom.b / N
        x = -geom.a + 0.5 * dx + np.arange(N) * dx
        y = -geom.b + 0.5 * dy + np.arange(N) * dy
    else:
        x = np.linspace(-geom.a, geom.a, N)
        y = np.linspace(-geom.b, geom.b, N)
    Y, X = np.meshgrid(y, x, indexing="ij")
    fluid = np.asarray(geom.implicit(X, Y) < 0.0)
    return X, Y, fluid


def _is_on_lid(geom, x: float, y: float, mode: str, param: Optional[float]) -> bool:
    if isinstance(geom, SquareGeometry):
        return (y >= geom.b - 1e-12) and (abs(x) <= geom.a + 1e-12)
    m = geom.lid_mask(np.array([[x, y]], dtype=np.float64), mode=mode, param=param)
    return bool(m[0])


def _wall_velocity(geom, wx: float, wy: float, U_lid: float):
    if isinstance(geom, SquareGeometry):
        nx, ny = 0.0, 1.0
    else:
        nx, ny = geom.normal(np.array([wx]), np.array([wy]))
        nx, ny = float(nx[0]), float(ny[0])
    # tangente de reference t = (n_y, -n_x) (aucun flip point-par-point : sur
    # l'arc superieur n_y>0 donc t_x>0 ; sur toute la frontiere cela donne une
    # rotation rigide coherente).
    tx, ty = ny, -nx
    return U_lid * tx, U_lid * ty, nx, ny


def _bisect_wall_t(geom, xf, yf, xn, yn, iters: int) -> float:
    """Fraction t in [0,1] du point de paroi sur le segment fluide->solide."""
    f0 = float(geom.implicit(xf, yf))
    f1 = float(geom.implicit(xn, yn))
    if f0 >= 0.0:
        return 0.0
    if f1 <= 0.0:
        return 1.0
    lo, hi = 0.0, 1.0
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        if float(geom.implicit(xf + mid * (xn - xf), yf + mid * (yn - yf))) < 0.0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def build_links(
    geom,
    X: np.ndarray,
    Y: np.ndarray,
    fluid: np.ndarray,
    U_lid: float = U_REF,
    lid_mode: str = "upper",
    lid_param: Optional[float] = None,
    iters: int = 60,
):
    """Liens fluide->solide avec fraction de paroi exacte et vitesse de paroi."""
    Ni, Nj = fluid.shape
    # pas de grille deduit des coordonnees (gere grilles standard et demi-cellule)
    dx = float(X[0, 1] - X[0, 0]) if Nj > 1 else (2.0 * geom.a) / max(Ni - 1, 1)
    dy = float(Y[1, 0] - Y[0, 0]) if Ni > 1 else (2.0 * geom.b) / max(Nj - 1, 1)

    li, lj, lq, lfrac = [], [], [], []
    lwx, lwy, lux, luy, llid = [], [], [], [], []

    for i in range(Ni):
        for j in range(Nj):
            if not fluid[i, j]:
                continue
            xf, yf = X[i, j], Y[i, j]
            for q in range(1, 9):
                ni = i + int(EY[q])
                nj = j + int(EX[q])
                outside = (ni < 0) or (ni >= Ni) or (nj < 0) or (nj >= Nj)
                if not outside and fluid[ni, nj]:
                    continue
                if outside:
                    xn = xf + EX[q] * dx
                    yn = yf + EY[q] * dy
                else:
                    xn, yn = X[ni, nj], Y[ni, nj]
                t = _bisect_wall_t(geom, xf, yf, xn, yn, iters)
                wx = xf + t * (xn - xf)
                wy = yf + t * (yn - yf)
                on_lid = _is_on_lid(geom, wx, wy, lid_mode, lid_param)
                if on_lid:
                    ux, uy, _, _ = _wall_velocity(geom, wx, wy, U_lid)
                else:
                    ux, uy = 0.0, 0.0
                li.append(i); lj.append(j); lq.append(q); lfrac.append(t)
                lwx.append(wx); lwy.append(wy)
                lux.append(ux); luy.append(uy); llid.append(on_lid)

    return {
        "i": np.asarray(li, dtype=np.int32),
        "j": np.asarray(lj, dtype=np.int32),
        "q": np.asarray(lq, dtype=np.int32),
        "frac": np.asarray(lfrac, dtype=np.float64),
        "wx": np.asarray(lwx, dtype=np.float64),
        "wy": np.asarray(lwy, dtype=np.float64),
        "ux": np.asarray(lux, dtype=np.float64),
        "uy": np.asarray(luy, dtype=np.float64),
        "on_lid": np.asarray(llid, dtype=bool),
        "N": Ni,
        "dx": dx,
        "dy": dy,
        "n": (np.inf if isinstance(geom, SquareGeometry) else geom.n),
        "a": geom.a,
        "b": geom.b,
        "n_fluid": int(fluid.sum()),
        "n_solid": int(fluid.size - fluid.sum()),
        "n_links": len(li),
        "n_lid": int(np.sum(llid)),
    }


class SuperellipseGrid:
    """Masque + liens de paroi pour une geometrie donnee a une resolution N."""

    def __init__(self, geom, N: int, U_lid: float = U_REF,
                 lid_mode: str = "upper", lid_param: Optional[float] = None):
        self.geom = geom
        self.N = N
        self.U_lid = U_lid
        self.lid_mode = lid_mode
        self.lid_param = lid_param
        self.X, self.Y, self.fluid = build_grid(geom, N)
        links = build_links(geom, self.X, self.Y, self.fluid, U_lid=U_lid,
                            lid_mode=lid_mode, lid_param=lid_param)
        self.links = links
        self.area_est = links["n_fluid"] * links["dx"] * links["dy"]
        self.area_exact = geom.area()

    # ------------------------------------------------------------------
    def area_error(self) -> float:
        return abs(self.area_est - self.area_exact) / self.area_exact

    def is_fluid(self, i: int, j: int) -> bool:
        return bool(self.fluid[i, j])

    def summary(self) -> str:
        L = self.links
        return (
            f"n={L['n']:g} N={self.N}  fluid={L['n_fluid']} solid={L['n_solid']}  "
            f"links={L['n_links']} (lid={L['n_lid']})  "
            f"A_est={self.area_est:.6f} A_exact={self.area_exact:.6f} "
            f"err={100.0 * self.area_error():.3f}%"
        )

    # ------------------------------------------------------------------
    def save_npz(self, path: str) -> None:
        L = self.links
        np.savez_compressed(
            path,
            fluid=self.fluid, X=self.X, Y=self.Y,
            i=L["i"], j=L["j"], q=L["q"], frac=L["frac"],
            wx=L["wx"], wy=L["wy"], ux=L["ux"], uy=L["uy"],
            on_lid=L["on_lid"], N=self.N, dx=L["dx"], dy=L["dy"],
            n=L["n"], a=L["a"], b=L["b"],
        )


# ----------------------------------------------------------------------
#  Auto-test
# ----------------------------------------------------------------------
def _symmetry_error(grid: SuperellipseGrid) -> float:
    """Ecart max de symetrie du masque sous x->-x et y->-y (geometries sym.)."""
    f = grid.fluid
    ex = np.abs(f.astype(np.int8) - f[:, ::-1].astype(np.int8)).max()
    ey = np.abs(f.astype(np.int8) - f[::-1, :].astype(np.int8)).max()
    return float(max(ex, ey))


def _wall_symmetry_error(grid: SuperellipseGrid) -> float:
    """Ecart max (en pas de grille) des points de paroi sous x->-x."""
    wx, wy = grid.links["wx"], grid.links["wy"]
    if wx.size == 0:
        return 0.0
    pts = np.column_stack([wx, wy])
    ref = pts.copy()
    ref[:, 0] = -ref[:, 0]
    err = 0.0
    for p in ref:
        d = np.hypot(pts[:, 0] - p[0], pts[:, 1] - p[1]).min()
        err = max(err, d)
    return err / grid.links["dx"]


def main() -> int:
    import argparse

    parser = argparse.ArgumentParser(description="Grille LBM superellipse (masque + liens)")
    parser.add_argument("--n", type=float, nargs="*", default=[1.25, 2.0, 3.0, 8.0])
    parser.add_argument("--ood", type=float, nargs="*", default=[1.0, 0.9, 0.75, 0.6])
    parser.add_argument("--N", type=int, nargs="*", default=[64, 128, 256])
    parser.add_argument("--U", type=float, default=U_REF)
    parser.add_argument("--plot", action="store_true")
    args = parser.parse_args()

    print("=== Convergence de l'aire (masque) ===")
    for n in args.n:
        g = Superellipse(n=n)
        line = [f"n={n:>5g}"]
        for N in args.N:
            gr = SuperellipseGrid(g, N, U_lid=args.U)
            line.append(f"N={N}:{100 * gr.area_error():.3f}%")
        print("  " + "  ".join(line))

    print("\n=== Fractions de paroi q et symetrie ===")
    for n in args.n:
        g = Superellipse(n=n)
        for N in (128, 256):
            gr = SuperellipseGrid(g, N, U_lid=args.U)
            q = gr.links["frac"]
            sym_f = _symmetry_error(gr)
            sym_w = _wall_symmetry_error(gr)
            print(
                f"  n={n:>5g} N={N:>3d}  q in [{q.min():.4f},{q.max():.4f}]  "
                f"sym_mask={sym_f:.0f}  sym_wall={sym_w:.3e} pas"
            )
            assert q.size > 0, "aucun lien de paroi"
            assert q.min() > 0.0 and q.max() <= 1.0, "q hors de (0,1]"
            assert sym_f == 0.0, "masque non symetrique"

    print("\n=== Couvercle tangent (impermeabilite) ===")
    for n in args.n:
        g = Superellipse(n=n)
        gr = SuperellipseGrid(g, 256, U_lid=args.U)
        L = gr.links
        lid = L["on_lid"]
        assert lid.sum() > 0, "aucun lien de couvercle detecte"
        nx, ny = g.normal(L["wx"][lid], L["wy"][lid])
        un = L["ux"][lid] * nx + L["uy"][lid] * ny
        speed = np.hypot(L["ux"][lid], L["uy"][lid])
        print(
            f"  n={n:>5g}  liens couvercle={int(lid.sum()):>5d}  "
            f"|u.n|max={np.abs(un).max():.2e}  |u| in "
            f"[{speed.min():.4f},{speed.max():.4f}]  "
            f"mean(u_x)={L['ux'][lid].mean():+.5f}"
        )
        assert np.abs(un).max() < 1e-9
        assert L["ux"][lid].mean() > 0.0

    print("\n=== OOD (non convexe, cuspides) : robustesse ===")
    for n in args.ood:
        g = Superellipse(n=n)
        gr = SuperellipseGrid(g, 256, U_lid=args.U)
        q = gr.links["frac"]
        ok = np.isfinite(q).all() and np.isfinite(gr.links["wx"]).all()
        print(f"  n={n:>5g}  links={gr.links['n_links']:>6d}  "
              f"q in [{q.min():.4f},{q.max():.4f}]  finite={ok}")
        assert ok

    print("\n=== Carre (geometrie distincte) ===")
    sq = SquareGeometry()
    gr = SuperellipseGrid(sq, 256, U_lid=args.U)
    L = gr.links
    print(f"  square  fluid={L['n_fluid']}/{gr.fluid.size}  "
          f"links={L['n_links']}  lid={L['n_lid']}  err_aire={100*gr.area_error():.3f}%")

    if args.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(1, len(args.n), figsize=(4 * len(args.n), 4.2))
        if len(args.n) == 1:
            axes = [axes]
        outdir = os.path.join(os.path.dirname(_HERE), "08_figures")
        os.makedirs(outdir, exist_ok=True)
        for ax, n in zip(axes, args.n):
            g = Superellipse(n=n)
            gr = SuperellipseGrid(g, 128, U_lid=args.U)
            ax.imshow(gr.fluid.T, origin="lower", cmap="Blues", extent=[-g.a, g.a, -g.b, g.b])
            L = gr.links
            ax.scatter(L["wx"][~L["on_lid"]], L["wy"][~L["on_lid"]], s=1, c="k")
            ax.scatter(L["wx"][L["on_lid"]], L["wy"][L["on_lid"]], s=1, c="r")
            ax.set_title(f"n={n:g}, N=128")
            ax.set_aspect("equal")
        fpath = os.path.join(outdir, "grid_mask_links.png")
        fig.tight_layout()
        fig.savefig(fpath, dpi=150)
        print(f"\nFigure sauvee : {fpath}")

    print("\nOK - tous les controles passes.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
