# -*- coding: utf-8 -*-
"""
01_geometry/superellipse.py

Geometrie superelliptique pour l'etude LBM / operateur physique.

Domaine :
    Omega_n = { (x,y) : |x/a|^n + |y/b|^n < 1 }

Choix par defaut (article principal) :
    a = b = 1  ->  AR = b/a = 1,  L_c = 2a = 2 (boite englobante fixe)
    n = 2      ->  cercle
    n = 1      ->  losange
    n > 1      ->  convexe, tend vers le carre quand n -> infini
    n < 1      ->  non convexe (pointes/cusps aux points cardinaux), OOD

Aire exacte :
    A(n) = 4ab * Gamma(1 + 1/n)^2 / Gamma(1 + 2/n)

=== Convention de reference verrouillee (couvercle) ===

Pour n < infini (superellipses) :
    couvercle = tout l'arc superieur  Gamma_lid = { (x,y) in dOmega_n : y >= 0 }
    (mode='upper', xfrac=1)
    u_lid = U_lid * t(x,y),  t = (n_y, -n_x) orientee vers +x
    -> |u_lid| = U_lid et u_lid . n = 0 partout.

Pour n = infini (carre, benchmark Ghia) :
    couvercle = arete superieure horizontale, u = (U_lid, 0) (flat-horizontal).
    Voir la fonction square_lid(). Ne pas extrapoler n -> infini pour
    reproduire les coins : le carre est traite comme une geometrie distincte.

Modes de selection du couvercle : 'upper' (reference), 'xfrac', 'yfrac',
'normal', 'normal45'. Le mode 'flat' (selection) est distinct de la
convention 'flat-horizontal' (vitesse) reservee au carre.

Fournit : equation implicite, gradient/normale/tangente, aire, echantillonnage
de la frontiere (uniforme en abscisse curviligne), masque interieur, SDF
(approchee et numerique), et definition du couvercle tangent.

Dependances : numpy (matplotlib optionnel pour le trace).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Callable, Optional, Tuple

import numpy as np

# ----------------------------------------------------------------------
#  Grilles recommandees (cf. plan de campagne)
# ----------------------------------------------------------------------
N_GRID_MAIN = (1.25, 1.5, 1.75, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0)
N_GRID_COARSE = (1.25, 1.5, 2.0, 3.0, 4.0, 8.0)
N_GRID_OOD = (1.0, 0.9, 0.75, 0.6)

RE_GRID = (100.0, 400.0, 1000.0)
RE_GRID_FINE = (100.0, 200.0, 400.0, 600.0, 800.0, 1000.0)

EPS = 1e-12


def _spow(t: np.ndarray, p: float) -> np.ndarray:
    """Puissance signee : sign(t) * |t|^p."""
    return np.sign(t) * np.power(np.abs(t), p)


@dataclass
class Superellipse:
    """Geometrie Omega_n = {(x,y): |x/a|^n + |y/b|^n < 1}."""

    a: float = 1.0
    b: float = 1.0
    n: float = 2.0
    _cache: dict = field(default_factory=dict, repr=False, compare=False)

    def __post_init__(self) -> None:
        if self.a <= 0.0 or self.b <= 0.0:
            raise ValueError("a et b doivent etre strictement positifs.")
        if self.n <= 0.0:
            raise ValueError("n doit etre strictement positif.")

    # ------------------------------------------------------------------
    #  Grandeurs globales
    # ------------------------------------------------------------------
    @property
    def aspect_ratio(self) -> float:
        return self.b / self.a

    @property
    def characteristic_length(self) -> float:
        """L_c = 2a, identique pour toutes les geometries (boite englobante)."""
        return 2.0 * self.a

    @property
    def bounding_box(self) -> Tuple[float, float, float, float]:
        return (-self.a, self.a, -self.b, self.b)

    def area(self) -> float:
        """Aire exacte A(n) = 4ab Gamma(1+1/n)^2 / Gamma(1+2/n)."""
        return (
            4.0
            * self.a
            * self.b
            * math.gamma(1.0 + 1.0 / self.n) ** 2
            / math.gamma(1.0 + 2.0 / self.n)
        )

    def area_numeric(self, num: int = 20001) -> float:
        """Aire par integration numerique de la frontiere (controle)."""
        x, y = self.boundary_points(num).T
        x = np.append(x, x[0])
        y = np.append(y, y[0])
        return 0.5 * abs(np.sum(x[:-1] * y[1:] - x[1:] * y[:-1]))

    def equivalent_diameter(self) -> float:
        """Diametre hydraulique equivalent D = 2 sqrt(A/pi) (indicatif)."""
        return 2.0 * math.sqrt(self.area() / math.pi)

    # ------------------------------------------------------------------
    #  Equation implicite et champs
    # ------------------------------------------------------------------
    def implicit(self, x, y) -> np.ndarray:
        """F(x,y) = |x/a|^n + |y/b|^n - 1 (< 0 interieur, > 0 exterieur)."""
        x = np.asarray(x, dtype=np.float64)
        y = np.asarray(y, dtype=np.float64)
        return np.power(np.abs(x) / self.a, self.n) + np.power(
            np.abs(y) / self.b, self.n
        ) - 1.0

    def is_inside(self, x, y, tol: float = 0.0) -> np.ndarray:
        return self.implicit(x, y) <= tol

    def gradient(self, x, y) -> Tuple[np.ndarray, np.ndarray]:
        """Gradient de F. Tolerant aux cusps (n < 1) via garde sur |x|,|y|."""
        x = np.asarray(x, dtype=np.float64)
        y = np.asarray(y, dtype=np.float64)
        n, a, b = self.n, self.a, self.b
        xa = np.abs(x) / a
        ya = np.abs(y) / b
        xa_safe = np.where(xa > EPS, xa, EPS)
        ya_safe = np.where(ya > EPS, ya, EPS)
        gx = (n / a) * np.sign(x) * np.power(xa_safe, n - 1.0)
        gy = (n / b) * np.sign(y) * np.power(ya_safe, n - 1.0)
        return gx, gy

    def normal(self, x, y) -> Tuple[np.ndarray, np.ndarray]:
        """Normale exterieure unitaire n = grad F / |grad F|."""
        gx, gy = self.gradient(x, y)
        norm = np.hypot(gx, gy)
        norm = np.where(norm > EPS, norm, np.inf)
        return gx / norm, gy / norm

    def tangent(self, x, y) -> Tuple[np.ndarray, np.ndarray]:
        """Tangente unitaire t = (n_y, -n_x) (couvercle : orientation gauche -> droite)."""
        nx, ny = self.normal(x, y)
        return ny, -nx

    # ------------------------------------------------------------------
    #  Frontiere
    # ------------------------------------------------------------------
    def boundary_points(self, num: int = 4096, arc_uniform: bool = True) -> np.ndarray:
        """Points sur la frontiere, uniformement repartis en abscisse curviligne."""
        if num < 3:
            raise ValueError("num >= 3 requis.")
        if not arc_uniform:
            t = np.linspace(0.0, 2.0 * np.pi, num, endpoint=False)
            return self._param_curve(t)
        dense = max(8 * num, 8192)
        t = np.linspace(0.0, 2.0 * np.pi, dense, endpoint=False)
        pts = self._param_curve(t)
        seg = np.hypot(np.diff(pts[:, 0]), np.diff(pts[:, 1]))
        s = np.concatenate(([0.0], np.cumsum(seg)))
        total = s[-1]
        target = np.linspace(0.0, total, num + 1)[:-1]
        closed = np.vstack([pts, pts[0]])
        s_closed = np.concatenate((s, [total]))
        x = np.interp(target, s_closed, closed[:, 0])
        y = np.interp(target, s_closed, closed[:, 1])
        return np.column_stack([x, y])

    def _param_curve(self, t: np.ndarray) -> np.ndarray:
        """Forme parametrique exacte : sign(c)|c|^(2/n), sign(s)|s|^(2/n)."""
        p = 2.0 / self.n
        x = self.a * _spow(np.cos(t), p)
        y = self.b * _spow(np.sin(t), p)
        return np.column_stack([x, y])

    def boundary_normals(self, num: int = 4096) -> Tuple[np.ndarray, np.ndarray]:
        pts = self.boundary_points(num)
        nx, ny = self.normal(pts[:, 0], pts[:, 1])
        return pts, np.column_stack([nx, ny])

    # ------------------------------------------------------------------
    #  Couvercle mobile (tangent)
    # ------------------------------------------------------------------
    def lid_mask(
        self,
        pts: np.ndarray,
        mode: str = "upper",
        param: Optional[float] = None,
    ) -> np.ndarray:
        """
        Selection des points du couvercle mobile parmi des points de frontiere.

        mode = 'upper'   : y > 0 (arc superieur complet)          [reference]
               'all'     : toute la frontiere (couvercle tournant /
                           validation de Couette circulaire rigide,
                           u = U_lid * t partout)
               'xfrac'   : |x| <= param * a
               'yfrac'   : y   >= param * b
               'normal'  : angle de la normale a la verticale <= param (rad),
                           equivalent a |n_x| <= tan(param) |n_y|
               'normal45': |n_x| <= |n_y| (convenience, param ignore)
        """
        x, y = pts[:, 0], pts[:, 1]
        if mode == "all":
            return np.ones(x.shape, dtype=bool)
        if mode == "upper":
            return y > 0.0
        if mode == "xfrac":
            return (y > 0.0) & (np.abs(x) <= (1.0 if param is None else param) * self.a)
        if mode == "yfrac":
            return y >= (0.0 if param is None else param) * self.b
        if mode == "normal":
            nx, ny = self.normal(x, y)
            theta = 0.0 if param is None else param
            return (y > 0.0) & (np.abs(nx) <= abs(math.tan(theta)) * np.abs(ny))
        if mode == "normal45":
            nx, ny = self.normal(x, y)
            return (y > 0.0) & (np.abs(nx) <= np.abs(ny))
        raise ValueError(f"mode couvercle inconnu : {mode!r}")

    def lid(
        self,
        num: int = 4096,
        U_lid: float = 1.0,
        mode: str = "upper",
        param: Optional[float] = None,
        modulation: Optional[Callable[[np.ndarray], np.ndarray]] = None,
        tangential: bool = True,
    ) -> dict:
        """
        Couvercle : u_lid = U_lid * modulation(x) * direction imposee.

        tangential=True  (reference pour n < infini) : direction = t = (n_y,-n_x),
            orientee vers +x, d'ou u_lid . n = 0 partout.
        tangential=False (convention flat-horizontal, reservee au carre n=infini) :
            direction = (1, 0) ; utiliser square_lid() pour le carre exact.
        """
        allpts = self.boundary_points(num)
        mask = self.lid_mask(allpts, mode=mode, param=param)
        pts = allpts[mask]
        nx, ny = self.normal(pts[:, 0], pts[:, 1])
        scale = np.ones(pts.shape[0])
        if modulation is not None:
            scale = np.asarray(modulation(pts[:, 0]), dtype=np.float64)
        speed = U_lid * scale
        if tangential:
            tx, ty = ny, -nx
            vel = np.column_stack([speed * tx, speed * ty])
            if pts.shape[0] and np.mean(vel[:, 0]) < 0.0:
                tx, ty = -tx, -ty
                vel = np.column_stack([speed * tx, speed * ty])
        else:
            tx = np.ones_like(nx)
            ty = np.zeros_like(ny)
            vel = np.column_stack([speed, np.zeros_like(speed)])
        return {
            "points": pts,
            "normals": np.column_stack([nx, ny]),
            "tangents": np.column_stack([tx, ty]),
            "velocity": vel,
            "speed": speed,
            "mode": mode,
            "param": param,
            "tangential": tangential,
        }

    # ------------------------------------------------------------------
    #  SDF
    # ------------------------------------------------------------------
    def sdf_approx(self, x, y) -> np.ndarray:
        """SDF au premier ordre : F / |grad F| (signe : negatif a l'interieur)."""
        F = self.implicit(x, y)
        gx, gy = self.gradient(x, y)
        gnorm = np.hypot(gx, gy)
        gnorm = np.where(gnorm > EPS, gnorm, np.inf)
        return F / gnorm

    def nearest_boundary(
        self, x, y, boundary: Optional[np.ndarray] = None, chunk: int = 1024
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Distance signee et point de frontiere le plus proche (recherche exacte)."""
        if boundary is None:
            boundary = self.boundary_points(4096)
        X = np.asarray(x, dtype=np.float64)
        Y = np.asarray(y, dtype=np.float64)
        shape = X.shape
        xf, yf = X.ravel(), Y.ravel()
        bx, by = boundary[:, 0], boundary[:, 1]
        dist = np.empty(xf.size)
        px = np.empty(xf.size)
        py = np.empty(xf.size)
        for start in range(0, xf.size, chunk):
            sl = slice(start, start + chunk)
            dx = xf[sl, None] - bx[None, :]
            dy = yf[sl, None] - by[None, :]
            d2 = dx * dx + dy * dy
            idx = np.argmin(d2, axis=1)
            rows = np.arange(idx.size)
            dist[sl] = np.sqrt(d2[rows, idx])
            px[sl] = bx[idx]
            py[sl] = by[idx]
        sign = np.where(self.implicit(X, Y) <= 0.0, -1.0, 1.0).ravel()
        return (
            (sign * dist).reshape(shape),
            px.reshape(shape),
            py.reshape(shape),
        )

    def sdf(self, x, y, boundary: Optional[np.ndarray] = None) -> np.ndarray:
        """SDF numerique signee (distances exactes a la frontiere echantillonnee)."""
        d, _, _ = self.nearest_boundary(x, y, boundary=boundary)
        return d

    # ------------------------------------------------------------------
    #  Grille
    # ------------------------------------------------------------------
    def grid(
        self, nx: int = 256, ny: Optional[int] = None, pad: float = 1.0
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Grille cartesienne sur la boite englobante ; masque interieur booleen."""
        ny = nx if ny is None else ny
        x = np.linspace(-self.a * pad, self.a * pad, nx)
        y = np.linspace(-self.b * pad, self.b * pad, ny)
        X, Y = np.meshgrid(x, y, indexing="ij")
        return X, Y, self.is_inside(X, Y)

    # ------------------------------------------------------------------
    #  Metadonnees (cf. fichier par simulation)
    # ------------------------------------------------------------------
    def metadata(self) -> dict:
        return {
            "geometry": "superellipse",
            "n": self.n,
            "a": self.a,
            "b": self.b,
            "aspect_ratio": self.aspect_ratio,
            "characteristic_length": self.characteristic_length,
            "area": self.area(),
            "convex": self.n >= 1.0,
            "regular_boundary": self.n >= 1.0,
        }

    def summary(self) -> str:
        return (
            f"Superellipse n={self.n:g}  a={self.a:g}  b={self.b:g}  "
            f"AR={self.aspect_ratio:g}  L_c={self.characteristic_length:g}  "
            f"A={self.area():.6f}  convex={self.n >= 1.0}"
        )

    # ------------------------------------------------------------------
    #  Trace (optionnel)
    # ------------------------------------------------------------------
    def plot(self, num: int = 1024, ax=None, **kwargs):
        import matplotlib.pyplot as plt

        if ax is None:
            _, ax = plt.subplots(figsize=(4, 4))
        pts = self.boundary_points(num)
        ax.plot(pts[:, 0], pts[:, 1], **kwargs)
        ax.set_aspect("equal")
        return ax


def superellipse_points(n: float, num: int = 4096, a: float = 1.0, b: float = 1.0):
    """Raccourci : points de frontiere pour un n donne."""
    return Superellipse(a=a, b=b, n=n).boundary_points(num)


def square_lid(
    a: float = 1.0,
    b: float = 1.0,
    num: int = 4096,
    U_lid: float = 1.0,
    modulation: Optional[Callable[[np.ndarray], np.ndarray]] = None,
) -> dict:
    """
    Convention flat-horizontal pour le carre (n = infini, benchmark Ghia).

    Arete superieure y = b, x in [-a, a] ; u = (U_lid*modulation(x), 0).
    Geometrie distincte : on ne passe pas par la limite numerique n -> infini.
    """
    x = np.linspace(-a, a, num)
    y = np.full_like(x, b)
    scale = np.ones_like(x)
    if modulation is not None:
        scale = np.asarray(modulation(x), dtype=np.float64)
    vel = np.column_stack([U_lid * scale, np.zeros_like(x)])
    return {
        "points": np.column_stack([x, y]),
        "normals": np.tile([0.0, 1.0], (num, 1)),
        "tangents": np.tile([1.0, 0.0], (num, 1)),
        "velocity": vel,
        "speed": U_lid * scale,
        "mode": "flat-horizontal",
        "param": None,
        "tangential": False,
    }


def square_geometry(a: float = 1.0, b: float = 1.0) -> "SquareGeometry":
    """Retourne la geometrie carree explicite (n = infini)."""
    return SquareGeometry(a=a, b=b)


class SquareGeometry:
    """Carre [-a,a] x [-b,b] comme geometrie distincte (benchmark Ghia)."""

    def __init__(self, a: float = 1.0, b: float = 1.0):
        self.a = a
        self.b = b

    @property
    def aspect_ratio(self) -> float:
        return self.b / self.a

    @property
    def characteristic_length(self) -> float:
        return 2.0 * self.a

    def area(self) -> float:
        return 4.0 * self.a * self.b

    def implicit(self, x, y) -> np.ndarray:
        x = np.asarray(x, dtype=np.float64)
        y = np.asarray(y, dtype=np.float64)
        return np.maximum(np.abs(x) / self.a, np.abs(y) / self.b) - 1.0

    def is_inside(self, x, y, tol: float = 0.0) -> np.ndarray:
        return self.implicit(x, y) <= tol

    def lid(self, num: int = 4096, U_lid: float = 1.0, **kwargs) -> dict:
        return square_lid(
            a=self.a, b=self.b, num=num, U_lid=U_lid,
            modulation=kwargs.get("modulation"),
        )

    def metadata(self) -> dict:
        return {
            "geometry": "square",
            "n": float("inf"),
            "a": self.a,
            "b": self.b,
            "aspect_ratio": self.aspect_ratio,
            "characteristic_length": self.characteristic_length,
            "area": self.area(),
            "convex": True,
            "regular_boundary": True,
        }

    def summary(self) -> str:
        return (
            f"Square a={self.a:g} b={self.b:g} AR={self.aspect_ratio:g} "
            f"L_c={self.characteristic_length:g} A={self.area():.6f}"
        )


def geometry_table(n_values=N_GRID_MAIN, a: float = 1.0, b: float = 1.0) -> str:
    header = f"{'n':>6} {'A(n)':>10} {'AR':>5} {'L_c':>5} {'convex':>7}"
    lines = [header, "-" * len(header)]
    for n in n_values:
        g = Superellipse(a=a, b=b, n=n)
        lines.append(
            f"{g.n:>6g} {g.area():>10.6f} {g.aspect_ratio:>5g} "
            f"{g.characteristic_length:>5g} {str(g.n >= 1.0):>7}"
        )
    return "\n".join(lines)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Superellipse geometry module")
    parser.add_argument("--n", type=float, nargs="*", default=None)
    parser.add_argument("--a", type=float, default=1.0)
    parser.add_argument("--b", type=float, default=1.0)
    parser.add_argument("--num", type=int, default=4096)
    parser.add_argument("--plot", action="store_true")
    args = parser.parse_args()

    ns = args.n if args.n else list(N_GRID_MAIN)
    print("=== Grille principale (convexe) ===")
    print(geometry_table(ns, a=args.a, b=args.b))
    print("\n=== OOD (non convexe) ===")
    print(geometry_table(N_GRID_OOD, a=args.a, b=args.b))

    print("\n=== Controle aire (analytique vs numerique) ===")
    for n in ns:
        g = Superellipse(a=args.a, b=args.b, n=n)
        print(f"n={n:>5g}  A_exact={g.area():.6f}  A_num={g.area_numeric():.6f}")

    g = Superellipse(a=args.a, b=args.b, n=2.0)
    print("\n=== Couvercle tangent (mode=upper, n=2, cercle) ===")
    lid = g.lid(num=args.num, U_lid=0.05, mode="upper")
    pts, nrm, vel = lid["points"], lid["normals"], lid["velocity"]
    un = np.sum(vel * nrm, axis=1)
    print(f"points couvercle = {pts.shape[0]}")
    print(f"y min/max = {pts[:,1].min():.4f} / {pts[:,1].max():.4f}")
    print(f"|u| min/max = {np.linalg.norm(vel,axis=1).min():.4f} / "
          f"{np.linalg.norm(vel,axis=1).max():.4f}")
    print(f"vx min/max = {vel[:,0].min():.4f} / {vel[:,0].max():.4f} (attendu >= 0)")
    print(f"|u.n| max = {np.abs(un).max():.3e} (attendu ~ 0, impermeabilite)")

    print("\n=== Carre flat-horizontal (n=infini, Ghia) ===")
    sq = square_lid(a=args.a, b=args.b, num=args.num, U_lid=0.05)
    print(f"points = {sq['points'].shape[0]}  y={sq['points'][0,1]:.2f}  "
          f"|u|={np.linalg.norm(sq['velocity'],axis=1).max():.4f}  "
          f"u=(ux,uy)=({sq['velocity'][0,0]:.3f},{sq['velocity'][0,1]:.3f})")

    if args.plot:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(5, 5))
        for n in ns:
            Superellipse(n=n).plot(num=1024, ax=ax, label=f"n={n:g}")
        for n in N_GRID_OOD:
            Superellipse(n=n).plot(num=1024, ax=ax, ls="--", label=f"n={n:g} (OOD)")
        ax.legend(fontsize=8, ncol=2)
        ax.set_title("Superellipses Omega_n (a=b=1)")
        out = "01_geometry/superellipses_preview.png"
        fig.tight_layout()
        fig.savefig(out, dpi=150)
        print(f"\nFigure sauvee : {out}")
