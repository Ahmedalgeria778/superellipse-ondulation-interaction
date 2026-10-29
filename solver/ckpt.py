# -*- coding: utf-8 -*-
"""
solver/ckpt.py
==============
Checkpoint et reprise d'un run LBM **sans toucher a la provenance**.

Principe
--------
`WavyMRT.run()` resout `mrt_step_wavy` comme un nom GLOBAL du module
`lbm_wavy` : on installe donc un hook sur `lbm_wavy.mrt_step_wavy` (simple
fonction Python qui rappelle le noyau Numba d'origine). Toutes les `every`
iterations, l'etat `f` est ecrit de facon atomique dans
`results/ckpt/<case_id>.npz`.

A l'execution d'un run, `load()` recharge `f` s'il existe et si le case_id
correspond ; `run_case.run_one()` le passe alors via `solver.f_init` :
`WavyMRT.run()` repart de cet etat (la rampe du lid se refait, soit ~5 t*
= quelques centaines d'iterations, puis la boucle continue).

Le critere de convergence est un critere sur `dU` (et non sur le compteur
d'iterations) : la reprise ne modifie donc ni la solution ni le verdict, elle
ne fait que rebattre les cartes de la rampe. Le fichier est supprime lorsque
le cas est `converged`.

L'ecriture est atomique (`.npz` temporaire + `os.replace`) : un processus
tue pendant l'ecriture ne laisse jamais un fichier corrompu.
"""

from __future__ import annotations

import os

import numpy as np

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CKPT_DIR = os.path.join(_ROOT, "results", "ckpt")

_state = {"it": 0, "cid": None, "path": None, "every": 0, "active": False}


def install(every=200000):
    """Installe (une seule fois) le hook sur le noyau de la provenance."""
    import lbm_wavy as lw
    _state["every"] = int(every)
    if getattr(lw, "_ckpt_hook", None) is not None:
        return
    orig = lw.mrt_step_wavy

    def wrapper(f, *args, **kwargs):
        out = orig(f, *args, **kwargs)
        st = _state
        if st["active"]:
            st["it"] += 1
            if st["every"] > 0 and st["it"] % st["every"] == 0:
                save(out[0], st["cid"], st["path"], st["it"])
        return out

    wrapper.__doc__ = getattr(orig, "__doc__", None)
    lw.mrt_step_wavy = wrapper
    lw._ckpt_hook = wrapper


def activate(cid, path):
    _state["cid"] = cid
    _state["path"] = path
    _state["it"] = 0
    _state["active"] = True


def deactivate():
    _state["active"] = False


def save(f, cid, path, it):
    try:
        os.makedirs(os.path.dirname(path), exist_ok=True)
        tmp = path + ".tmp.npz"
        np.savez_compressed(tmp, f=f, it=np.int64(it), cid=np.str_(cid))
        os.replace(tmp, path)
    except Exception as exc:          # un echec de checkpoint ne doit pas
        print("[ckpt  ] ecriture impossible: %s" % exc, flush=True)  # tuer le run


def load(path, cid):
    """Renvoie f du checkpoint s'il existe et correspond, sinon None."""
    if not os.path.exists(path):
        return None
    try:
        with np.load(path) as d:
            if str(d["cid"]) != str(cid):
                print("[ckpt  ] case_id different -> demarrage a zero", flush=True)
                return None
            f = np.ascontiguousarray(d["f"])
            it = int(d["it"])
    except Exception as exc:
        print("[ckpt  ] fichier illisible (%s) -> demarrage a zero" % exc,
              flush=True)
        return None
    print("[ckpt  ] reprise depuis %s (it=%d)" % (os.path.basename(path), it),
          flush=True)
    return f


def clear(path):
    try:
        if os.path.exists(path):
            os.remove(path)
    except OSError:
        pass
