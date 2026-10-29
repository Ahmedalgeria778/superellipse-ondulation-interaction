# -*- coding: utf-8 -*-
"""
solver/cases.py
Liste des cas de campagne. Charge geometry/configs.py par chemin de fichier
(pour eviter le conflit de nom avec provenance/.../geometry.py du Projet A).
"""

from __future__ import annotations

import importlib.util
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod  # requis par dataclasses
    spec.loader.exec_module(mod)
    return mod


_cfg = _load("geo_configs", os.path.join(_ROOT, "geometry", "configs.py"))
CaseConfig = _cfg.CaseConfig
PILOT_4CASES = _cfg.PILOT_4CASES


def pilot_cases():
    return list(PILOT_4CASES)


def by_name(name):
    for c in PILOT_4CASES:
        if c.name == name or c.case_id == name:
            return c
    raise KeyError(name)
