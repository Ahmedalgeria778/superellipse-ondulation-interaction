# -*- coding: utf-8 -*-
"""
solver/check_health.py
======================
Controle de sante des calculs en cours : un verdict PASS / WARN / FAIL par cas.

Verifie, pour chaque cas LBM vivant :
  - exactement un processus de calcul actif (pas de doublon, pas de zombie) ;
  - consommation CPU reelle sur une fenetre glissante ;
  - log de progression valide et a jour (chemin, derniere iteration, dU) ;
  - decroissance effective de dU ;
  - checkpoint present et son age ;
  - coherence entre le case_id lance et le nom de sortie attendu ;
  - chemins de sortie existants et inscriptibles ;
  - espace disque.

Usage :
  python -m solver.check_health            # controle complet
  python -m solver.check_health --quiet    # seulement les problemes
"""

from __future__ import annotations

import glob
import json
import os
import re
import shutil
import subprocess
import sys
import time

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
LOGS = os.path.join(_ROOT, "results", "logs")
SCAL = os.path.join(_ROOT, "results", "scalars")
CKPT = os.path.join(_ROOT, "results", "ckpt")

WINDOW = 20.0          # s, fenetre de mesure CPU
MIN_CPU_FRACTION = 0.25  # un enfant doit consommer au moins 25 % de la fenetre
STALE_LOG_S = 3600.0   # au-dela, un log muet est suspect
TREND_SPAN = 120000    # it, recul minimal pour juger la tendance de dU
TREND_FACTOR = 3.0     # dU final / dU de reference au-dela duquel on alerte
PAUSE_JOURNAL = os.path.join(_ROOT, "results", "pause_journal.jsonl")
HISTORY = os.path.join(_ROOT, "results", "health_history.csv")

FLAGS = ("n", "alpha", "Nw", "Re", "N", "tol", "max_iter")


def ps_json(query):
    out = subprocess.run(["powershell", "-NoProfile", "-Command", query],
                         capture_output=True, text=True).stdout
    d = json.loads(out.strip() or "[]")
    return [d] if isinstance(d, dict) else d


def snapshot():
    q = ("Get-CimInstance Win32_Process -Filter \"Name='python.exe'\" | "
         "Select-Object ProcessId,UserModeTime,KernelModeTime,WorkingSetSize,"
         "CreationDate,CommandLine | ConvertTo-Json -Compress")
    out = {}
    for p in ps_json(q):
        cmd = " ".join((p.get("CommandLine") or "").split())
        out[int(p["ProcessId"])] = dict(
            cpu=(int(p.get("UserModeTime") or 0)
                 + int(p.get("KernelModeTime") or 0)) / 1e7,
            mem=int(p.get("WorkingSetSize") or 0) / 2 ** 20,
            cmd=cmd, start=p.get("CreationDate"))
    return out


def flags(cmd):
    got = {}
    for f in FLAGS:
        m = re.search(r"--%s\s+(\S+)" % f, cmd)
        if m:
            got[f] = m.group(1)
    got["log"] = None
    m = re.search(r'--log\s+(?:"([^"]+)"|(\S+))', cmd)
    if m:
        got["log"] = (m.group(1) or m.group(2)).rstrip('"')
    return got


def case_id(f):
    n = float(f["n"])
    tag = "se_n%g_w%.2f_Nw%s_Re%d_N%s" % (
        n, float(f["alpha"]), f["Nw"], float(f["Re"]), f["N"])
    tol = float(f.get("tol", 1e-9))
    return tag + ("" if tol == 1e-9 else "_tol%.0e" % tol)


def pause_windows():
    """case_id -> liste d'intervalles [t0, t1] de suspension enregistres."""
    if not os.path.exists(PAUSE_JOURNAL):
        return {}
    out = {}
    open_ = {}
    for line in open(PAUSE_JOURNAL, encoding="utf-8"):
        line = line.strip()
        if not line:
            continue
        try:
            ev = json.loads(line)
        except ValueError:
            continue
        for cid in ev.get("cases", []):
            if ev["event"] == "pause":
                open_[cid] = float(ev["at"])
            elif cid in open_:
                out.setdefault(cid, []).append((open_.pop(cid), float(ev["at"])))
    now = time.time()
    for cid, t0 in open_.items():        # pause encore en cours
        out.setdefault(cid, []).append((t0, now))
    return out


def suspended_overlap(windows, t_from, t_to):
    return sum(max(0.0, min(t1, t_to) - max(t0, t_from)) for t0, t1 in windows)


def trend_verdict(prog):
    """Compare dU final a un point de reference situe au moins TREND_SPAN avant.

    Un transitoire de reprise (rampe rejouee) fait remonter dU sur quelques
    dizaines de milliers d'iterations puis le fait redescendre : on ne conclut
    qu'avec un recul suffisant et un facteur franc.
    """
    if len(prog) < 3:
        return "ok", None
    it_last, du_last = prog[-1][0], prog[-1][1]
    refs = [p[1] for p in prog if p[0] <= it_last - TREND_SPAN]
    if not refs:
        return "jeune", None
    ref = min(refs)
    if du_last > TREND_FACTOR * ref:
        return "remonte", (ref, du_last)
    return "ok", ref


def last_progress(log):
    if not log or not os.path.exists(log):
        return None
    txt = open(log, encoding="utf-8", errors="replace").read()
    pat = r"it=\s*(\d+)\s+dU=([\d.eE+-]+)(?:\s+t=([\d.]+)s)?"
    return [(int(i), float(d), float(t) if t else None)
            for i, d, t in re.findall(pat, txt)]


def main():
    quiet = "--quiet" in sys.argv
    a = snapshot()
    time.sleep(WINDOW)
    b = snapshot()

    cases = {}
    for pid, p in b.items():
        if "solver.run_case" not in p["cmd"]:
            continue
        f = flags(p["cmd"])
        if "n" not in f or "Re" not in f:
            continue
        cid = case_id(f)
        e = cases.setdefault(cid, dict(procs=[], flags=f))
        e["procs"].append((pid, b[pid]["cpu"] - a.get(pid, dict(cpu=b[pid]["cpu"]))["cpu"],
                           b[pid]["mem"]))

    verdicts = []
    rows = []
    pw = pause_windows()
    stamps = {}
    for cid, e in cases.items():
        lp = e["flags"].get("log")
        if lp and os.path.exists(lp):
            stamps[cid] = os.path.getmtime(lp)

    print("fenetre CPU = %.0f s | %d cas LBM detectes\n" % (WINDOW, len(cases)))
    print("rappel : le compteur it= des logs repart de zero a chaque reprise ; "
          "le ckpt porte le compte cumule\n")
    hdr = "%-40s %-6s %s" % ("cas", "etat", "detail")
    print(hdr)
    print("-" * len(hdr))
    for cid in sorted(cases):
        e = cases[cid]
        issues, notes = [], []
        busy = [p for p in e["procs"] if p[1] > WINDOW * MIN_CPU_FRACTION]
        idle = [p for p in e["procs"] if p not in busy]

        if not busy:
            issues.append("aucun processus actif (CPU=0)")
        if len(busy) > 1:
            issues.append("%d processus actifs en parallele" % len(busy))
        if len(idle) > 1:
            issues.append("%d processus dormants (attendu: 1)" % len(idle))

        log = e["flags"].get("log")
        prog = last_progress(log)
        if log is None:
            issues.append("--log absent")
        elif not log.lower().endswith(".log"):
            issues.append("chemin de log tronque (--log non protege)")
        elif os.path.dirname(os.path.abspath(log)) != os.path.abspath(LOGS):
            issues.append("log hors de results/logs")
        if prog:
            it, du, _ = prog[-1]
            notes.append("it=%d dU=%.2e" % (it, du))
            verdict, ref = trend_verdict(prog)
            if verdict == "remonte":
                issues.append("dU remonte (%.2e -> %.2e) sur %d it"
                              % (ref, du, TREND_SPAN))
            elif verdict == "jeune":
                notes.append("tendance dU non jugeable (recul < %d it)" % TREND_SPAN)
            elif ref:
                notes.append("dU / ref = x%.1f" % (du / ref))
            else:
                notes.append("dU decroit")
            stamp = stamps.get(cid)
            if stamp is not None and busy:
                silence = time.time() - stamp
                dt_line = (prog[-1][2] - prog[-2][2]) if len(prog) >= 2 else 0.0
                limit = max(STALE_LOG_S, 4.0 * dt_line)
                if silence > limit:
                    explained = suspended_overlap(pw.get(cid, []), stamp, time.time())
                    if explained >= 0.8 * silence:
                        notes.append("silence de %.1f h explique par la pause"
                                     % (silence / 3600))
                    else:
                        issues.append("log muet depuis %.1f h (cadence %.0f min)"
                                      % (silence / 3600, dt_line / 60))
        elif log and log.lower().endswith(".log") and os.path.exists(log):
            notes.append("log sans ligne it= (demarrage)")
        elif log and not log.lower().endswith(".log"):
            pass
        else:
            issues.append("pas de log de progression")

        ck = os.path.join(CKPT, cid + ".npz")
        if os.path.exists(ck):
            age_h = (time.time() - os.path.getmtime(ck)) / 3600
            it_ck = ""
            try:
                import numpy as np
                with np.load(ck) as d:
                    it_ck = " it=%d" % int(d["it"])
            except Exception:
                pass
            notes.append("ckpt%.1f h%s" % (age_h, it_ck))
        else:
            notes.append("pas encore de ckpt")

        dup = os.path.exists(os.path.join(SCAL, cid + ".json"))
        if dup:
            issues.append("JSON deja present (cas relance a tort ?)")

        cpu = sum(p[1] for p in e["procs"])
        mem = sum(p[2] for p in e["procs"])
        if prog:
            rows.append((time.time(), cid, prog[-1][0], prog[-1][1], cpu))
        st = "PASS" if not issues else ("FAIL" if any(
            "aucun processus" in i or "parallele" in i or "remonte" in i for i in issues) else "WARN")
        verdicts.append((st, cid, issues))
        if quiet and st == "PASS":
            continue
        print("%-40s %-6s %s" % (cid, st, " | ".join(notes + issues)))
        print("%-40s %-6s cpu=%.0f s / %.0f s  ram=%.0f Mo" %
              ("", "", cpu, WINDOW * max(1, len(busy)), mem))

    # cas interdits / inattendus
    banned = [c for c in cases if re.search(r"Re100_N512", c)]
    print("\ncas Re=100 a N=512 (interdits) :", banned or "aucun")

    # pilotes de campagne
    live = " ".join(p["cmd"] for p in b.values())
    for nm, pat in (("P3.2 (alpha)", r"run_campaign.*alpha-list 0\.05"),
                    ("watchdog P2.1", r"watch_p21")):
        print("%-18s %s" % (nm, "vivant" if re.search(pat, live) else "absent"))

    # ressources
    try:
        du = shutil.disk_usage(_ROOT)
        print("disque : %.1f Go libres sur %.1f Go" % (du.free / 2 ** 30, du.total / 2 ** 30))
    except OSError:
        pass
    for p in (SCAL, LOGS, CKPT, os.path.join(_ROOT, "results", "fields")):
        if not (os.path.isdir(p) and os.access(p, os.W_OK)):
            print("ATTENTION : dossier de sortie absent ou non inscriptible : %s" % p)

    n_fail = sum(1 for s, _, _ in verdicts if s == "FAIL")
    n_warn = sum(1 for s, _, _ in verdicts if s == "WARN")
    new = not os.path.exists(HISTORY)
    with open(HISTORY, "a", encoding="utf-8") as fh:
        if new:
            fh.write("epoch,case,it,dU,cpu_s\n")
        for t, cid, it, du, cpu in rows:
            fh.write("%.0f,%s,%d,%.6e,%.0f\n" % (t, cid, it, du, cpu))
    print("\nBILAN : %d PASS, %d WARN, %d FAIL sur %d cas"
          % (len(verdicts) - n_fail - n_warn, n_warn, n_fail, len(verdicts)))
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
