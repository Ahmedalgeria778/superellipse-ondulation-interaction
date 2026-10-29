# -*- coding: utf-8 -*-
"""
solver/pause_all.py
===================
Suspend / reprend les calculs LBM en cours sans les tuer.

  python -m solver.pause_all            # suspend tous les cas solver.run_case
  python -m solver.pause_all --resume   # reprend les processus suspendus
  python -m solver.pause_all --status   # etat courant

NtSuspendProcess gele le processus : CPU libere, etat (champs, iteration en
cours) conserve en memoire, reprise exacte au meme point. Les checkpoints sur
disque restent la securite en cas de redemarrage machine.
"""

from __future__ import annotations

import ctypes
import io
import json
import os
import re
import subprocess
import sys
import time

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
STATE = os.path.join(_ROOT, "results", "paused_cases.json")
JOURNAL = os.path.join(_ROOT, "results", "pause_journal.jsonl")

PROCESS_SUSPEND_RESUME = 0x0800
PROCESS_QUERY_LIMITED = 0x1000
MIN_CHILD_MB = 60.0   # le vrai calcul ; en dessous c'est le lanceur du venv

_ntdll = ctypes.WinDLL("ntdll")
_k32 = ctypes.WinDLL("kernel32", use_last_error=True)


def procs():
    q = ("Get-CimInstance Win32_Process -Filter \"Name='python.exe'\" | "
         "Select-Object ProcessId,WorkingSetSize,CommandLine | "
         "ConvertTo-Json -Compress")
    out = subprocess.run(["powershell", "-NoProfile", "-Command", q],
                         capture_output=True, text=True).stdout
    d = json.loads(out.strip() or "[]")
    return [d] if isinstance(d, dict) else d


def canon(n, alpha, Nw, Re, N, tol):
    """case_id identique a run_case.case_id (cle commune avec check_health)."""
    cid = ("se_n%g_w%.2f_Nw%d_Re%d_N%d"
           % (float(n), float(alpha), int(Nw), int(float(Re)), int(float(N))))
    return cid + ("" if float(tol) == 1e-9 else "_tol%.0e" % float(tol))


def cases():
    """case_id -> liste de (pid, memoire_Mo) pour les processus solver.run_case."""
    found = {}
    for p in procs():
        cmd = " ".join((p.get("CommandLine") or "").split())
        if "solver.run_case" not in cmd:
            continue
        m = re.search(r"--n\s+(\S+).*?--alpha\s+(\S+).*?--Re\s+(\S+).*?--N\s+(\S+)", cmd)
        if not m:
            continue
        n, alpha, Re, N = m.groups()
        mw = re.search(r"--Nw\s+(\S+)", cmd)
        mt = re.search(r"--tol\s+(\S+)", cmd)
        cid = canon(n, alpha, mw.group(1) if mw else 1, Re, N,
                    mt.group(1) if mt else 1e-9)
        found.setdefault(cid, []).append(
            (int(p["ProcessId"]), int(p.get("WorkingSetSize") or 0) / 2 ** 20))
    return found


def journal(event, case_ids):
    """Historique durable des pauses (check_health s'en sert pour les silences)."""
    rec = dict(event=event, at=time.time(), stamp=time.strftime("%Y-%m-%d %H:%M:%S"),
               cases=sorted(case_ids))
    with io.open(JOURNAL, "a", encoding="utf-8") as fh:
        fh.write(json.dumps(rec) + "\n")


def _open(pid):
    h = _k32.OpenProcess(PROCESS_SUSPEND_RESUME | PROCESS_QUERY_LIMITED, False, pid)
    return h if h else None


def suspend(pid):
    h = _open(pid)
    if not h:
        return False
    try:
        return _ntdll.NtSuspendProcess(h) == 0
    finally:
        _k32.CloseHandle(h)


def resume(pid):
    h = _open(pid)
    if not h:
        return False
    try:
        return _ntdll.NtResumeProcess(h) == 0
    finally:
        _k32.CloseHandle(h)


def cpu_of():
    q = ("Get-CimInstance Win32_Process -Filter \"Name='python.exe'\" | "
         "Select-Object ProcessId,UserModeTime,KernelModeTime | "
         "ConvertTo-Json -Compress")
    out = subprocess.run(["powershell", "-NoProfile", "-Command", q],
                         capture_output=True, text=True).stdout
    d = json.loads(out.strip() or "[]")
    d = [d] if isinstance(d, dict) else d
    return {int(p["ProcessId"]): (int(p.get("UserModeTime") or 0)
                                  + int(p.get("KernelModeTime") or 0)) / 1e7
            for p in d}


def main():
    found = cases()
    if not found:
        print("aucun calcul solver.run_case en cours")
        return 0

    if "--resume" in sys.argv:
        if not os.path.exists(STATE):
            print("aucun etat de pause enregistre (%s absent)" % STATE)
            return 1
        st = json.load(io.open(STATE, encoding="utf-8"))
        live = cpu_of()
        n, done = 0, []
        for rec in st["paused"]:
            if rec["pid"] in live and resume(rec["pid"]):
                print("[reprend] %-34s pid=%d" % (rec["case"], rec["pid"]))
                done.append(rec["case"])
                n += 1
        if done:
            journal("resume", done)
        os.remove(STATE)
        print("%d processus repris" % n)
        return 0

    if "--status" in sys.argv:
        st = json.load(io.open(STATE, encoding="utf-8")) if os.path.exists(STATE) else {"paused": []}
        known = {r["pid"] for r in st["paused"]}
        for cid in sorted(found):
            for pid, mb in found[cid]:
                tag = "SUSPENDU" if pid in known else "actif"
                print("  %-24s pid=%-7d %6.0f Mo  %s" % (cid, pid, mb, tag))
        return 0

    # suspension : uniquement le processus de calcul (le plus gros de la paire)
    paused, total_mb = [], 0.0
    for cid in sorted(found):
        kids = [x for x in found[cid] if x[1] >= MIN_CHILD_MB]
        if not kids:
            print("[attention] %s : pas de processus de calcul identifiable" % cid)
            continue
        pid, mb = max(kids, key=lambda x: x[1])
        if suspend(pid):
            paused.append(dict(case=cid, pid=pid, mem_mb=mb,
                               at=time.strftime("%Y-%m-%d %H:%M:%S")))
            total_mb += mb
            print("[suspend] %-24s pid=%-7d %6.0f Mo" % (cid, pid, mb))
        else:
            print("[ECHEC  ] %-24s pid=%d" % (cid, pid))

    before = cpu_of()
    time.sleep(12)
    after = cpu_of()
    busy = [(pid, after[pid] - before.get(pid, after[pid]))
            for pid in {r["pid"] for r in paused} if pid in after
            and after[pid] - before.get(pid, after[pid]) > 0.5]
    print("\nCPU residuel sur 12 s : %s"
          % (", ".join("pid=%d %.1fs" % b for b in busy) if busy else "aucun (0 s)"))
    print("memoire encore immobilisee : %.0f Mo" % total_mb)
    if paused:
        journal("pause", [r["case"] for r in paused])
    json.dump({"paused": paused, "at": time.strftime("%Y-%m-%d %H:%M:%S")},
              io.open(STATE, "w", encoding="utf-8"), indent=1)
    print("etat enregistre dans results/paused_cases.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
