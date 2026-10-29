# -*- coding: utf-8 -*-
"""Diagnostic lecture seule de tous les calculs en cours."""
import io
import json
import os
import re
import subprocess
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCAL = os.path.join(ROOT, "results", "scalars")
LOGS = os.path.join(ROOT, "results", "logs")


def ps(query):
    out = subprocess.run(["powershell", "-NoProfile", "-Command", query],
                         capture_output=True, text=True).stdout
    return out


def procs():
    q = ("Get-CimInstance Win32_Process -Filter \"Name='python.exe'\" | "
         "Select-Object ProcessId,UserModeTime,KernelModeTime,CommandLine | "
         "ConvertTo-Json -Compress")
    d = json.loads(ps(q).strip() or "[]")
    if isinstance(d, dict):
        d = [d]
    return d


def flag(cmd, name):
    m = re.search(r"--%s\s+(\S+)" % name, cmd)
    return m.group(1) if m else None


def case_name(n, alpha, nw, Re, N):
    return "se_n%s_w%.2f_Nw%s_Re%s_N%s_tol1e-08" % (n, float(alpha), nw, Re, N)


def last_iter(tag):
    p = os.path.join(LOGS, tag + ".log")
    if not os.path.exists(p):
        return None, p
    t = io.open(p, encoding="utf-8", errors="replace").read()
    m = re.findall(r"it=\s*(\d+)\s+dU=([\d.eE+-]+)\s+t=([\d.]+)s", t)
    if not m:
        return ("lance, aucune ligne it=", os.path.getsize(p)), p
    it, dU, cpu = m[-1]
    return ("it=%9s  dU=%-9s  cpu_ours=%.2f  (%d lignes)"
            % (it, dU, float(cpu) / 3600.0, len(m))), p


def json_status(tag):
    for suf in ("_tol1e-08", ""):
        p = os.path.join(SCAL, tag + suf + ".json")
        if os.path.exists(p):
            try:
                d = json.load(io.open(p, encoding="utf-8"))
            except ValueError:
                return "JSON illisible"
            return "%s  res=%s" % (d.get("status"), d.get("residual"))
    return "pas de JSON"


def main():
    a = {int(p["ProcessId"]): int(p.get("UserModeTime") or 0)
         + int(p.get("KernelModeTime") or 0) for p in procs()}
    print("... mesure CPU sur 25 s ...", flush=True)
    time.sleep(25)
    rows = procs()
    now = {int(p["ProcessId"]): int(p.get("UserModeTime") or 0)
           + int(p.get("KernelModeTime") or 0) for p in rows}

    cases = []
    for p in rows:
        cmd = (p.get("CommandLine") or "").replace("\r", " ")
        if "solver.run_case" not in cmd:
            continue
        pid = int(p["ProcessId"])
        delta = (now.get(pid, 0) - a.get(pid, now.get(pid, 0))) / 1e7
        cases.append((pid, delta, cmd))

    # regrouper par nom de cas (stub + enfant partagent la ligne de commande)
    by = {}
    for pid, delta, cmd in cases:
        n = flag(cmd, "n")
        alpha = flag(cmd, "alpha")
        nw = flag(cmd, "Nw")
        Re = flag(cmd, "Re")
        N = flag(cmd, "N")
        tag = case_name(n, alpha, nw, Re, N)
        e = by.setdefault(tag, {"cpu": 0.0, "pids": [], "started": None})
        e["cpu"] += delta
        e["pids"].append(pid)
        m = re.search(r"--log\s+(\S+)", cmd)
        e["started"] = e["started"] or (m.group(1) if m else None)

    print("\n=== CAS LBM EN VOL (%d) ===" % len(by))
    for tag in sorted(by):
        e = by[tag]
        tot = sum(1 for x in by.values())
        info, path = last_iter(tag)
        print("%-42s pids=%s  cpu_25s=%5.1fs" % (tag, e["pids"], e["cpu"]))
        print("      log %s -> %s" % (os.path.basename(path), info))
        print("      sortie JSON : %s" % json_status(tag))
    print("\n=== AUTRES PROCESSUS python ===")
    for p in rows:
        cmd = (p.get("CommandLine") or "").replace("\r", " ")
        if "solver.run_case" in cmd:
            continue
        if not re.search(r"run_campaign|orchestrat|watch|fem|FEM", cmd):
            continue
        pid = int(p["ProcessId"])
        d = (now.get(pid, 0) - a.get(pid, now.get(pid, 0))) / 1e7
        print("  %-7s cpu_25s=%5.1fs  %s" % (pid, d, " ".join(cmd.split())[:150]))

    print("\n=== PILOTES DE CAMPAGNE ATTENDUS ===")
    live = " ".join((p.get("CommandLine") or "") for p in rows)
    for nm, pat in [("P3.1 (Nw=2)", r"Re-list 400,1000"),
                    ("P3.2 (alpha=0.05/0.20)", r"--n-list 2,4"),
                    ("P2.1 N=512 (interdit Re=100)", r"N 512.*run_campaign|"
                     r"run_campaign.*N 512"),
                    ("watchdog Re400 n2 N512", r"watch_p21")]:
        print("  %-28s %s" % (nm, "VIVANT" if re.search(pat, live) else "absent"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
