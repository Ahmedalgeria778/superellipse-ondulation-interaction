# -*- coding: utf-8 -*-
"""Campagne FEM P1.2 v2 — avec saut des cas deja calcules.

Reprend la liste des 12 cas P1.2 mais :
  * saute tout cas dont le JSON existe deja dans results/fem/ ;
  * attend la fin du cas nx=256 deja en cours (processus orphelin) sans le
    relancer ; si le processus orphelin est mort sans produire de JSON,
    relance le cas.
"""
import json, os, subprocess, sys, time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PY = sys.executable
DRV = os.path.join(ROOT, "solver", "run_fem_xcheck.py")
FEM_DIR = os.path.join(ROOT, "results", "fem")

JOBS = [
    (100, 2, 64), (100, 2, 128), (100, 2, 256),
    (100, 4, 64), (100, 4, 128), (100, 4, 256),
    (1000, 2, 64), (1000, 2, 128),
    (1000, 4, 64), (1000, 4, 128),
]

# nx=256 abandonne a Re=1000 : la reference n2/Re=100 en nx=256 a coute 11,5 h
# contre 1493 s en nx=128, avec seulement ~1,3 GB de RAM libre pendant les
# campagnes LBM. La convergence maillage est etablie a Re=100 (64/128/256) et le
# cross-check du contraste (2,4) a Re=1000 est couvert par nx=64/128.

# Cas deja lance par un driver precedemment tue -> attendre son JSON.
INFLIGHT = (100, 4, 256)

PS_CHECK = (
    "Get-CimInstance Win32_Process -Filter \"Name='python.exe'\" | "
    "Where-Object { $_.CommandLine -match 'run_fem_xcheck.*nx 256' } | "
    "Measure-Object | Select-Object -ExpandProperty Count"
)


def json_of(Re, n, nx):
    return os.path.join(FEM_DIR, "FEMX-n%d-a0.00-R%d-nx%d.json" % (n, Re, nx))


def inflight_alive():
    try:
        out = subprocess.run(["powershell", "-NoProfile", "-Command", PS_CHECK],
                             capture_output=True, text=True, timeout=60).stdout
        return int(out.strip() or 0) > 0
    except Exception:
        return True  # en cas de doute, on attend


def wait_inflight():
    target = json_of(*INFLIGHT)
    print("[v2] attend la fin du cas inflight %s..." % os.path.basename(target),
          flush=True)
    while not os.path.exists(target):
        if not inflight_alive():
            print("[v2] processus inflight mort sans JSON -> relance", flush=True)
            return False
        time.sleep(180)
    print("[v2] cas inflight termine (JSON present)", flush=True)
    return True


def main():
    relance_inflight = not os.path.exists(json_of(*INFLIGHT))
    todo = []
    for job in JOBS:
        if os.path.exists(json_of(*job)):
            print("[v2] SAUT (deja calcule)  Re=%d n=%d nx=%d" % job, flush=True)
        elif job == INFLIGHT and relance_inflight:
            ok = wait_inflight()
            if ok:
                print("[v2] SAUT (inflight termine)  Re=%d n=%d nx=%d" % job,
                      flush=True)
            else:
                todo.append(job)
        else:
            todo.append(job)

    print("[v2] %d cas restant(s) a lancer" % len(todo), flush=True)
    for i, (Re, n, nx) in enumerate(todo, 1):
        cmd = [PY, DRV, "--Re", str(Re), "--n", str(n), "--alpha", "0.0",
               "--nx", str(nx), "--ny", str(nx), "--N", "128"]
        print("[v2 %d/%d] Re=%d n=%d nx=%d" % (i, len(todo), Re, n, nx),
              flush=True)
        r = subprocess.run(cmd, cwd=ROOT)
        print("  rc=%d %s" % (r.returncode,
                              "OK" if r.returncode == 0 else "ERREUR"),
              flush=True)

    print("[v2] CAMPAGNE FEM P1.2 TERMINEE", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
