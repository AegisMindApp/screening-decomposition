#!/usr/bin/env python3
"""Re-dock the Mpro panel against the REPAIRED (protonated) receptor, keeping all 5 modes.

Protocol is read from the recorded bundle manifests rather than retyped, and the receptor md5
is asserted before any ligand is docked. Rules fixed in PREREGISTRATION_PROTONATED.md
(6bac07fd9) before this ran. Resumable: an existing non-empty output pose file is skipped.
"""
import json, glob, hashlib, os, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor

ROOT = "analysis/retrospective_benchmark"
BUN  = f"{ROOT}/kaggle/bundle_CHEMBL4523582_PROTONATED"
VINA = os.path.abspath(f"{ROOT}/bin/vina")
OUT  = "analysis/pose_ensemble/protonated_poses"
NPROC = int(os.environ.get("REDOCK_NPROC", "4"))
os.makedirs(OUT, exist_ok=True)

mans = sorted(glob.glob(f"{BUN}/shard*/manifest.json"))
assert mans, "no shard manifests"
proto = json.load(open(mans[0]))["protocol"]
for m in mans[1:]:
    assert json.load(open(m))["protocol"] == proto, "shards disagree on protocol"

def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""): h.update(b)
    return h.hexdigest()

# GATE: the receptor must be the repaired one, not the defective one.
for m in mans:
    r = os.path.join(os.path.dirname(m), "receptor.pdbqt")
    got = md5(r)
    assert got == proto["receptor_md5"], (
        f"{r}: md5 {got} != recorded {proto['receptor_md5']} - refusing to dock")
print(f"receptor md5 OK ({proto['receptor_md5'][:8]}...), "
      f"exh={proto['exhaustiveness']} modes={proto['num_modes']} "
      f"box={proto['box_center']} {proto['box_size']}")

jobs, labels = [], {}
for m in mans:
    d = os.path.dirname(m)
    for c in json.load(open(m))["compounds"]:
        lig = f"{d}/ligands/bench_{c['name']}.pdbqt"
        if os.path.exists(lig):
            jobs.append((c["name"], lig))
            labels[c["name"]] = c["label"]
print(f"{len(jobs)} ligands, {sum(labels.values())} active, {NPROC} workers")

ctr, size = proto["box_center"], proto["box_size"]
def dock(job):
    name, lig = job
    out = f"{OUT}/{name}.pdbqt"
    if os.path.exists(out) and os.path.getsize(out) > 0:
        return name, "cached"
    cmd = [VINA, "--receptor", os.path.join(os.path.dirname(os.path.dirname(lig)), "receptor.pdbqt"),
           "--ligand", lig,
           "--center_x", str(ctr[0]), "--center_y", str(ctr[1]), "--center_z", str(ctr[2]),
           "--size_x", str(size[0]), "--size_y", str(size[1]), "--size_z", str(size[2]),
           "--exhaustiveness", str(proto["exhaustiveness"]),
           "--num_modes", str(proto["num_modes"]),
           "--cpu", "1", "--out", out]
    try:
        subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                       timeout=600, check=True)
        return name, "ok"
    except Exception as e:
        if os.path.exists(out): os.remove(out)
        return name, f"FAIL {type(e).__name__}"

t0 = time.time(); done = fails = 0
with ThreadPoolExecutor(NPROC) as ex:
    for name, status in ex.map(dock, jobs):
        done += 1
        if status.startswith("FAIL"): fails += 1
        if done % 25 == 0 or done == len(jobs):
            el = time.time() - t0
            print(f"  {done}/{len(jobs)}  {el/60:.1f} min elapsed  "
                  f"ETA {(el/done)*(len(jobs)-done)/60:.1f} min  fails={fails}", flush=True)

json.dump(labels, open("analysis/pose_ensemble/protonated_labels.json", "w"), indent=0)
print(f"DONE {done-fails}/{len(jobs)} scored, {fails} failed, {(time.time()-t0)/60:.1f} min")
