#!/usr/bin/env python3
"""Minimise and score the recovered DiffDock poses against the REPAIRED receptor.

Rules fixed in PREREGISTRATION_DIFFDOCK_PROT.md (ce636ea0a) before any score existed.
No GPU: DiffDock's poses were recovered from the original Kaggle kernel outputs, and
protonation does not change pose generation, only the Vina treatment of the poses.
"""
import json, glob, os, re, subprocess, hashlib, sys
from concurrent.futures import ThreadPoolExecutor

K = "analysis/retrospective_benchmark/kaggle"
BUN = f"{K}/bundle_CHEMBL4523582_PROTONATED/shard0"
# Use the binary the bundle ships and whose md5 the protocol records, not the repo copy.
VINA = os.path.abspath(f"{BUN}/vina")
POSES = "analysis/three_arm_docking/diffdock_poses"
OUT = "analysis/three_arm_docking/diffdock_prot_scored"
os.makedirs(OUT, exist_ok=True)
proto = json.load(open(f"{BUN}/manifest.json"))["protocol"]
REC = f"{BUN}/receptor.pdbqt"

h = hashlib.md5(open(REC, "rb").read()).hexdigest()
assert h == proto["receptor_md5"], f"receptor md5 {h} != {proto['receptor_md5']}"
vh = hashlib.md5(open(VINA, "rb").read()).hexdigest()
assert vh == proto["vina_md5"], f"vina md5 {vh} != recorded {proto['vina_md5']}"
print(f"receptor md5 OK ({h[:8]}...)  box={proto['box_center']} {proto['box_size']}")

ctr, size = proto["box_center"], proto["box_size"]
# Vina 1.2.x --local_only reports the score as "Estimated Free Energy of Binding".
SCORE = re.compile(r"Estimated Free Energy of Binding\s*:\s*(-?\d+\.\d+)"
                   r"|Affinity:\s*(-?\d+\.\d+)", re.M)

def one(sdf):
    cid = os.path.basename(sdf)[:-4]
    outp = f"{OUT}/{cid}.json"
    if os.path.exists(outp): return json.load(open(outp)).get("score")
    pq = f"/tmp/dd_{cid}.pdbqt"
    try:
        r = subprocess.run(["obabel", sdf, "-O", pq, "--partialcharge", "gasteiger"],
                           capture_output=True, timeout=120)
        if not os.path.exists(pq) or os.path.getsize(pq) == 0: return None
        r = subprocess.run([VINA, "--receptor", REC, "--ligand", pq,
                            "--center_x", str(ctr[0]), "--center_y", str(ctr[1]),
                            "--center_z", str(ctr[2]), "--size_x", str(size[0]),
                            "--size_y", str(size[1]), "--size_z", str(size[2]),
                            "--local_only", "--cpu", "1"],
                           capture_output=True, text=True, timeout=600)
        m = SCORE.search(r.stdout)
        if not m:
            blob = r.stdout + r.stderr
            why = ("outside_box" if "outside the grid box" in blob
                   or "not in the search space" in blob else "no_score")
            json.dump({"score": None, "why": why}, open(outp, "w"))
            return None
        s = float(m.group(1) or m.group(2))
        json.dump({"score": s}, open(outp, "w"))
        return s
    except Exception:
        return None
    finally:
        if os.path.exists(pq): os.remove(pq)

sdfs = sorted(glob.glob(f"{POSES}/*.sdf"))
print(f"{len(sdfs)} DiffDock poses to minimise and score")
ok = 0
with ThreadPoolExecutor(4) as ex:
    for i, s in enumerate(ex.map(one, sdfs), 1):
        if s is not None: ok += 1
        if i % 100 == 0: print(f"  {i}/{len(sdfs)} scored={ok}", flush=True)
print(f"DONE {ok}/{len(sdfs)} scored ({100*ok/len(sdfs):.1f}%)")
raw = {os.path.basename(p)[:-5]: json.load(open(p)) for p in glob.glob(f"{OUT}/*.json")}
res = {k: v["score"] for k, v in raw.items() if v.get("score") is not None}
import collections
why = collections.Counter(v.get("why") for v in raw.values() if v.get("score") is None)
print("failure reasons:", dict(why))
assert len(res) >= 0.75 * len(sdfs), (   # Amendment 1: published arm achieved 80.4%
    f"coverage {len(res)}/{len(sdfs)} below the amended 75% floor - run is void")
json.dump(res, open("analysis/three_arm_docking/diffdock_prot_scores.json", "w"), indent=0)
print(f"wrote diffdock_prot_scores.json ({len(res)})")
