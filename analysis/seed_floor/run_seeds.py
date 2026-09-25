#!/usr/bin/env python3
"""Dock the Mpro panel repeatedly under an UNCHANGED protocol, varying only --seed.

Implements PREREGISTRATION.md. The spread of AUROC across seeds is the run-to-run resolution
limit that `ci-2026-03151v` lacked.

Reuses the existing bundle rather than rebuilding one: same receptor, same vina binary, same
box, same compounds, same labels. Both md5s are asserted against the manifest before anything
runs — a resolution limit measured with a different binary would be measuring the wrong thing,
and this is the one experiment where "everything else identical" IS the experiment.

    python analysis/seed_floor/run_seeds.py --exh 4 --seeds 1 2 3 --jobs 4
    python analysis/seed_floor/run_seeds.py --exh 4 --seeds 1 2 3 --shard 0 --n-shards 3
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import subprocess
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

HERE = Path(__file__).resolve().parent
BUNDLE = HERE.parents[0] / "retrospective_benchmark" / "kaggle"
SRC = ["bundle_CHEMBL4523582_7VU6/shard0", "bundle_CHEMBL4523582_7VU6/shard1",
       "bundle_CHEMBL4523582_7VU6_overflow/shard0"]
# NOT `(-?\d+\.\d+)`. Vina prints a whole-number affinity without a decimal point -- observed
# live: CHEMBL4975056 scored "   1           -7" and the decimal-requiring pattern discarded a
# perfectly good result as "rc=0 no-score". The same pattern is in the original
# kaggle_shard_worker.py:32, so the manuscript's runs dropped those compounds too and counted
# them as failures. Rare (1 in ~250 here) and unrelated to activity, so unlikely to bias the
# AUROC -- but it is silent data loss, which is the failure class this whole exercise is about.
SCORE = re.compile(r"^\s*1\s+(-?\d+(?:\.\d+)?)\b", re.M)


def md5(p) -> str:
    return hashlib.md5(Path(p).read_bytes()).hexdigest()


def shard_dirs(root: Path):
    """Which shard directories to read.

    The Mpro paths stay hardcoded and take precedence, because the published stage-1 floor is
    computed from exactly those three and a discovery rule that happened to pick up a fourth
    would silently change a number the manuscript quotes. Only when none of them is present --
    i.e. this is some other target's bundle -- fall back to discovering `*/shardN` or `shardN`
    under the root.
    """
    named = [s for s in SRC if (root / s / "manifest.json").exists()]
    if named:
        return named
    found = sorted(q.parent.relative_to(root).as_posix()
                   for q in root.glob("*/shard*/manifest.json"))
    found += sorted(q.parent.relative_to(root).as_posix()
                    for q in root.glob("shard*/manifest.json"))
    return found


def collect(root: Path):
    """Every (name, label, ligand path) in the panel, plus the protocol they share."""
    recs, proto = {}, None
    for s in shard_dirs(root):
        d = root / s
        if not (d / "manifest.json").exists():
            continue
        m = json.loads((d / "manifest.json").read_text())
        proto = proto or dict(m["protocol"])
        for c in m["compounds"]:
            if c["name"] in recs:
                continue
            for stem in (f"bench_{c['name']}", c["name"]):
                for ext in ("pdbqt", "pdbqt.gz"):
                    p = d / "ligands" / f"{stem}.{ext}"
                    if p.exists():
                        recs[c["name"]] = {"label": int(c["label"]), "lig": str(p),
                                           "receptor": str(d / "receptor.pdbqt"),
                                           "vina": str(d / "vina")}
                        break
                if c["name"] in recs:
                    break
    return recs, proto


def dock(a):
    name, rec, exh, seed, ctr, size, modes, tmo = a
    cmd = [rec["vina"], "--receptor", rec["receptor"], "--ligand", rec["lig"],
           "--center_x", str(ctr[0]), "--center_y", str(ctr[1]), "--center_z", str(ctr[2]),
           "--size_x", str(size[0]), "--size_y", str(size[1]), "--size_z", str(size[2]),
           "--exhaustiveness", str(exh), "--num_modes", str(modes),
           "--cpu", "1", "--seed", str(seed),
           "--out", (os.path.join(POSES, f"{name}_out.pdbqt") if POSES else os.devnull)]
    t0 = time.time()
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=tmo)
        m = SCORE.search(r.stdout)
        if m:
            return name, float(m.group(1)), None
        why = (r.stderr or r.stdout or "").strip().splitlines()
        return name, None, f"rc={r.returncode} no-score: {(why[-1] if why else '')[:160]}"
    except subprocess.TimeoutExpired:
        return name, None, f"TIMEOUT {time.time()-t0:.0f}s exh={exh} seed={seed}"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--exh", type=int, required=True)
    ap.add_argument("--seeds", type=int, nargs="+", required=True)
    ap.add_argument("--jobs", type=int, default=os.cpu_count() or 4)
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--poses", default=None,
                    help="directory to KEEP the docked poses in, one {name}_out.pdbqt per "
                         "compound holding all --num_modes modes. Default discards them to "
                         "/dev/null, which is what the seed-floor runs want and is why no "
                         "Factor Xa pose set existed. Needed for the pose-ensemble and gnina "
                         "replications (analysis/replication_fxa/PREREGISTRATION.md).")
    ap.add_argument("--n-shards", type=int, default=1)
    ap.add_argument("--timeout", type=int, default=3600)
    ap.add_argument("--panel", type=int, default=0,
                    help="restrict to a FIXED random subset of N compounds, chosen "
                         "deterministically from the sorted panel with a fixed RNG so every "
                         "seed and both exhaustiveness settings see exactly the same "
                         "compounds. Stage 2 needs this: at exh=32 the full panel costs ~251 "
                         "core-hours per seed.")
    ap.add_argument("--bundle", default=str(BUNDLE))
    ap.add_argument("--out-dir", default=str(HERE / "results"))
    a = ap.parse_args()
    global POSES
    POSES = a.poses
    if POSES:
        os.makedirs(POSES, exist_ok=True)

    recs, proto = collect(Path(a.bundle))
    if not recs:
        raise SystemExit(f"no compounds found under {a.bundle}")
    # The protocol is the experiment. A different binary or receptor measures a different floor.
    vr, rr = next(iter(recs.values()))["vina"], next(iter(recs.values()))["receptor"]
    assert md5(vr) == proto["vina_md5"], f"vina md5 {md5(vr)} != manifest {proto['vina_md5']}"
    assert md5(rr) == proto["receptor_md5"], "receptor md5 differs from the manifest"
    print(f"panel {len(recs)} compounds | vina+receptor md5 match the manifest | "
          f"exh={a.exh} (manifest {proto['exhaustiveness']}) | seeds {a.seeds}", flush=True)

    allnames = sorted(recs)
    panel_md5 = ""
    if a.panel and a.panel < len(allnames):
        # Deterministic, and independent of shard and seed, so every shard of every seed sees
        # the identical subset. The ANALYSIS does not reimplement this selection -- it
        # intersects the compounds actually scored, so the two cannot drift apart.
        import random as _r
        allnames = sorted(_r.Random(20260919).sample(allnames, a.panel))
        panel_md5 = hashlib.md5(",".join(allnames).encode()).hexdigest()[:12]
        print(f"[panel] fixed {a.panel}-compound subset (rng 20260919), md5 {panel_md5}",
              flush=True)
    names = allnames[a.shard::a.n_shards]
    ctr, size, modes = proto["box_center"], proto["box_size"], proto["num_modes"]
    out_dir = Path(a.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    for seed in a.seeds:
        out = out_dir / f"seed_exh{a.exh}_s{seed}_shard{a.shard}of{a.n_shards}.json"
        done = {}
        if out.exists():                      # resume: a killed session costs nothing
            done = json.loads(out.read_text()).get("scores", {})
        todo = [n for n in names if n not in done]
        print(f"\nseed {seed}: {len(todo)} to dock, {len(done)} cached", flush=True)
        t0, fails = time.time(), {}
        args = [(n, recs[n], a.exh, seed, ctr, size, modes, a.timeout) for n in todo]
        # as_completed, NOT ex.map. map() yields results IN SUBMISSION ORDER, so a single
        # pathological ligand blocks every later result from being recorded -- measured here:
        # one ligand ran 13+ minutes at exh=4 while three workers finished dozens more whose
        # scores were neither written to the resume file nor counted. A crash in that window
        # loses all of them, and the printed rate is meaningless. Recording on completion fixes
        # the durability hole and makes the ETA honest.
        with ProcessPoolExecutor(a.jobs) as ex:
            futs = [ex.submit(dock, x) for x in args]
            for i, fut in enumerate(as_completed(futs), 1):
                name, score, err = fut.result()
                if score is not None:
                    done[name] = score
                else:
                    fails[name] = err
                if i % 25 == 0 or i == len(args):
                    el = time.time() - t0
                    print(f"  {i}/{len(args)}  {el/i:.1f}s/lig  "
                          f"eta {(len(args)-i)*el/i/60:.0f}min  fails {len(fails)}", flush=True)
                    out.write_text(json.dumps(
                        {"exh": a.exh, "seed": seed, "shard": a.shard, "n_shards": a.n_shards,
                         "protocol": proto, "labels": {n: recs[n]["label"] for n in done},
                         "scores": done, "failures": fails,
                         "complete": False}, indent=1))
        out.write_text(json.dumps(
            {"exh": a.exh, "seed": seed, "shard": a.shard, "n_shards": a.n_shards,
             "panel": a.panel or len(allnames), "panel_md5": panel_md5,
             "protocol": proto, "labels": {n: recs[n]["label"] for n in done},
             "scores": done, "failures": fails,
             # Over NAMES, not len(done): this shard resumes from a cache that may hold
             # compounds from an earlier, larger panel, and counting those would mark the
             # shard complete while compounds in its actual name list were still unscored.
             "complete": all(n in done or n in fails for n in names),
             "wall_s": round(time.time() - t0, 1)}, indent=1))
        print(f"  -> {out.name}  {len(done)} scored, {len(fails)} failed, "
              f"{(time.time()-t0)/60:.0f}min", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
