#!/usr/bin/env python3
"""Arm 1: receptor-preparation repair on Factor Xa (unprotonated vs obabel -p 7.4 protonated).

Paired within seed on the compounds common to both receptors and all three seeds, mirroring the
Mpro arm (analysis/receptor_prep/RESULT.md, +0.0451 [+0.0236, +0.0671]).
Floor: Factor Xa's OWN exh-4 95% bound, rescaled to the panel that survives.
Rules: analysis/replication_fxa/PREREGISTRATION.md.
"""
import json, glob, math
import numpy as np
from sklearn.metrics import roc_auc_score

FXA_EXH4_FLOOR = 0.00682   # 95% bound, n=12, 886 compounds -- analysis/seed_floor/fxa_floor_n12.json
# Was 0.00703 (n=10). Read from the file rather than hardcoded would be better still; pinned here
# with its provenance so the paper and the code cannot drift, which they had.
K = "analysis/retrospective_benchmark/kaggle/bundle_CHEMBL244_PROT"
lab = {c["name"]: c["label"] for m in glob.glob(f"{K}/shard*/manifest.json")
       for c in json.load(open(m))["compounds"]}
load = lambda pat: {json.load(open(f))["seed"]: json.load(open(f))["scores"]
                    for f in glob.glob(pat)}
raw  = load("analysis/seed_floor/results_fxa_raw/seed_exh4_s*.json")
prot = load("analysis/seed_floor/results_fxa/seed_exh4_s*.json")
SEEDS = sorted(set(raw) & set(prot))
names = sorted(set.intersection(*[set(raw[s]) for s in SEEDS],
                                *[set(prot[s]) for s in SEEDS]) & set(lab))
y = np.array([lab[n] for n in names])

# DROPOUT CONTROL -- the unprotonated seeds were cut off by Kaggle's 12h limit, so the panel is
# whatever finished, not a random sample. Class-selected dropout would void the comparison.
miss = sorted(set(lab) - set(names))
dd = (np.mean([lab[n] for n in miss]) - y.mean()) if miss else 0.0
print(f"seeds {SEEDS}  panel {len(names)}  dropped {len(miss)}  active-fraction diff {dd:+.1%}")
assert abs(dd) <= 0.10, "VOID: dropout is class-selected beyond 10pp"

# Paired within seed: AUROC(protonated) - AUROC(unprotonated), same compounds, same seed.
per = [roc_auc_score(y, [-prot[s][n] for n in names]) - roc_auc_score(y, [-raw[s][n] for n in names])
       for s in SEEDS]
for s, d in zip(SEEDS, per):
    print(f"  seed {s}: prot {roc_auc_score(y,[-prot[s][n] for n in names]):.4f}  "
          f"unprot {roc_auc_score(y,[-raw[s][n] for n in names]):.4f}  delta {d:+.4f}")
delta = float(np.mean(per))

# Compound bootstrap, holding the seed pairing.
rng = np.random.default_rng(20260925); bs = []
for _ in range(10000):
    i = rng.integers(0, len(y), len(y))
    if 0 < y[i].sum() < len(i):
        bs.append(np.mean([roc_auc_score(y[i], [-prot[s][names[j]] for j in i])
                           - roc_auc_score(y[i], [-raw[s][names[j]] for j in i]) for s in SEEDS]))
lo, hi = np.percentile(bs, [2.5, 97.5])
floor = FXA_EXH4_FLOOR * math.sqrt(886 / len(names))
verdict = ("REPLICATES" if lo > floor else
           "DOES NOT REPLICATE" if hi < -floor else "INCONCLUSIVE")
print(f"\ndelta (protonated - unprotonated) {delta:+.4f}  95% CI [{lo:+.4f}, {hi:+.4f}]")
print(f"floor (FXa exh-4, n={len(names)}) {floor:.5f}")
print(f"Mpro reference: +0.0451 [+0.0236, +0.0671]")
print(f"PRE-REGISTERED READING: {verdict}")
json.dump(dict(seeds=SEEDS, n=len(names), dropped=len(miss), dropout_diff=round(float(dd),4),
               per_seed=[round(d,4) for d in per], delta=round(delta,4),
               ci95=[round(float(lo),4), round(float(hi),4)], floor=round(floor,5),
               verdict=verdict, mpro_reference=[0.0451,0.0236,0.0671]),
          open("analysis/replication_fxa/RECEPTOR_PREP_FXA.json","w"), indent=1)
