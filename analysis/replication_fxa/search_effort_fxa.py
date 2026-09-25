#!/usr/bin/env python3
"""Arm 2: eight-fold search effort on Factor Xa (exh 32 vs exh 4), paired within seed.

Mirrors the Mpro estimator (analysis/seed_floor/paired_exhaustiveness.py): the difference is
taken WITHIN each seed on identical compounds, so what the two arms share cancels.
Floor: Factor Xa's own exh-32 95% bound, measured from these same six seeds -- not Mpro's, and
not Factor Xa's exh-4 floor. Rules: analysis/replication_fxa/PREREGISTRATION.md.
"""
import json, glob, math
import numpy as np
from scipy import stats
from sklearn.metrics import roc_auc_score

K = "analysis/retrospective_benchmark/kaggle/bundle_CHEMBL244_PROT"
lab = {c["name"]: c["label"] for m in glob.glob(f"{K}/shard*/manifest.json")
       for c in json.load(open(m))["compounds"]}
load = lambda pat: {json.load(open(f))["seed"]: json.load(open(f))["scores"] for f in glob.glob(pat)}
e32 = load("analysis/replication_fxa/arm2/seed_exh32_s*_shard0of1.json")
e4  = load("analysis/seed_floor/results_fxa/seed_exh4_s*_shard0of1.json")
SEEDS = sorted(set(e32) & set(e4))
names = sorted(set.intersection(*[set(e32[s]) for s in SEEDS],
                                *[set(e4[s]) for s in SEEDS]) & set(lab))
y = np.array([lab[n] for n in names])
print(f"seeds {SEEDS}  panel {len(names)}  ({y.mean():.1%} active)")

a32 = np.array([roc_auc_score(y, [-e32[s][n] for n in names]) for s in SEEDS])
a4  = np.array([roc_auc_score(y, [-e4[s][n]  for n in names]) for s in SEEDS])
d = a32 - a4
for s, x, z, dd in zip(SEEDS, a4, a32, d):
    print(f"  seed {s}: exh4 {x:.4f}  exh32 {z:.4f}  delta {dd:+.4f}")

# Factor Xa's OWN exh-32 floor, from these six seeds, chi2 one-sided 95% upper bound.
sd32 = float(a32.std(ddof=1)); df = len(a32) - 1
floor = sd32 * math.sqrt(df / stats.chi2.ppf(0.05, df))
n_ = len(SEEDS); se = d.std(ddof=1) / math.sqrt(n_)
t95 = stats.t.ppf(0.975, n_ - 1)
lo, hi = d.mean() - t95 * se, d.mean() + t95 * se
verdict = ("REPLICATES" if lo > floor else
           "DOES NOT REPLICATE" if hi < -floor else "INCONCLUSIVE")
print(f"\nexh-32 AUROC sd over {n_} seeds {sd32:.5f} -> 95% floor {floor:.5f} (Factor Xa's own)")
print(f"delta (exh32 - exh4) {d.mean():+.5f}  95% CI [{lo:+.5f}, {hi:+.5f}]")
print(f"Mpro reference: +0.0082 [+0.0016, +0.0148], INCONCLUSIVE there too")
print(f"PRE-REGISTERED READING: {verdict}")
json.dump(dict(seeds=SEEDS, n=len(names), auroc_exh4=[round(float(x),4) for x in a4],
               auroc_exh32=[round(float(x),4) for x in a32],
               per_seed_delta=[round(float(x),4) for x in d], delta=round(float(d.mean()),5),
               ci95=[round(float(lo),5), round(float(hi),5)], exh32_sd=round(sd32,5),
               floor=round(floor,5), verdict=verdict, mpro_reference=[0.0082,0.0016,0.0148]),
          open("analysis/replication_fxa/SEARCH_EFFORT_FXA.json","w"), indent=1)
