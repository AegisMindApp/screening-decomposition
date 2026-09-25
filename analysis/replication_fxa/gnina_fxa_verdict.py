"""gnina replication on Factor Xa. Rules: analysis/replication_fxa/PREREGISTRATION.md (amendment)."""
import json, glob, math
import numpy as np
from sklearn.metrics import roc_auc_score

FXA_EXH4_FLOOR = 0.00682   # 95% bound, n=12, 886 compounds -- analysis/seed_floor/fxa_floor_n12.json
# Was 0.00703 (n=10). Read from the file rather than hardcoded would be better still; pinned here
# with its provenance so the paper and the code cannot drift, which they had.
lab = {c["name"]: c["label"]
       for m in glob.glob("analysis/retrospective_benchmark/kaggle/bundle_CHEMBL244_PROT/shard*/manifest.json")
       for c in json.load(open(m))["compounds"]}
vina = json.load(open("analysis/replication_fxa/gcp_seed_exh4_s1.json"))["scores"]
g = json.load(open("analysis/replication_fxa/gnina_fxa.json"))
names = sorted(set(vina) & set(g["scores"]) & set(lab))
cov = len(names) / len(lab)

y = np.array([lab[n] for n in names])
v = -np.array([vina[n] for n in names])                        # higher = better
c = np.array([g["scores"][n]["cnn_affinity"] for n in names])  # higher = better
ctrl = np.corrcoef([g["scores"][n]["vina_affinity"] for n in names],
                   [vina[n] for n in names])[0, 1]

rng = np.random.default_rng(20260924); d = []
for _ in range(10000):
    i = rng.integers(0, len(y), len(y))
    if 0 < y[i].sum() < len(i):
        d.append(roc_auc_score(y[i], c[i]) - roc_auc_score(y[i], v[i]))
lo, hi = np.percentile(d, [2.5, 97.5])
a_v, a_c = roc_auc_score(y, v), roc_auc_score(y, c)
floor = FXA_EXH4_FLOOR * math.sqrt(886 / len(names))

verdict = ("VOID" if cov < 0.90 or ctrl <= 0.90 else
           "REPLICATES" if lo > floor else
           "DOES NOT REPLICATE" if hi < -floor else "INCONCLUSIVE")
out = dict(n=len(names), coverage=round(cov, 4), failures=len(g["failures"]),
           control_r=round(float(ctrl), 4), auroc_vina=round(a_v, 4), auroc_gnina=round(a_c, 4),
           delta=round(a_c - a_v, 4), ci95=[round(float(lo), 4), round(float(hi), 4)],
           floor=round(floor, 5), verdict=verdict, mpro_reference=dict(delta=0.1397, ci95=[0.0913, 0.1867]))
json.dump(out, open("analysis/replication_fxa/GNINA_FXA_RESULT.json", "w"), indent=1)
for k, val in out.items(): print(f"{k:15} {val}")
