#!/usr/bin/env python3
"""Does restoring the receptor's hydrogen-bond donors change the Mpro result?

Reading rule pre-registered at 1f68b20c4 BEFORE these numbers existed.
Paired within compound: same 753, same labels, same box, same vina binary,
exhaustiveness 4. Protonation is the only difference.
"""
import json, glob, random, statistics as st, hashlib

K = "analysis/retrospective_benchmark/kaggle"
lab = {}
for m in glob.glob(f"{K}/bundle_CHEMBL4523582_7VU6*/shard*/manifest.json"):
    for c in json.load(open(m))["compounds"]: lab[c["name"]] = c["label"]
sp = lambda k: k[6:] if k.startswith("bench_") else k
old = {sp(k): v for k, v in json.load(open(f"{K}/exh4_scores_merged.json")).items()}
new = {sp(k): v for k, v in json.load(open(f"{K}/protonated_scores_merged.json")).items()}

# The receptors must actually differ, and only in protonation.
h_old = hashlib.md5(open(f"{K}/bundle_CHEMBL4523582_7VU6/shard0/receptor.pdbqt","rb").read()).hexdigest()
h_new = hashlib.md5(open(f"{K}/mpro_receptor_H.pdbqt","rb").read()).hexdigest()
assert h_old != h_new, "receptors identical — nothing was tested"
print(f"receptor md5  unprotonated {h_old[:8]}   protonated {h_new[:8]}   differ: True")

both = sorted(set(old) & set(new) & set(lab))
cov = len(both) / 753
print(f"unprotonated {len(old)}  protonated {len(new)}  paired with labels: {len(both)}  ({cov:.1%})")
assert cov >= 0.95, "below the pre-registered 95% floor — comparison void"

def auroc(s, y):
    p = sorted(zip(s, y)); n1 = sum(y); n0 = len(y) - n1
    r = {}; i = 0
    while i < len(p):
        j = i
        while j + 1 < len(p) and p[j+1][0] == p[i][0]: j += 1
        for k in range(i, j+1): r[k] = (i + j) / 2 + 1
        i = j + 1
    return (sum(r[i] for i, (_, yy) in enumerate(p) if yy == 1) - n1*(n1+1)/2) / (n1*n0)

y = [lab[n] for n in both]
a_old = [-old[n] for n in both]      # Vina: more negative = better, so negate
a_new = [-new[n] for n in both]
A0, A1 = auroc(a_old, y), auroc(a_new, y)
print(f"\nAUROC unprotonated (0 donors) : {A0:.4f}")
print(f"AUROC protonated  (727 donors): {A1:.4f}")
print(f"dAUROC = {A1-A0:+.4f}")

rnd = random.Random(20260831); N = len(both); ds = []
for _ in range(10000):
    idx = [rnd.randrange(N) for _ in range(N)]
    yy = [y[i] for i in idx]
    if not (0 < sum(yy) < len(yy)): continue
    ds.append(auroc([a_new[i] for i in idx], yy) - auroc([a_old[i] for i in idx], yy))
ds.sort(); lo, hi = ds[int(.025*len(ds))], ds[int(.975*len(ds))]
print(f"95% CI: [{lo:+.4f}, {hi:+.4f}]  ({len(ds)} resamples)")

if lo > 0:   v = "PREPARATION WAS A REAL PART OF THE FAILURE — every docking number this project has published is affected"
elif hi < 0: v = "THE UNPROTONATED RECEPTOR WAS BETTER — evidence the H-bond term is not doing useful work here"
else:        v = "DEFECT REAL BUT NOT THE CAUSE — Vina's H-bond term is not what is wrong on this target"
print(f"\nPRE-REGISTERED READING: {v}")
print(f"  clears chance (0.5) : {A1 >= 0.5}")
print(f"  clears the 0.763 descriptor bar: {A1 >= 0.763}")

def spearman(a, b):
    def rk(v):
        o = sorted(range(len(v)), key=lambda i: v[i]); r = [0]*len(v); i = 0
        while i < len(o):
            j = i
            while j+1 < len(o) and v[o[j+1]] == v[o[i]]: j += 1
            for k in range(i, j+1): r[o[k]] = (i+j)/2 + 1
            i = j+1
        return r
    ra, rb = rk(a), rk(b); n = len(a); ma, mb = sum(ra)/n, sum(rb)/n
    return (sum((x-ma)*(z-mb) for x, z in zip(ra, rb)) /
            ((sum((x-ma)**2 for x in ra)*sum((z-mb)**2 for z in rb))**.5))
rho = spearman([old[n] for n in both], [new[n] for n in both])
top = lambda d, k=10: set(sorted(both, key=lambda n: d[n])[:k])
print(f"\nSpearman rho between the two score vectors: {rho:+.4f}")
print(f"mean |score change|: {st.mean(abs(old[n]-new[n]) for n in both):.3f} kcal/mol")
print(f"mean unprotonated {st.mean([old[n] for n in both]):+.3f}   protonated {st.mean([new[n] for n in both]):+.3f}")
for k in (10, 50): print(f"top-{k} overlap: {len(top(old,k)&top(new,k))}/{k}")
json.dump({"n": len(both), "auroc_unprotonated": A0, "auroc_protonated": A1,
           "delta": A1-A0, "ci95": [lo, hi], "verdict": v, "spearman": rho,
           "clears_chance": A1 >= 0.5, "clears_descriptor_bar": A1 >= 0.763},
          open("analysis/receptor_prep/RESULT.json","w"), indent=1)
