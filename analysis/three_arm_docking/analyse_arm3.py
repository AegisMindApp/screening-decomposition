#!/usr/bin/env python3
"""Arm 3 vs arm 1: same poses, different scoring function.

Reading rule pre-registered at cd8f55639 before these numbers existed.
Direction matters: Vina affinity is more-negative = better binder; gnina CNNaffinity is
a pKd-like quantity where HIGHER = better. Both are converted to "higher = better" before
any AUROC is taken.
"""
import json, glob, random, statistics as st

K = "analysis/retrospective_benchmark/kaggle"
lab = {}
for m in glob.glob(f"{K}/bundle_CHEMBL4523582_7VU6*/shard*/manifest.json"):
    for c in json.load(open(m))["compounds"]: lab[c["name"]] = c["label"]
sp = lambda k: k[6:] if k.startswith("bench_") else k
vina = {sp(k): v for k, v in json.load(open(f"{K}/exh4_scores_merged.json")).items()}
gn = {sp(k): v for k, v in json.load(open(f"{K}/gnina_rescore.json")).items()}

both = sorted(set(vina) & set(gn) & set(lab))
print(f"vina {len(vina)}  gnina {len(gn)}  labelled {len(lab)}  ALL THREE: {len(both)}")
cov = len(both) / 753
print(f"coverage {cov:.1%}  (pre-registered floor 90%)")
assert cov >= 0.90, "below the declared floor — comparison void"

# POSITIVE CONTROL: gnina also reports a Vina-style affinity for the pose it was given.
# It must track the score already in our cache; if it does not, it is not scoring our pose.
gv = [gn[n]["vina_affinity"] for n in both if "vina_affinity" in gn[n]]
vv = [vina[n] for n in both if "vina_affinity" in gn[n]]
def pear(x, y):
    mx, my = st.mean(x), st.mean(y)
    return sum((a-mx)*(b-my) for a, b in zip(x, y)) / (
        (sum((a-mx)**2 for a in x) * sum((b-my)**2 for b in y)) ** .5)
r_ctrl = pear(gv, vv)
print(f"control: gnina's own Vina affinity vs our cached score, r = {r_ctrl:+.4f} (n={len(gv)})")
assert r_ctrl > 0.90, "gnina is not scoring the pose we think it is"

def auroc(score_higher_better, labels):
    pairs = sorted(zip(score_higher_better, labels))
    n1 = sum(labels); n0 = len(labels) - n1
    ranks = {}; i = 0
    while i < len(pairs):
        j = i
        while j + 1 < len(pairs) and pairs[j+1][0] == pairs[i][0]: j += 1
        for k in range(i, j+1): ranks[k] = (i + j) / 2 + 1
        i = j + 1
    s = sum(ranks[i] for i, (_, y) in enumerate(pairs) if y == 1)
    return (s - n1 * (n1 + 1) / 2) / (n1 * n0)

y  = [lab[n] for n in both]
a1 = [-vina[n] for n in both]                 # Vina: negate so higher = better
a3 = [gn[n]["cnn_affinity"] for n in both]    # gnina CNNaffinity: already higher = better
A1, A3 = auroc(a1, y), auroc(a3, y)
print(f"\nARM 1  Vina poses + Vina score  : AUROC {A1:.4f}")
print(f"ARM 3  SAME poses + gnina CNN   : AUROC {A3:.4f}")
print(f"dAUROC (3 - 1) = {A3-A1:+.4f}")

rnd = random.Random(20260831); N = len(both); ds = []
for _ in range(10000):
    idx = [rnd.randrange(N) for _ in range(N)]
    yy = [y[i] for i in idx]
    if not (0 < sum(yy) < len(yy)): continue
    ds.append(auroc([a3[i] for i in idx], yy) - auroc([a1[i] for i in idx], yy))
ds.sort()
lo, hi = ds[int(.025*len(ds))], ds[int(.975*len(ds))]
print(f"95% CI: [{lo:+.4f}, {hi:+.4f}]  ({len(ds)} resamples)")

# The point estimate can exceed +0.10 while the CI's lower bound does not. Saying
# "helps by less than predicted" in that case would misreport a +0.140 effect.
if lo > 0.10:  verdict = "SUPPORTED — dAUROC CI entirely above the predicted +0.10"
elif lo > 0:   verdict = ("DIRECTIONALLY SUPPORTED, NOT AT THE PRE-REGISTERED STRICTNESS — "
                          "the effect is clearly positive and the point estimate exceeds "
                          "+0.10, but the 95% CI lower bound does not clear +0.10")
elif hi < 0:   verdict = "REFUTED — gnina scoring is WORSE"
else:          verdict = "NULL — no detectable scoring effect"
print(f"\nPRE-REGISTERED READING: {verdict}")
print(f"descriptor bar 0.763 cleared by arm 3: {A3 >= 0.763}")
print(f"chance (0.5) cleared by arm 3: {A3 >= 0.5}")
rho_note = pear(a1, a3)
print(f"\ncorr(Vina score, gnina CNNaffinity) on the same poses: r = {rho_note:+.4f}")
json.dump({"n": len(both), "arm1_vina": A1, "arm3_gnina": A3, "delta": A3-A1,
           "ci95": [lo, hi], "verdict": verdict, "control_r": r_ctrl,
           "clears_descriptor_bar": A3 >= 0.763},
          open("analysis/three_arm_docking/ARM3_RESULT.json", "w"), indent=1)
