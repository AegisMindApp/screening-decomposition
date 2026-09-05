#!/usr/bin/env python3
"""Bootstrap the resolution floor -- and separate the two things it was conflating.

MARGINAL_VALUE.md 3 reports the floor as 0.039 AUROC and describes it as gnina's
--score_only re-evaluation of "the identical Vina poses". Reproducing it shows the cached
side is the exhaustiveness-32 run, while gnina rescored the exhaustiveness-4 poses. Those
are two independent stochastic searches, not identical poses.

Both quantities are real and both are useful, but they measure different things:

  PROTOCOL floor  (exh32 cache vs gnina on exh4 poses) -- run-to-run reproducibility of the
                  whole docking protocol. The right floor for comparing two published
                  protocols, which is what most reported improvements are.
  SCORING floor   (exh4 cache vs gnina on the SAME exh4 poses) -- pure numerical/scoring
                  noise with placement held fixed. The right floor for comparing two
                  scoring functions on shared poses.

Positive control: the protocol arm must reproduce MARGINAL_VALUE.md 3 before either number
is read.
"""
import json, glob, random, statistics as st

K = "analysis/retrospective_benchmark/kaggle"
lab = {}
for m in glob.glob(f"{K}/bundle_CHEMBL4523582_7VU6*/shard*/manifest.json"):
    for c in json.load(open(m))["compounds"]: lab[c["name"]] = c["label"]
sp = lambda k: k[6:] if k.startswith("bench_") else k
gn = {sp(k): v for k, v in json.load(open(f"{K}/gnina_rescore.json")).items()}
exh4 = {sp(k): v for k, v in json.load(open(f"{K}/exh4_scores_merged.json")).items()}
exh32 = {}
for f in sorted(glob.glob(f"{K}/vina_cache_CHEMBL4523582_7VU6_shard*.json")):
    for k, v in json.load(open(f)).items(): exh32[sp(k)] = v

def auroc(s, y):
    p = sorted(zip(s, y)); n1 = sum(y); n0 = len(y) - n1
    if n1 == 0 or n0 == 0: return float("nan")
    r = {}; i = 0
    while i < len(p):
        j = i
        while j + 1 < len(p) and p[j+1][0] == p[i][0]: j += 1
        for k in range(i, j+1): r[k] = (i + j) / 2 + 1
        i = j + 1
    return (sum(r[i] for i, (_, t) in enumerate(p) if t == 1) - n1*(n1+1)/2) / (n1*n0)

def pear(x, y):
    mx, my = st.mean(x), st.mean(y)
    return sum((a-mx)*(b-my) for a, b in zip(x, y)) / (
        (sum((a-mx)**2 for a in x) * sum((b-my)**2 for b in y)) ** .5)

def arm(name, cache, seed):
    both = sorted(n for n in set(cache) & set(gn) & set(lab) if "vina_affinity" in gn[n])
    y = [lab[n] for n in both]
    c = [cache[n] for n in both]
    g = [gn[n]["vina_affinity"] for n in both]
    A_c, A_g = auroc([-v for v in c], y), auroc([-v for v in g], y)
    r = pear(c, g); md = st.mean([b - a for a, b in zip(c, g)])
    gt2 = sum(abs(a - b) > 2 for a, b in zip(c, g))
    rnd = random.Random(seed); N = len(both); ds = []
    for _ in range(10000):
        idx = [rnd.randrange(N) for _ in range(N)]
        yy = [y[i] for i in idx]
        if not (0 < sum(yy) < len(yy)): continue
        ds.append(auroc([-c[i] for i in idx], yy) - auroc([-g[i] for i in idx], yy))
    ds.sort(); lo, hi = ds[int(.025*len(ds))], ds[int(.975*len(ds))]
    print(f"\n=== {name}   n = {len(both)}")
    print(f"  r = {r:.4f}   mean diff = {md:+.3f} kcal/mol   |diff| > 2 kcal: {gt2}/{len(both)}")
    print(f"  AUROC cached Vina = {A_c:.4f}   AUROC gnina Vina term = {A_g:.4f}")
    print(f"  GAP = {A_c-A_g:+.4f}   95% CI [{lo:+.4f}, {hi:+.4f}]  ({len(ds)} resamples)")
    return dict(n=len(both), r=r, mean_diff_kcal=md, n_gt_2kcal=gt2,
                auroc_cached=A_c, auroc_gnina_vina_term=A_g, gap=A_c-A_g, ci95=[lo, hi])

print("PROTOCOL floor -- exh32 cache vs gnina --score_only on exh4 poses.")
print("This is what MARGINAL_VALUE.md 3 actually measured. Control targets:")
print("  r 0.9369 | mean -0.043 | 0/745 >2kcal | AUROC 0.4182 / 0.3793 | gap 0.039")
p = arm("PROTOCOL floor (different pose sets, same scoring function)", exh32, 20260906)
ok = (abs(p["r"]-0.9369) < 2e-3 and abs(p["auroc_cached"]-0.4182) < 2e-3
      and abs(p["auroc_gnina_vina_term"]-0.3793) < 2e-3 and p["n"] == 745)
print(f"  CONTROL: {'PASS' if ok else 'FAIL'}")
assert ok, "protocol arm does not reproduce the published numbers - nothing below is readable"

print("\nSCORING floor -- exh4 cache vs gnina --score_only on THE SAME exh4 poses.")
s = arm("SCORING floor (identical poses, two implementations)", exh4, 20260907)

json.dump({"protocol_floor": p, "scoring_floor": s,
           "receptor": "unprotonated (donor-defective)",
           "interpretation":
             "The published 0.039 floor compares the exhaustiveness-32 Vina cache against "
             "gnina's --score_only re-evaluation of the exhaustiveness-4 poses. Those are two "
             "independent stochastic searches, so 0.039 is the run-to-run reproducibility of "
             "the whole docking PROTOCOL, not of a scoring function on fixed poses. With "
             "placement genuinely held fixed the gap is about half that. Both are real floors "
             "and each gates a different comparison; the paper's qualitative readings are "
             "unchanged under either, but receptor repair (+0.045) moves from marginal to "
             "clearly resolvable under the scoring floor."},
          open("analysis/docking_value/FLOOR_INTERVAL.json", "w"), indent=1)
print("\nwrote analysis/docking_value/FLOOR_INTERVAL.json")
