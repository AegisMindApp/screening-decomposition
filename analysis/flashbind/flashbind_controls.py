#!/usr/bin/env python3
"""The pre-registered controls for the FlashBind arm, plus the checks a below-chance raw
AUROC demands before it is reported as a finding rather than a bug.

PREREGISTRATION.md requires, before the residual is readable:
  1. descriptor baseline on these compounds reproduces 0.7654 +/- 0.005
  2. coverage >= 90% of 751          (measured on Kaggle: 750/751 = 0.999)
  3. shuffled labels give AUROC in [0.45, 0.55]
"""
import json, numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
from sklearn.model_selection import StratifiedKFold
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score

K = "analysis/retrospective_benchmark/kaggle"
SEEDS = [20260906 + i for i in range(12)]

def descriptors(smi):
    m = Chem.MolFromSmiles(smi)
    if m is None: return None
    return [Descriptors.MolWt(m), Descriptors.MolLogP(m),
            rdMolDescriptors.CalcNumHBD(m), rdMolDescriptors.CalcNumHBA(m),
            rdMolDescriptors.CalcNumRotatableBonds(m), rdMolDescriptors.CalcTPSA(m),
            Chem.GetFormalCharge(m)]

smi = json.load(open(f"{K}/mpro_smiles.json"))
lab = json.load(open("analysis/docking_value/mpro_labels.json"))
raw = json.load(open("analysis/flashbind/results/flashbind_scores_751.json"))
fb  = {k.split("_", 1)[1]: v for k, v in raw["ensemble"].items()}
pm  = [{k.split("_", 1)[1]: v for k, v in m.items()} for m in raw["per_model"]]

rows = []
for n, s in smi.items():
    if n not in fb or n not in lab: continue
    d = descriptors(s)
    if d is None: continue
    rows.append((n, fb[n], lab[n], d, [m.get(n) for m in pm]))
if len(rows) < 100:
    raise SystemExit(f"only {len(rows)} rows joined - an empty join is not a result")

names = [r[0] for r in rows]
sc = np.array([r[1] for r in rows], float)
y  = np.array([r[2] for r in rows], int)
X  = np.array([r[3] for r in rows], float)
print(f"joined {len(rows)} compounds, {y.sum()} active ({y.mean():.3f})")

def oof(Z, seed):
    """The published procedure, verbatim from analysis/pose_ensemble/reconstruct.py:
    StandardScaler + LogisticRegression under StratifiedKFold(5). My first attempt used
    gradient boosting under plain KFold and returned 0.8982 against a published 0.7651 - the
    control failed because the control was wrong, not because the panel was. A positive
    control only works if it reproduces the procedure as well as the data."""
    Z = np.asarray(Z, float); p = np.zeros(len(y))
    for tr, te in StratifiedKFold(5, shuffle=True, random_state=seed).split(Z, y):
        m = make_pipeline(StandardScaler(), LogisticRegression(max_iter=5000))
        m.fit(Z[tr], y[tr]); p[te] = m.predict_proba(Z[te])[:, 1]
    return p

def oof_desc(seed):        return oof(X, seed)
def oof_desc_plus(e, seed): return oof(np.column_stack([X, e]), seed)

print("\n=== CONTROL 1  descriptor baseline (pre-registered 0.7654 +/- 0.005; published 0.7651)")
d_auc = [roc_auc_score(y, oof_desc(s)) for s in SEEDS]
print(f"  {np.mean(d_auc):.4f}  [{min(d_auc):.4f}, {max(d_auc):.4f}] over {len(SEEDS)} seeds")
c1 = abs(np.mean(d_auc) - 0.7654) <= 0.005
print(f"  CONTROL 1: {'PASS' if c1 else 'FAIL'}  (delta {abs(np.mean(d_auc)-0.7654):.4f})")

print("\n=== CONTROL 3  permutation (pre-registered [0.45, 0.55])")
rng = np.random.default_rng(20260908)
perm = [roc_auc_score(rng.permutation(y), sc) for _ in range(200)]
print(f"  shuffled-label AUROC mean {np.mean(perm):.4f}  "
      f"2.5-97.5% [{np.percentile(perm,2.5):.4f}, {np.percentile(perm,97.5):.4f}]")
c3 = 0.45 <= np.mean(perm) <= 0.55
print(f"  CONTROL 3: {'PASS' if c3 else 'FAIL'}")

print("\n=== RAW AUROC, and whether the inversion is real")
print(f"  ensemble        {roc_auc_score(y, sc):.4f}")
for i, m in enumerate(pm):
    v = np.array([m[n] for n in names], float)
    print(f"  checkpoint {i+1}    {roc_auc_score(y, v):.4f}")
# Is it significantly below chance? Bootstrap the AUROC itself.
bs = [roc_auc_score(y[i], sc[i]) for i in
      (rng.integers(0, len(y), len(y)) for _ in range(2000))
      if len(set(y[i])) > 1]
print(f"  bootstrap 2.5-97.5%   [{np.percentile(bs,2.5):.4f}, {np.percentile(bs,97.5):.4f}]")
print(f"  0.5 inside that interval: {np.percentile(bs,2.5) <= 0.5 <= np.percentile(bs,97.5)}")
print(f"  inverted score (1-p)  {roc_auc_score(y, -sc):.4f}")

print("\n=== What the score is tracking instead")
for j, nm in enumerate(["MolWt","LogP","HBD","HBA","RotB","TPSA","charge"]):
    print(f"  spearman(score, {nm:6s}) = "
          f"{np.corrcoef(np.argsort(np.argsort(sc)), np.argsort(np.argsort(X[:,j])))[0,1]:+.3f}")
print(f"  spearman(label, MolWt)   = "
      f"{np.corrcoef(y, np.argsort(np.argsort(X[:,0])))[0,1]:+.3f}")

print("\n=== MARGINAL VALUE over descriptors (Boltz-2 reference +0.093 [+0.062, +0.125])")
# Swept over all 12 fold seeds, not quoted at one. residual_band.py's own note applies here:
# residual and increment AUROCs move with the fold assignment by 0.01-0.03, which is the size
# of the effect being measured. The first attempt quoted SEEDS[0], where the descriptor
# baseline happens to sit at the MINIMUM of its own 12-seed range - i.e. the seed most
# flattering to any increment measured against it.
deltas, bases, pluses = [], [], []
for sd in SEEDS:
    b = oof_desc(sd); q = oof_desc_plus(sc, sd)
    ab, aq = roc_auc_score(y, b), roc_auc_score(y, q)
    bases.append(ab); pluses.append(aq); deltas.append(aq - ab)
deltas = np.array(deltas)
print(f"  descriptors      {np.mean(bases):.4f}  [{min(bases):.4f}, {max(bases):.4f}]")
print(f"  desc+FlashBind   {np.mean(pluses):.4f}  [{min(pluses):.4f}, {max(pluses):.4f}]")
print(f"  delta            {deltas.mean():+.4f}  seed range [{deltas.min():+.4f}, "
      f"{deltas.max():+.4f}]")
print(f"  seed range straddles zero: {deltas.min() < 0 < deltas.max()}")
print(f"  clear of Boltz-2's +0.062 lower bound: {deltas.max() < 0.062}")
a0, a1 = float(np.mean(bases)), float(np.mean(pluses))

print("\n=== IS THE SCORE JUST A SIZE PROXY?")
mw_auc = roc_auc_score(y, -X[:, 0])
rb_auc = roc_auc_score(y, -X[:, 4])
print(f"  AUROC of -MolWt alone          {mw_auc:.4f}")
print(f"  AUROC of -rotatable bonds alone {rb_auc:.4f}")
print(f"  FlashBind raw                  {roc_auc_score(y, sc):.4f}")
print(f"  size alone explains it: {abs(mw_auc - roc_auc_score(y, sc)) < 0.03}")

print("\n=== COMPARATOR CONTEXT: raw AUROC of the docking arms on this same panel")
import os as _os
_rb = "analysis/docking_value/RESIDUAL_BAND.json"
if _os.path.exists(_rb):
    for nm, o in json.load(open(_rb))["arms"].items():
        print(f"  {nm:38s} raw {o['raw']:.4f}")

json.dump({"n": len(rows), "n_active": int(y.sum()),
           "control_1_descriptor_baseline": {"mean": float(np.mean(d_auc)),
               "range": [float(min(d_auc)), float(max(d_auc))],
               "target": 0.7654, "pass": bool(c1)},
           "control_3_permutation": {"mean": float(np.mean(perm)), "pass": bool(c3)},
           "raw_auroc": {"ensemble": float(roc_auc_score(y, sc)),
               "per_checkpoint": [float(roc_auc_score(y, np.array([m[n] for n in names], float)))
                                  for m in pm],
               "bootstrap_ci": [float(np.percentile(bs,2.5)), float(np.percentile(bs,97.5))]},
           "marginal_value": {"base": a0, "plus": a1, "delta": float(deltas.mean()),
               "seed_range": [float(deltas.min()), float(deltas.max())],
               "straddles_zero": bool(deltas.min() < 0 < deltas.max())},
           "size_proxy": {"auroc_neg_molwt": float(mw_auc),
                          "auroc_neg_rotb": float(rb_auc)}},
          open("analysis/flashbind/results/CONTROLS.json", "w"), indent=1)
