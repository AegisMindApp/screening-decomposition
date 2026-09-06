#!/usr/bin/env python3
"""Proper reconstruction of the pose-ensemble arm.

Rules fixed in PREREGISTRATION_RECONSTRUCTION.md (824270441) before any estimate was read.
Tests H1: the effect is unstable under cross-validation fold assignment, and the published
+0.0194 and the reimplementation's +0.0393 are both draws from its sampling distribution.

Every seed drawn is reported. No seed is selected for reproducing a target value.
"""
import json, glob, os, re, sys
import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, rdMolDescriptors
from sklearn.model_selection import StratifiedKFold
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.metrics import roc_auc_score
RDLogger.DisableLog("rdApp.*")

K = "analysis/retrospective_benchmark/kaggle"
SEEDS = [20260907 + i for i in range(40)]
FLOOR = 0.0201
RT = 0.001987 * 298.15
RES = re.compile(r"REMARK VINA RESULT:\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)")
strip = lambda n: n[6:] if n.startswith("bench_") else n

def modes(p):
    out = []
    for line in open(p):
        m = RES.search(line)
        if m: out.append((float(m.group(1)), float(m.group(2))))
    return out

def feats(ms):
    a = np.array([x[0] for x in ms], float); r = np.array([x[1] for x in ms], float)
    top = a.min(); w = np.exp(-(a - top) / RT); w /= w.sum()
    gap = float(np.sort(a)[1] - top) if len(a) > 1 else 0.0
    near = a <= top + 1.0
    well = float(r[near].std()) if near.sum() > 1 else 0.0
    return [float(top), float((w * a).sum()), float(a.mean()), float(a.std()),
            gap, float(near.sum()), well, float(len(a))]

def desc(s):
    m = Chem.MolFromSmiles(s) if s else None
    if m is None: return None
    return [Descriptors.MolWt(m), Descriptors.MolLogP(m), rdMolDescriptors.CalcNumHBD(m),
            rdMolDescriptors.CalcNumHBA(m), rdMolDescriptors.CalcNumRotatableBonds(m),
            rdMolDescriptors.CalcTPSA(m), Chem.GetFormalCharge(m)]

LAB = {}
for m in glob.glob(f"{K}/bundle_CHEMBL4523582_PROTONATED/shard*/manifest.json"):
    for c in json.load(open(m))["compounds"]: LAB[c["name"]] = c["label"]
SMI = json.load(open(f"{K}/mpro_smiles.json"))

def build(d):
    rows = []
    for p in glob.glob(f"{d}/*.pdbqt"):
        n = strip(os.path.basename(p).replace("_out.pdbqt", "").replace(".pdbqt", ""))
        if n not in LAB or n not in SMI: continue
        ms = modes(p)
        if not ms: continue
        dd = desc(SMI[n])
        if dd: rows.append((feats(ms), dd, LAB[n]))
    return rows

def oof(X, y, seed):
    X = np.asarray(X, float); y = np.asarray(y, int); p = np.zeros(len(y))
    for tr, te in StratifiedKFold(5, shuffle=True, random_state=seed).split(X, y):
        mdl = make_pipeline(StandardScaler(), LogisticRegression(max_iter=5000))
        mdl.fit(X[tr], y[tr]); p[te] = mdl.predict_proba(X[te])[:, 1]
    return p

def run(rows, tag):
    F = [r[0] for r in rows]; D = [r[1] for r in rows]; y = [r[2] for r in rows]
    TOP = [[f[0]] for f in F]
    deltas, descs, degen, perml, permf = [], [], [], [], []
    rng0 = np.random.default_rng(12345)
    for sd in SEEDS:
        a_top = roc_auc_score(y, oof(TOP, y, sd))
        a_ens = roc_auc_score(y, oof(F, y, sd))
        deltas.append(a_ens - a_top)
        descs.append(roc_auc_score(y, oof(D, y, sd)))
        # control 1: ensemble restricted to its top-score feature IS the top-pose arm
        degen.append(abs(roc_auc_score(y, oof([[f[0]] for f in F], y, sd)) - a_top))
        # control 3: permuted labels
        yp = list(rng0.permutation(y))
        perml.append(roc_auc_score(yp, oof(F, yp, sd)))
        # control 4: permuted feature rows
        Fp = np.array(F)[rng0.permutation(len(F))]
        permf.append(roc_auc_score(y, oof(Fp, y, sd)) - a_top)
    d = np.array(deltas)
    lo, hi = float(np.percentile(d, 2.5)), float(np.percentile(d, 97.5))
    out = dict(tag=tag, n=len(rows), mean=float(d.mean()), sd=float(d.std()),
               lo=lo, hi=hi, min=float(d.min()), max=float(d.max()),
               desc_min=float(min(descs)), desc_max=float(max(descs)),
               degen_max=float(max(degen)), perm_label_mean=float(np.mean(perml)),
               perm_feat_mean=float(np.mean(permf)), deltas=[float(x) for x in d])
    print(f"\n=== {tag}  (n={out['n']}, {len(SEEDS)} fold seeds) ===")
    print(f"  ensemble - top pose : mean {out['mean']:+.4f}  sd {out['sd']:.4f}")
    print(f"                        95% seed range [{lo:+.4f}, {hi:+.4f}]   "
          f"full [{out['min']:+.4f}, {out['max']:+.4f}]")
    print(f"  descriptor baseline : [{out['desc_min']:.4f}, {out['desc_max']:.4f}]  "
          f"(published 0.7651)")
    print(f"  CONTROL degenerate identity : max |diff| {out['degen_max']:.2e}  "
          f"{'PASS' if out['degen_max'] < 1e-9 else 'FAIL'}")
    print(f"  CONTROL permuted labels     : {out['perm_label_mean']:.4f}  "
          f"{'PASS' if 0.45 <= out['perm_label_mean'] <= 0.55 else 'FAIL'}")
    print(f"  CONTROL permuted features   : {out['perm_feat_mean']:+.4f}  "
          f"{'PASS' if out['perm_feat_mean'] <= 0 else 'FAIL'}")
    print(f"  CONTROL desc range spans 0.7651 : "
          f"{'PASS' if out['desc_min'] <= 0.7651 <= out['desc_max'] else 'FAIL'}")
    return out

rep = run(build("analysis/pose_ensemble/kaggle_poses"), "REPAIRED receptor")
dfc = run(build(f"{K}/exh4_poses"), "DEFECTIVE receptor")

def verdict(o):
    if o["lo"] > FLOOR:  return "RESOLVABLE, above floor"
    if o["hi"] < FLOOR:  return "RESOLVABLE, below floor"
    return "UNRESOLVED (H1 supported)"

print("\n=== VERDICTS against the pre-registered rules ===")
for o in (rep, dfc):
    print(f"  {o['tag']:20s} {verdict(o)}")
print(f"\n  Do the published +0.0194 and the reimplementation's +0.0393 both fall inside the")
print(f"  DEFECTIVE seed range [{dfc['min']:+.4f}, {dfc['max']:+.4f}]?  "
      f"+0.0194: {dfc['min'] <= 0.0194 <= dfc['max']}   "
      f"+0.0393: {dfc['min'] <= 0.0393 <= dfc['max']}")
json.dump({"seeds": SEEDS, "floor": FLOOR, "repaired": rep, "defective": dfc,
           "verdict_repaired": verdict(rep), "verdict_defective": verdict(dfc)},
          open("analysis/pose_ensemble/RECONSTRUCTION.json", "w"), indent=1)
print("\nwrote analysis/pose_ensemble/RECONSTRUCTION.json")
