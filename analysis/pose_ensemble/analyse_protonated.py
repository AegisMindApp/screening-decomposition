#!/usr/bin/env python3
"""Pose-ensemble scoring on the REPAIRED receptor.

Rules fixed in PREREGISTRATION_PROTONATED.md (6bac07fd9) before the re-dock started.

Both pose sets are run through IDENTICAL code: the defective-receptor poses (exh4_poses, still
on disk) and the repaired-receptor poses just returned from Kaggle. Rule 2 is therefore a real
paired comparison rather than a comparison against a published number produced by other code.
"""
import json, glob, os, re, math
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
SEED, NBOOT, RT = 20260906, 10000, 0.001987 * 298.15   # kcal/mol at 298.15 K
FLOOR_SCORING = 0.0201
RES = re.compile(r"REMARK VINA RESULT:\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)")
strip = lambda n: n[6:] if n.startswith("bench_") else n

def modes(path):
    """(affinity, rmsd_lb) per MODEL, in file order."""
    out = []
    with open(path) as f:
        for line in f:
            m = RES.search(line)
            if m: out.append((float(m.group(1)), float(m.group(2))))
    return out

def feats(ms):
    """The 8 ensemble features described in RESULT.md."""
    a = np.array([x[0] for x in ms], float)
    r = np.array([x[1] for x in ms], float)
    top = a.min()
    w = np.exp(-(a - top) / RT); w /= w.sum()
    boltz = float((w * a).sum())
    gap = float(np.sort(a)[1] - top) if len(a) > 1 else 0.0
    near = a <= top + 1.0
    well = float(r[near].std()) if near.sum() > 1 else 0.0
    return [float(top), boltz, float(a.mean()), float(a.std()), gap,
            float(near.sum()), well, float(len(a))]

def descriptors(smi):
    m = Chem.MolFromSmiles(smi) if smi else None
    if m is None: return None
    return [Descriptors.MolWt(m), Descriptors.MolLogP(m),
            rdMolDescriptors.CalcNumHBD(m), rdMolDescriptors.CalcNumHBA(m),
            rdMolDescriptors.CalcNumRotatableBonds(m), rdMolDescriptors.CalcTPSA(m),
            Chem.GetFormalCharge(m)]

LAB = {}
for m in glob.glob(f"{K}/bundle_CHEMBL4523582_PROTONATED/shard*/manifest.json"):
    for c in json.load(open(m))["compounds"]: LAB[c["name"]] = c["label"]
SMI = json.load(open(f"{K}/mpro_smiles.json"))

def oof(X, y, seed=SEED):
    X = np.asarray(X, float); y = np.asarray(y, int); p = np.zeros(len(y))
    for tr, te in StratifiedKFold(5, shuffle=True, random_state=seed).split(X, y):
        mdl = make_pipeline(StandardScaler(), LogisticRegression(max_iter=5000))
        mdl.fit(X[tr], y[tr]); p[te] = mdl.predict_proba(X[te])[:, 1]
    return p

def boot(pa, pb, y, n=NBOOT, seed=SEED):
    rng = np.random.default_rng(seed); y = np.asarray(y); N = len(y); out = []
    for _ in range(n):
        i = rng.integers(0, N, N)
        if len(set(y[i].tolist())) < 2: continue
        out.append(roc_auc_score(y[i], pa[i]) - roc_auc_score(y[i], pb[i]))
    o = np.sort(out)
    return float(np.mean(o)), float(o[int(.025*len(o))]), float(o[int(.975*len(o))])

def build(posedir, tag):
    rows = []
    for p in glob.glob(f"{posedir}/*.pdbqt"):
        name = strip(os.path.basename(p).replace("_out.pdbqt", "").replace(".pdbqt", ""))
        if name not in LAB or name not in SMI: continue
        ms = modes(p)
        if not ms: continue
        d = descriptors(SMI[name])
        if d is None: continue
        rows.append((name, feats(ms), d, LAB[name]))
    print(f"{tag}: {len(rows)} compounds with poses, labels and descriptors")
    return rows

def analyse(rows, tag):
    y = [r[3] for r in rows]
    F = [r[1] for r in rows]; D = [r[2] for r in rows]
    top = [[r[1][0]] for r in rows]
    p_top, p_ens, p_desc = oof(top, y), oof(F, y), oof(D, y)
    p_dens = oof([d + f for d, f in zip(D, F)], y)
    rng = np.random.default_rng(SEED)
    Fp = np.array(F)[rng.permutation(len(F))]
    p_perm = oof(Fp, y)
    a = lambda p: roc_auc_score(y, p)
    r = dict(tag=tag, n=len(rows), auroc_top=a(p_top), auroc_ens=a(p_ens),
             auroc_desc=a(p_desc), auroc_desc_ens=a(p_dens), auroc_perm=a(p_perm),
             rule1=boot(p_ens, p_top, y), rule_perm=boot(p_perm, p_top, y),
             rule3=boot(p_dens, p_desc, y),
             raw_top_auroc=roc_auc_score(y, [-t[0] for t in top]))
    print(f"\n=== {tag} (n={r['n']}) ===")
    print(f"  raw Vina top-pose AUROC   {r['raw_top_auroc']:.4f}")
    print(f"  top pose (logistic)       {r['auroc_top']:.4f}")
    print(f"  ensemble (8 features)     {r['auroc_ens']:.4f}")
    print(f"  permuted control          {r['auroc_perm']:.4f}")
    print(f"  descriptors               {r['auroc_desc']:.4f}")
    print(f"  descriptors + ensemble    {r['auroc_desc_ens']:.4f}")
    for k, lbl in (("rule1", "1. ensemble - top pose"), ("rule_perm", "   permuted control"),
                   ("rule3", "3. desc+ens - descriptors")):
        d, lo, hi = r[k]
        print(f"  {lbl:28s} {d:+.4f} [{lo:+.4f}, {hi:+.4f}]")
    return r

prot = build("analysis/pose_ensemble/kaggle_poses", "REPAIRED receptor")
defe = build(f"{K}/exh4_poses", "DEFECTIVE receptor")
rp, rd = analyse(prot, "REPAIRED receptor"), analyse(defe, "DEFECTIVE receptor")

print("\n=== POSITIVE CONTROL (pre-registered) ===")
pub = 0.4530; got = rp["raw_top_auroc"]
ok = abs(got - pub) <= 0.010
print(f"  repaired top-mode raw AUROC {got:.4f} vs published {pub:.4f} "
      f"(tol 0.010) -> {'PASS' if ok else 'FAIL'}")
cov = rp["n"] / 751
print(f"  coverage {cov:.1%} (floor 90%) -> {'PASS' if cov >= 0.90 else 'FAIL'}")

d, lo, hi = rp["rule1"]
verdict = "SUPPORTED" if (d > FLOOR_SCORING and lo > 0) else "REFUTED"
print(f"\n=== PRIMARY VERDICT ===")
print(f"  rule 1: {d:+.4f} [{lo:+.4f}, {hi:+.4f}] vs scoring floor {FLOOR_SCORING} -> {verdict}")
print(f"  rule 2 (repair effect on this arm): repaired {rp['rule1'][0]:+.4f} "
       f"vs defective {rd['rule1'][0]:+.4f}")

json.dump({"control_pass": bool(ok and cov >= 0.90), "coverage": cov,
           "floor_scoring": FLOOR_SCORING, "verdict": verdict,
           "repaired": rp, "defective": rd},
          open("analysis/pose_ensemble/PROTONATED_RESULT.json", "w"), indent=1, default=float)
print("\nwrote analysis/pose_ensemble/PROTONATED_RESULT.json")
