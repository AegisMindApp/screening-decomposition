#!/usr/bin/env python3
"""Pose-ensemble replication on Factor Xa. Mirrors analysis/pose_ensemble/analyse_protonated.py:
same 8 features, same RT, same 5-fold OOF logistic, same 10k paired bootstrap, same permuted
control. Only the target, the SMILES source and the floor change.

Floor: Factor Xa's OWN exh-4 95% bound (0.00703 on 886), rescaled to the panel that survives --
NOT Mpro's, and not the 0.0201 scoring floor the Mpro script hardcodes.
Rules: analysis/replication_fxa/PREREGISTRATION.md (amendment, 24 Sep).
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

FXA_EXH4_FLOOR = 0.00682   # 95% bound, n=12, 886 compounds -- analysis/seed_floor/fxa_floor_n12.json
# Was 0.00703 (n=10). Read from the file rather than hardcoded would be better still; pinned here
# with its provenance so the paper and the code cannot drift, which they had.
K = "analysis/retrospective_benchmark/kaggle"
SEED, NBOOT, RT = 20260924, 10000, 0.001987 * 298.15
RES = re.compile(r"REMARK VINA RESULT:\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)")

def modes(path):
    out = []
    for line in open(path):
        m = RES.search(line)
        if m: out.append((float(m.group(1)), float(m.group(2))))
    return out

def feats(ms):
    a = np.array([x[0] for x in ms], float); r = np.array([x[1] for x in ms], float)
    top = a.min()
    w = np.exp(-(a - top) / RT); w /= w.sum()
    gap = float(np.sort(a)[1] - top) if len(a) > 1 else 0.0
    near = a <= top + 1.0
    well = float(r[near].std()) if near.sum() > 1 else 0.0
    return [float(top), float((w * a).sum()), float(a.mean()), float(a.std()), gap,
            float(near.sum()), well, float(len(a))]

def descriptors(smi):
    m = Chem.MolFromSmiles(smi) if smi else None
    if m is None: return None
    return [Descriptors.MolWt(m), Descriptors.MolLogP(m), rdMolDescriptors.CalcNumHBD(m),
            rdMolDescriptors.CalcNumHBA(m), rdMolDescriptors.CalcNumRotatableBonds(m),
            rdMolDescriptors.CalcTPSA(m), Chem.GetFormalCharge(m)]

LAB = {c["name"]: c["label"]
       for m in glob.glob(f"{K}/bundle_CHEMBL244_PROT/shard*/manifest.json")
       for c in json.load(open(m))["compounds"]}
SMI = json.load(open(f"{K}/fxa_smiles.json"))

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
        if len(set(y[i].tolist())) > 1:
            out.append(roc_auc_score(y[i], pa[i]) - roc_auc_score(y[i], pb[i]))
    o = np.sort(out)
    return float(np.mean(o)), float(o[int(.025*len(o))]), float(o[int(.975*len(o))])

rows = []
for p in glob.glob("analysis/replication_fxa/poses_exh4_s1/*_out.pdbqt"):
    name = os.path.basename(p).replace("_out.pdbqt", "")
    if name not in LAB or name not in SMI: continue
    ms = modes(p)
    d = descriptors(SMI[name])
    if ms and d: rows.append((name, feats(ms), d, LAB[name]))
print(f"Factor Xa: {len(rows)} compounds with poses, labels and descriptors")

y = [r[3] for r in rows]; F = [r[1] for r in rows]; D = [r[2] for r in rows]
top = [[r[1][0]] for r in rows]
p_top, p_ens, p_desc = oof(top, y), oof(F, y), oof(D, y)
p_dens = oof([d + f for d, f in zip(D, F)], y)
rng = np.random.default_rng(SEED)
p_perm = oof(np.array(F)[rng.permutation(len(F))], y)

# DEGENERATE IDENTITY CONTROL (as on Mpro): ensemble restricted to its top-score feature
# must reproduce the top-pose arm exactly. If it does not, the two arms are not comparable.
p_id = oof([[f[0]] for f in F], y)
ident = float(np.max(np.abs(p_id - p_top)))
print(f"degenerate-identity control: max |Δ| = {ident:.2e}  ({'PASS' if ident < 1e-9 else 'FAIL'})")

a = lambda p: roc_auc_score(y, p)
d1, lo1, hi1 = boot(p_ens, p_top, y)
dp, lop, hip = boot(p_perm, p_top, y)
d3, lo3, hi3 = boot(p_dens, p_desc, y)
floor = FXA_EXH4_FLOOR * math.sqrt(886 / len(rows))
verdict = ("REPLICATES" if lo1 > floor else
           "DOES NOT REPLICATE" if hi1 < -floor else "INCONCLUSIVE")

print(f"\n  raw Vina top-pose AUROC   {roc_auc_score(y, [-t[0] for t in top]):.4f}")
print(f"  top pose (logistic)       {a(p_top):.4f}")
print(f"  ensemble (8 features)     {a(p_ens):.4f}")
print(f"  permuted control          {a(p_perm):.4f}")
print(f"  descriptors               {a(p_desc):.4f}")
print(f"  descriptors + ensemble    {a(p_dens):.4f}")
print(f"\n  1. ensemble - top pose      {d1:+.4f} [{lo1:+.4f}, {hi1:+.4f}]")
print(f"     permuted control         {dp:+.4f} [{lop:+.4f}, {hip:+.4f}]")
print(f"  3. desc+ens - descriptors   {d3:+.4f} [{lo3:+.4f}, {hi3:+.4f}]")
print(f"\n  floor (FXa exh-4, n={len(rows)}) = {floor:.5f}")
print(f"  Mpro reference: +0.0552 [+0.0434, +0.0680]")
print(f"  PRE-REGISTERED READING: {verdict}")

json.dump(dict(n=len(rows), identity_control=ident, auroc_top=a(p_top), auroc_ens=a(p_ens),
               auroc_perm=a(p_perm), auroc_desc=a(p_desc), auroc_desc_ens=a(p_dens),
               rule1=[d1, lo1, hi1], rule_perm=[dp, lop, hip], rule3=[d3, lo3, hi3],
               floor=floor, verdict=verdict, mpro_reference=[0.0552, 0.0434, 0.0680]),
          open("analysis/replication_fxa/POSE_ENSEMBLE_FXA.json", "w"), indent=1, default=float)
