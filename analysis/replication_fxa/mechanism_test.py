#!/usr/bin/env python3
"""Is gnina a constant-quality ranker, or does it track target difficulty?
Rules: analysis/replication_fxa/PREREGISTRATION_MECHANISM.md, fixed before this ran.

Regresses gnina's ABSOLUTE AUROC on Vina's across scaffold-binned sub-panels. Slope ~0 means
gnina ranks at fixed quality regardless of target difficulty (supports the broken-baseline
reading of the +0.140); slope ~1 means it tracks the target as a Vina re-run does (refutes it).
The delta-on-baseline regression is deliberately NOT used: its slope is -1 by construction.
"""
import json, glob, warnings
import numpy as np
from rdkit import Chem, RDLogger, DataStructs
from rdkit.Chem import AllChem
from rdkit.ML.Cluster import Butina
from sklearn.metrics import roc_auc_score
RDLogger.DisableLog("rdApp.*"); warnings.filterwarnings("ignore")
K = "analysis/retrospective_benchmark/kaggle"
MIN_BIN, SEED = 80, 20260925
strip = lambda n: n[6:] if n.startswith("bench_") else n

def bins_for(names, smi):
    fps = [AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(smi[n]), 2, 2048) for n in names]
    d = []
    for i in range(1, len(fps)):
        d += [1 - x for x in DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i])]
    cl = Butina.ClusterData(d, len(fps), 0.35, isDistData=True)
    cl = sorted(cl, key=len, reverse=True)
    out, cur = [], []
    for c in cl:                      # pack clusters into bins by chemistry only
        cur += list(c)
        if len(cur) >= MIN_BIN: out.append(cur); cur = []
    if cur and out: out[-1] += cur
    return [[names[i] for i in b] for b in out]

def panels(tag, lab, smi, vina, other, gnina):
    names = sorted(set(vina) & set(gnina) & set(lab) & set(smi) & set(other))
    rows = []
    rng = np.random.default_rng(SEED)
    perm = dict(zip(names, rng.permutation([gnina[n] for n in names])))
    for b in bins_for(names, smi):
        y = np.array([lab[n] for n in b])
        if not (0 < y.sum() < len(y)): continue
        rows.append(dict(tag=tag, n=len(b),
                         vina=roc_auc_score(y, [-vina[n] for n in b]),
                         gnina=roc_auc_score(y, [gnina[n] for n in b]),
                         ctrl=roc_auc_score(y, [-other[n] for n in b]),
                         perm=roc_auc_score(y, [perm[n] for n in b])))
    return rows

def slope(x, y, nboot=10000):
    x = np.asarray(x); y = np.asarray(y)
    b = np.polyfit(x, y, 1)[0]
    rng = np.random.default_rng(SEED); bs = []
    for _ in range(nboot):
        i = rng.integers(0, len(x), len(x))
        if np.ptp(x[i]) > 1e-9: bs.append(np.polyfit(x[i], y[i], 1)[0])
    return b, float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))

rows = []
# --- Mpro
lab = {c["name"]: c["label"] for m in glob.glob(f"{K}/bundle_CHEMBL4523582_7VU6*/shard*/manifest.json")
       for c in json.load(open(m))["compounds"]}
smi = json.load(open(f"{K}/mpro_smiles.json"))
gn = {strip(k): v["cnn_affinity"] for k, v in json.load(open(f"{K}/gnina_rescore.json")).items()}
vi = {strip(k): v for k, v in json.load(open(f"{K}/exh4_scores_merged.json")).items()}
s = {json.load(open(f))["seed"]: json.load(open(f))["scores"]
     for f in glob.glob("analysis/seed_floor/results/seed_exh4_s*_shard0of1.json")}
rows += panels("Mpro", lab, smi, vi, s[2], gn)
# --- Factor Xa
lab2 = {c["name"]: c["label"] for m in glob.glob(f"{K}/bundle_CHEMBL244_PROT/shard*/manifest.json")
        for c in json.load(open(m))["compounds"]}
smi2 = json.load(open(f"{K}/fxa_smiles.json"))
gn2 = {k: v["cnn_affinity"] for k, v in json.load(open("analysis/replication_fxa/gnina_fxa.json"))["scores"].items()}
s2 = {json.load(open(f))["seed"]: json.load(open(f))["scores"]
      for f in glob.glob("analysis/seed_floor/results_fxa/seed_exh4_s*_shard0of1.json")}
rows += panels("FXa", lab2, smi2, s2[1], s2[2], gn2)

for t in ("Mpro", "FXa"):
    r = [x for x in rows if x["tag"] == t]
    print(f"{t}: {len(r)} sub-panels, sizes {min(x['n'] for x in r)}-{max(x['n'] for x in r)}, "
          f"Vina AUROC {min(x['vina'] for x in r):.3f}-{max(x['vina'] for x in r):.3f}")
x = [r["vina"] for r in rows]
print(f"\npooled: {len(rows)} sub-panels, Vina AUROC spread {min(x):.3f}-{max(x):.3f}")

bc, lc, hc = slope(x, [r["ctrl"] for r in rows])
print(f"\nPOSITIVE CONTROL  Vina(seed 2) on Vina : slope {bc:+.3f} [{lc:+.3f}, {hc:+.3f}]")
print(f"  {'PASS -- the design can detect tracking' if bc >= 0.70 else '*** VOID: regression dilution, nothing is read ***'}")
bp, lp, hp = slope(x, [r["perm"] for r in rows])
print(f"NEGATIVE CONTROL  permuted gnina       : slope {bp:+.3f} [{lp:+.3f}, {hp:+.3f}]  "
      f"mean AUROC {np.mean([r['perm'] for r in rows]):.3f}")
bg, lg, hg = slope(x, [r["gnina"] for r in rows])
print(f"\nTEST              gnina on Vina        : slope {bg:+.3f} [{lg:+.3f}, {hg:+.3f}]")
verdict = ("VOID" if bc < 0.70 else
           "SUPPORTS broken-baseline" if hg < 0.5 else
           "REFUTES broken-baseline" if lg > 0.5 else "INCONCLUSIVE -> cut from the manuscript")
print(f"PRE-REGISTERED READING: {verdict}")
json.dump(dict(n_subpanels=len(rows), control_slope=[bc, lc, hc], perm_slope=[bp, lp, hp],
               gnina_slope=[bg, lg, hg], verdict=verdict, rows=rows),
          open("analysis/replication_fxa/MECHANISM_TEST.json", "w"), indent=1, default=float)
