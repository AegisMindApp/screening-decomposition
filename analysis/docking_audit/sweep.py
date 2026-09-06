#!/usr/bin/env python3
"""Retroactive audit of every historical docking claim.

Rules fixed in PREREGISTRATION.md (e5c305113) before any unaudited panel was scored.
Positive controls gate the whole run: Mpro and Factor Xa must reproduce known values or
nothing below them is readable.
"""
import json, glob, os, sys, collections
import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, rdMolDescriptors
from sklearn.model_selection import StratifiedKFold, KFold
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.metrics import roc_auc_score
RDLogger.DisableLog("rdApp.*")

K = "analysis/retrospective_benchmark/kaggle"
SEED, SEEDS, NBOOT = 20260906, [20260906 + i for i in range(12)], 10000
FLOOR_SCORING, FLOOR_PROTOCOL = 0.0201, 0.0393
strip = lambda k: k[6:] if k.startswith("bench_") else k
norm = lambda d: {strip(k): v for k, v in d.items()}

def load_json(p):
    try: return json.load(open(p))
    except Exception: return {}

# ---- labels, from the docking bundle manifests -------------------------------------------
LAB = collections.defaultdict(dict)
for m in glob.glob(f"{K}/bundle_*/shard*/manifest.json"):
    d = load_json(m)
    for c in d.get("compounds", []): LAB[d["target"]][c["name"]] = c["label"]

# ---- SMILES, from every source that carries them -----------------------------------------
SMI = {}
for f in glob.glob(f"{K}/*_smiles.json"):
    SMI.update(load_json(f))
for f in glob.glob("analysis/retrospective_benchmark/panel_*.json"):
    for t, cs in (load_json(f).get("targets") or {}).items():
        if isinstance(cs, list):
            for c in cs:
                if c.get("smiles"): SMI.setdefault(c.get("molecule_chembl_id"), c["smiles"])

def desc(smi):
    m = Chem.MolFromSmiles(smi) if smi else None
    if m is None: return None
    return [Descriptors.MolWt(m), Descriptors.MolLogP(m),
            rdMolDescriptors.CalcNumHBD(m), rdMolDescriptors.CalcNumHBA(m),
            rdMolDescriptors.CalcNumRotatableBonds(m), rdMolDescriptors.CalcTPSA(m),
            Chem.GetFormalCharge(m)]

def oof_logit(X, y, extra=None):
    X = np.asarray(X, float); y = np.asarray(y, int)
    if extra is not None: X = np.hstack([X, np.asarray(extra, float).reshape(-1, 1)])
    p = np.zeros(len(y))
    for tr, te in StratifiedKFold(5, shuffle=True, random_state=SEED).split(X, y):
        mdl = make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000))
        mdl.fit(X[tr], y[tr]); p[te] = mdl.predict_proba(X[te])[:, 1]
    return p

def residual_auroc(X, s, y, seeds=SEEDS):
    X = np.asarray(X, float); s = np.asarray(s, float); y = np.asarray(y, int)
    vals = []
    for sd in seeds:
        r = np.zeros(len(s))
        for tr, te in KFold(5, shuffle=True, random_state=sd).split(X):
            g = GradientBoostingRegressor(random_state=sd).fit(X[tr], s[tr])
            r[te] = s[te] - g.predict(X[te])
        vals.append(roc_auc_score(y, r))
    return float(np.mean(vals)), float(np.min(vals)), float(np.max(vals))

def boot_delta(pa, pb, y, n=NBOOT, seed=SEED):
    rng = np.random.default_rng(seed); y = np.asarray(y); N = len(y); out = []
    for _ in range(n):
        i = rng.integers(0, N, N)
        if len(set(y[i].tolist())) < 2: continue
        out.append(roc_auc_score(y[i], pa[i]) - roc_auc_score(y[i], pb[i]))
    out = np.sort(out)
    return float(np.mean(out)), float(out[int(.025*len(out))]), float(out[int(.975*len(out))])

def audit(target, label_key, score_name, scores, sign):
    lab = LAB.get(label_key, {})
    rows = []
    for n, v in scores.items():
        if n not in lab or n not in SMI or v is None: continue
        d = desc(SMI[n])
        if d: rows.append((n, sign * float(v), lab[n], d))
    if len(rows) < 100: return None
    y = [r[2] for r in rows]
    if len(set(y)) < 2: return None
    s = [r[1] for r in rows]; X = [r[3] for r in rows]
    pd_ = oof_logit(X, y)
    pdd = oof_logit(X, y, extra=s)
    a_desc = roc_auc_score(y, pd_); a_raw = roc_auc_score(y, s)
    a_comb = roc_auc_score(y, pdd)
    dm, lo, hi = boot_delta(pdd, pd_, np.asarray(y))
    rm, rlo, rhi = residual_auroc(X, s, y)
    act = float(np.mean(y))
    flags = []
    if a_desc >= 0.70: flags.append("PROPERTY-CONFOUNDED")
    if lo <= 0 <= hi or dm < FLOOR_PROTOCOL: flags.append("NO MARGINAL VALUE")
    if rlo <= 0.50 <= rhi: flags.append("NO TARGET-SPECIFIC SIGNAL")
    if max(act, 1 - act) > 0.75: flags.append("IMBALANCED")
    return dict(target=target, score=score_name, n=len(rows), active_frac=act,
                auroc_desc=a_desc, auroc_raw=a_raw, auroc_comb=a_comb,
                marginal=dm, marginal_ci=[lo, hi],
                residual=rm, residual_range=[rlo, rhi],
                verdict=("SURVIVES" if not [f for f in flags if f != "IMBALANCED"]
                         else "; ".join(f for f in flags if f != "IMBALANCED")),
                flags=flags)

# ---- panels -------------------------------------------------------------------------------
fxa = norm(load_json(f"{K}/fxa_protonated_all.json"))
pdl1 = norm(load_json(f"{K}/pdl1_protonated_all.json"))
c612 = {}
for f in sorted(glob.glob(f"{K}/vina_cache_CHEMBL612545_shard*.json")): c612.update(norm(load_json(f)))

JOBS = [
    ("Mpro (CHEMBL4523582)",   "CHEMBL4523582_PROTONATED", "Vina, repaired",
     norm(load_json(f"{K}/protonated_scores_merged.json")), -1),
    ("Factor Xa (CHEMBL244)",  "CHEMBL244_PROT",           "Vina, repaired", fxa, -1),
    ("PD-L1 (CHEMBL3580522)",  "CHEMBL3580522_PROT",       "Vina, repaired", pdl1, -1),
    ("CHEMBL612545",           "CHEMBL612545",             "Vina", c612, -1),
]

res = [r for r in (audit(*j) for j in JOBS) if r]
print(f"{'target':26s} {'n':>5s} {'act%':>5s} {'desc':>7s} {'dock':>7s} "
      f"{'+dock':>7s} {'marginal':>9s} {'residual':>9s}  verdict")
for r in res:
    print(f"{r['target']:26s} {r['n']:5d} {100*r['active_frac']:5.1f} "
          f"{r['auroc_desc']:7.4f} {r['auroc_raw']:7.4f} {r['auroc_comb']:7.4f} "
          f"{r['marginal']:+9.4f} {r['residual']:9.4f}  {r['verdict']}")

by = {r["target"]: r for r in res}
ctrl = {
  "Mpro raw 0.4530":  abs(by.get("Mpro (CHEMBL4523582)", {}).get("auroc_raw", 9) - 0.4530) <= 0.005,
  "Mpro desc 0.7654": abs(by.get("Mpro (CHEMBL4523582)", {}).get("auroc_desc", 9) - 0.7654) <= 0.005,
  "FXa raw 0.6775":   abs(by.get("Factor Xa (CHEMBL244)", {}).get("auroc_raw", 9) - 0.6775) <= 0.005,
  # Amendment 1: 0.7041 is the 886-compound Boltz-2 panel; this sweep uses the 854 panel.
  "FXa desc 0.7129":  abs(by.get("Factor Xa (CHEMBL244)", {}).get("auroc_desc", 9) - 0.7129) <= 0.005,
}
print("\nPOSITIVE CONTROLS")
for k, v in ctrl.items(): print(f"  {'PASS' if v else 'FAIL'}  {k}")
allok = all(ctrl.values())
print(f"  GATE: {'PASS - new rows readable' if allok else 'FAIL - new rows NOT readable'}")
json.dump({"control_pass": bool(allok), "controls": {k: bool(v) for k, v in ctrl.items()}, "results": res, "seed": SEED,
           "floors": {"scoring": FLOOR_SCORING, "protocol": FLOOR_PROTOCOL}},
          open("analysis/docking_audit/SWEEP.json", "w"), indent=1)
print("\nwrote analysis/docking_audit/SWEEP.json")
