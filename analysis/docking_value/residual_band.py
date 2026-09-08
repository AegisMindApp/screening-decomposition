#!/usr/bin/env python3
"""Recompute every docking residual under ONE procedure.

The 0.507-0.559 band in the manuscript mixes two residualisation procedures: MARGINAL_VALUE.md
used out-of-fold logistic regression, SIZE_BIAS.md used out-of-fold gradient boosting. Comparing
Boltz-2's gradient-boosted residual against a band built partly from logistic residuals is not
like for like. This recomputes all of them both ways so the comparison is internally consistent
and the choice of procedure is visible rather than buried.

Residual AUROC = AUROC of (score - E[score | seven descriptors]), predictions out-of-fold.

POSITIVE CONTROL: must reproduce SIZE_BIAS.md's gradient-boosted values (Vina repaired 0.5382,
Boltz-2 0.6418) before any new number is readable.
"""
import json, numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
from sklearn.model_selection import KFold
from sklearn.linear_model import LinearRegression
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.metrics import roc_auc_score

K = "analysis/retrospective_benchmark/kaggle"
SEED = 20260906
SEEDS = [20260906 + i for i in range(12)]

def descriptors(smi):
    m = Chem.MolFromSmiles(smi)
    if m is None: return None
    return [Descriptors.MolWt(m), Descriptors.MolLogP(m),
            rdMolDescriptors.CalcNumHBD(m), rdMolDescriptors.CalcNumHBA(m),
            rdMolDescriptors.CalcNumRotatableBonds(m), rdMolDescriptors.CalcTPSA(m),
            Chem.GetFormalCharge(m)]

def load(target):
    smi = json.load(open(f"{K}/{target}_smiles.json"))
    if target == "mpro":
        lab = json.load(open("analysis/docking_value/mpro_labels.json"))
    else:
        lab = {c["name"]: c["label"]
               for c in json.load(open("analysis/boltz2/bundle_fxa/manifest.json"))["compounds"]}
    return smi, lab

def residual_auroc(names, scores, labels, X, model, seed=SEED):
    """Out-of-fold residual: fit descriptors -> score on 4 folds, residualise the 5th."""
    X = np.asarray(X, float); y = np.asarray(scores, float); lab = np.asarray(labels, int)
    resid = np.zeros_like(y)
    for tr, te in KFold(5, shuffle=True, random_state=seed).split(X):
        m = (GradientBoostingRegressor(random_state=seed) if model == "gb"
             else LinearRegression())
        m.fit(X[tr], y[tr])
        resid[te] = y[te] - m.predict(X[te])
    # scores are "more negative = better" for Vina-family, handled by the caller via sign
    return roc_auc_score(lab, resid), roc_auc_score(lab, y)

def arm(target, name, getscore, sign):
    smi, lab = load(target)
    rows = []
    for n, s in smi.items():
        v = getscore(n)
        if v is None or n not in lab: continue
        d = descriptors(s)
        if d is None: continue
        rows.append((n, sign * v, lab[n], d))
    if len(rows) < 100:
        raise SystemExit(f"{name}: only {len(rows)} rows joined - check key formats, "
                         f"an empty join is not a result")
    names = [r[0] for r in rows]; sc = [r[1] for r in rows]
    y = [r[2] for r in rows];     X = [r[3] for r in rows]
    out = {}
    for model in ("gb", "lin"):
        # Residual AUROC moves with the fold assignment. Sweep seeds rather than quoting one.
        vals = [residual_auroc(names, sc, y, X, model, s)[0] for s in SEEDS]
        r, raw = residual_auroc(names, sc, y, X, model)
        out[model] = r
        out[model + "_mean"] = float(np.mean(vals))
        out[model + "_lo"] = float(np.min(vals)); out[model + "_hi"] = float(np.max(vals))
    out["raw"] = raw; out["n"] = len(rows); out["target"] = target
    return name, out

# Every score cache uses a "bench_" key prefix except the Boltz-2 output. Normalise all of
# them: an unnormalised join silently yields zero rows, which looks like a missing arm rather
# than a bug.
def norm(d):
    return {(k[6:] if k.startswith("bench_") else k): v for k, v in d.items()}

vina_prot = norm(json.load(open(f"{K}/protonated_scores_merged.json")))
vina_exh4 = norm(json.load(open(f"{K}/exh4_scores_merged.json")))
gn        = norm(json.load(open(f"{K}/gnina_rescore.json")))
fxa_vina  = norm(json.load(open(f"{K}/fxa_protonated_all.json")))
bz        = norm(json.load(open("analysis/boltz2/results/boltz2_scores.json")))

ARMS = [
    arm("mpro", "Vina, repaired receptor (Mpro)",  lambda n: vina_prot.get(n), -1),
    arm("mpro", "Vina, raw receptor (Mpro)",       lambda n: vina_exh4.get(n), -1),
    arm("mpro", "gnina Vina term (Mpro)",          lambda n: gn.get(n, {}).get("vina_affinity"), -1),
    arm("mpro", "gnina CNNaffinity (Mpro)",        lambda n: gn.get(n, {}).get("cnn_affinity"), +1),
    arm("mpro", "gnina CNNscore (Mpro)",           lambda n: gn.get(n, {}).get("cnn_score"), +1),
    arm("fxa",  "Vina, repaired receptor (FXa)",   lambda n: fxa_vina.get(n), -1),
]
BZ = arm("mpro", "Boltz-2 prob_binary (Mpro)",
         lambda n: bz.get(n, {}).get("affinity_probability_binary"), +1)

# FlashBind keys are "{prot_id}_{ligand_id}"; every other cache is keyed on the ligand alone.
import os
FB_PATH = "analysis/flashbind/results/flashbind_scores_751.json"
FB = None
if os.path.exists(FB_PATH):
    fb = {k.split("_", 1)[1]: v for k, v in
          json.load(open(FB_PATH))["ensemble"].items()}
    FB = arm("mpro", "FlashBind binary (Mpro)", lambda n: fb.get(n), +1)

print(f"{'score':38s} {'n':>4s} {'raw':>8s} {'resid GB (12-seed range)':>30s}")
for nm, o in ARMS + [BZ] + ([FB] if FB else []):
    print(f"{nm:38s} {o['n']:4d} {o['raw']:8.4f} "
          f"{o['gb_mean']:10.4f}  [{o['gb_lo']:.4f}, {o['gb_hi']:.4f}]")

ctrl_vina = dict(ARMS)["Vina, repaired receptor (Mpro)"]["gb"]
ctrl_bz   = BZ[1]["gb"]
print(f"\nPOSITIVE CONTROL vs SIZE_BIAS.md (gradient-boosted)")
print(f"  Vina repaired  {ctrl_vina:.4f}   published 0.5382   delta {abs(ctrl_vina-0.5382):.4f}")
print(f"  Boltz-2        {ctrl_bz:.4f}   published 0.6418   delta {abs(ctrl_bz-0.6418):.4f}")
ok = abs(ctrl_vina-0.5382) < 0.02 and abs(ctrl_bz-0.6418) < 0.02
print(f"  CONTROL: {'PASS' if ok else 'FAIL'}")

dock_gb  = [o["gb_mean"]  for nm, o in ARMS]
dock_lin = [o["lin_mean"] for nm, o in ARMS]
dock_gb_hi = max(o["gb_hi"] for nm, o in ARMS)
print(f"\nWorst case across seeds: highest docking residual {dock_gb_hi:.4f}, "
      f"lowest Boltz-2 residual {BZ[1]['gb_lo']:.4f}, "
      f"separated: {BZ[1]['gb_lo'] > dock_gb_hi}")
print(f"\nDOCKING BAND, single procedure")
print(f"  gradient-boosted : {min(dock_gb):.4f} - {max(dock_gb):.4f}")
print(f"  linear           : {min(dock_lin):.4f} - {max(dock_lin):.4f}")
print(f"  Boltz-2          : GB {BZ[1]['gb']:.4f}   lin {BZ[1]['lin']:.4f}")
print(f"  Boltz-2 outside the GB band?  {BZ[1]['gb_mean'] > max(dock_gb)}")
print(f"  Boltz-2 outside the lin band? {BZ[1]['lin_mean'] > max(dock_lin)}")

json.dump({"seed": SEED, "control_pass": bool(ok),
           "arms": {nm: o for nm, o in ARMS + [BZ] + ([FB] if FB else [])},
           "band_gb": [min(dock_gb), max(dock_gb)],
           "band_lin": [min(dock_lin), max(dock_lin)],
           "boltz2": {"gb": BZ[1]["gb_mean"], "lin": BZ[1]["lin_mean"],
                      "gb_range": [BZ[1]["gb_lo"], BZ[1]["gb_hi"]]},
           "seeds": SEEDS},
          open("analysis/docking_value/RESIDUAL_BAND.json", "w"), indent=1)

if FB:
    print(f"\nFLASHBIND vs the pre-registered reading rules")
    print(f"  raw AUROC          {FB[1]['raw']:.4f}   (Boltz-2 {BZ[1]['raw']:.4f})")
    print(f"  residual GB        {FB[1]['gb_mean']:.4f}  "
          f"[{FB[1]['gb_lo']:.4f}, {FB[1]['gb_hi']:.4f}]  n={FB[1]['n']}")
    print(f"  docking band top   {max(dock_gb):.4f}   (worst-case seed {dock_gb_hi:.4f})")
    print(f"  Boltz-2 residual   {BZ[1]['gb_mean']:.4f}  "
          f"[{BZ[1]['gb_lo']:.4f}, {BZ[1]['gb_hi']:.4f}]")
    above = FB[1]["gb_lo"] > 0.573
    within = 0.494 <= FB[1]["gb_mean"] <= 0.573
    print(f"  > 0.573 with the 12-seed range clear of it : {above}")
    print(f"  inside the pre-registered band [0.494, 0.573]: {within}")
