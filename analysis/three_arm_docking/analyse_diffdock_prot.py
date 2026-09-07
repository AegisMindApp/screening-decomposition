#!/usr/bin/env python3
"""DiffDock vs Vina poses, both scored by Vina on the REPAIRED receptor.

Rules in PREREGISTRATION_DIFFDOCK_PROT.md (ce636ea0a) + Amendment 1 (coverage floor 75%,
dropout-bias control). Gated against the PROTOCOL floor 0.0393: this arm changes placement.
"""
import json, glob, os
import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors
from sklearn.metrics import roc_auc_score
RDLogger.DisableLog("rdApp.*")

K = "analysis/retrospective_benchmark/kaggle"
FLOOR_PROTOCOL, SEED, NBOOT = 0.0393, 20260907, 10000
strip = lambda n: n[6:] if n.startswith("bench_") else n

dd = {os.path.basename(p)[:-5]: json.load(open(p))
      for p in glob.glob("analysis/three_arm_docking/diffdock_prot_scored/*.json")}
dd_ok = {k: v["score"] for k, v in dd.items() if v.get("score") is not None}
vina = {strip(k): v for k, v in json.load(open(f"{K}/protonated_scores_merged.json")).items()}
lab = {}
for m in glob.glob(f"{K}/bundle_CHEMBL4523582_PROTONATED/shard*/manifest.json"):
    for c in json.load(open(m))["compounds"]: lab[c["name"]] = c["label"]
smi = json.load(open(f"{K}/mpro_smiles.json"))

cov = len(dd_ok) / len(dd)
print(f"coverage {len(dd_ok)}/{len(dd)} = {cov:.1%}  (amended floor 75%)  "
      f"{'PASS' if cov >= 0.75 else 'VOID'}")
assert cov >= 0.75, "below the amended floor - run is void"
print(f"published arm coverage for comparison: 604/751 = 80.4%")

# ---- CONTROL: is the exclusion label-biased? --------------------------------------------
excl = [k for k in dd if k not in dd_ok and k in lab]
incl = [k for k in dd_ok if k in lab]
af_excl = np.mean([lab[k] for k in excl]); af_panel = np.mean(list(lab.values()))
mw = lambda ks: np.mean([Descriptors.MolWt(Chem.MolFromSmiles(smi[k]))
                         for k in ks if k in smi and Chem.MolFromSmiles(smi[k])])
print(f"\nDROPOUT CONTROL")
print(f"  excluded n={len(excl)}  active fraction {af_excl:.3f} vs panel {af_panel:.3f} "
      f"(|diff| {abs(af_excl-af_panel):.3f}, tol 0.10)  "
      f"{'PASS' if abs(af_excl-af_panel) <= 0.10 else 'FAIL - label-biased exclusion'}")
print(f"  mean MW excluded {mw(excl):.1f} vs included {mw(incl):.1f}")
assert abs(af_excl - af_panel) <= 0.10, "exclusion is label-biased - arm void"

# ---- CONTROL: comparator reproduces its known value on the shared subset -----------------
both = sorted(set(dd_ok) & set(vina) & set(lab))
y = np.array([lab[n] for n in both])
a_vina_full = roc_auc_score([lab[n] for n in vina if n in lab],
                            [-vina[n] for n in vina if n in lab])
a_vina = roc_auc_score(y, [-vina[n] for n in both])
a_dd = roc_auc_score(y, [-dd_ok[n] for n in both])
print(f"\nCOMPARATOR CONTROL")
print(f"  Vina poses, full panel {a_vina_full:.4f} vs published 0.4530 "
      f"(tol 0.005) {'PASS' if abs(a_vina_full-0.4530)<=0.005 else 'FAIL'}")

print(f"\nPRIMARY (paired on the {len(both)} shared compounds, repaired receptor)")
print(f"  Vina poses     {a_vina:.4f}")
print(f"  DiffDock poses {a_dd:.4f}")
d0 = a_dd - a_vina
rng = np.random.default_rng(SEED); N = len(y); out = []
dv = np.array([-vina[n] for n in both]); dm = np.array([-dd_ok[n] for n in both])
for _ in range(NBOOT):
    i = rng.integers(0, N, N)
    if len(set(y[i].tolist())) < 2: continue
    out.append(roc_auc_score(y[i], dm[i]) - roc_auc_score(y[i], dv[i]))
o = np.sort(out); lo, hi = float(o[int(.025*len(o))]), float(o[int(.975*len(o))])
print(f"  delta {d0:+.4f} [{lo:+.4f}, {hi:+.4f}]  vs protocol floor {FLOOR_PROTOCOL}")
resolvable = abs(d0) > FLOOR_PROTOCOL and (lo > 0 or hi < 0)
print(f"\n=== VERDICT: {'pose source MATTERS' if resolvable else 'null STANDS'} ===")
print(f"  published (defective receptor): +0.0173 [-0.0269, +0.0636]")
json.dump({"coverage": cov, "n_shared": len(both), "auroc_vina": a_vina, "auroc_diffdock": a_dd,
           "delta": d0, "ci": [lo, hi], "floor": FLOOR_PROTOCOL, "resolvable": bool(resolvable),
           "excluded_active_frac": float(af_excl), "panel_active_frac": float(af_panel),
           "mw_excluded": float(mw(excl)), "mw_included": float(mw(incl))},
          open("analysis/three_arm_docking/DIFFDOCK_PROT.json", "w"), indent=1)
json.dump(dd_ok, open("analysis/three_arm_docking/diffdock_prot_scores.json", "w"), indent=0)
print("\nwrote DIFFDOCK_PROT.json and diffdock_prot_scores.json")
