"""Analysis for the Boltz-2 Mpro run. WRITTEN BEFORE THE DATA EXISTS.

Every rule here is transcribed from PREREGISTRATION.md (9e3083d02 + amendments), which was
committed before any Boltz-2 output. Writing the analysis ahead of the results is the strongest
form of pre-registration: there is no opportunity to shape the test around what came back.

Usage:  python3 analyse_boltz2.py boltz2_scores.json [boltz2_failures.json]
"""
import json, sys, glob, warnings
warnings.filterwarnings("ignore")
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler

rng = np.random.default_rng(0)
K = "analysis/retrospective_benchmark/kaggle"
DESC_BAR = 0.7653          # pre-registered bar: the descriptor baseline
MARGIN = 0.04              # pre-registered required margin
sp = lambda k: k[6:] if k.startswith("bench_") else k

scores = json.load(open(sys.argv[1]))
fails = json.load(open(sys.argv[2])) if len(sys.argv) > 2 else {}
lab = {}
for m in glob.glob(f"{K}/bundle_CHEMBL4523582_7VU6*/shard*/manifest.json"):
    for c in json.load(open(m))["compounds"]:
        lab[c["name"]] = c["label"]
smi = json.load(open(f"{K}/mpro_smiles.json"))

# ---------------------------------------------------------------- dropout control (can VOID)
got = {sp(k) for k in scores}
miss = sorted(set(lab) - got)
if miss and got:
    am = np.mean([lab[n] for n in miss])
    ag = np.mean([lab[n] for n in got if n in lab])
    print(f"DROPOUT CONTROL  scored {len(got)}  missing {len(miss)}")
    print(f"  active fraction  missing {am:.1%}  scored {ag:.1%}  diff {am-ag:+.1%}")
    if abs(am - ag) > 0.10:
        print("  *** VOID: dropout is class-selected beyond 10pp (pre-registered) ***")
else:
    print(f"DROPOUT CONTROL  scored {len(got)}  missing 0 -- clean")

rows = []
for n in sorted(got & set(lab) & set(smi)):
    mol = Chem.MolFromSmiles(smi[n])
    if mol is None:
        continue
    v = scores.get(n) or scores.get("bench_" + n)
    pb, pv = v.get("affinity_probability_binary"), v.get("affinity_pred_value")
    if pb is None and pv is None:
        continue
    rows.append(([Descriptors.MolWt(mol), Crippen.MolLogP(mol), Descriptors.NumHDonors(mol),
                  Descriptors.NumHAcceptors(mol), Descriptors.NumRotatableBonds(mol),
                  Descriptors.TPSA(mol), Chem.GetFormalCharge(mol)], lab[n], pb, pv))
D = np.array([r[0] for r in rows]); Y = np.array([r[1] for r in rows])
PB = np.array([np.nan if r[2] is None else r[2] for r in rows], dtype=float)
PV = np.array([np.nan if r[3] is None else r[3] for r in rows], dtype=float)
print(f"\nn={len(Y)}  actives={int(Y.sum())} ({Y.mean():.1%})")

folds = list(StratifiedKFold(5, shuffle=True, random_state=0).split(D, Y))
def oof(X):
    p = np.zeros(len(Y))
    for tr, te in folds:
        s = StandardScaler().fit(X[tr])
        p[te] = LogisticRegression(max_iter=3000).fit(s.transform(X[tr]), Y[tr]
                ).predict_proba(s.transform(X[te]))[:, 1]
    return p
def boot(a, b, n=4000):
    d = []
    for _ in range(n):
        i = rng.integers(0, len(Y), len(Y))
        if len(set(Y[i])) > 1:
            d.append(roc_auc_score(Y[i], a[i]) - roc_auc_score(Y[i], b[i]))
    return np.percentile(d, [2.5, 97.5])

# ------------------------------------------------- positive control (must pass BEFORE reading)
desc = oof(D); a_desc = roc_auc_score(Y, desc)
print(f"\nPOSITIVE CONTROL  descriptors {a_desc:.4f}  (expect {DESC_BAR} +/- 0.005)")
if abs(a_desc - DESC_BAR) > 0.005:
    print("  *** pipeline does not reproduce a number we already know -- no null from it is "
          "trustworthy (pre-registered) ***")

perm = Y.copy(); rng.shuffle(perm)
a_perm = roc_auc_score(perm, desc)
print(f"PERMUTATION       shuffled labels {a_perm:.4f}  (expect [0.45, 0.55])")
if not 0.45 <= a_perm <= 0.55:
    print("  *** VOID: evaluation is leaking (pre-registered) ***")

# -------------------------------------------------------------------------- pre-registered rules
print("\n--- PRE-REGISTERED RULES ---")
for name, arr, primary in (("affinity_probability_binary", PB, True),
                           ("affinity_pred_value", PV, False)):
    if np.all(np.isnan(arr)):
        print(f"{name}: absent"); continue
    ok = ~np.isnan(arr)
    a_raw = roc_auc_score(Y[ok], arr[ok])
    tag = "PRIMARY" if primary else "secondary"
    print(f"\n[{tag}] {name}: AUROC {a_raw:.4f}  (Vina 0.4530, descriptors {a_desc:.4f})")
    lo, hi = boot(arr if not np.isnan(arr).any() else np.nan_to_num(arr, nan=np.nanmedian(arr)),
                  desc)
    d = a_raw - a_desc
    verdict = "SUPPORTED" if (d > MARGIN and lo > 0) else "NOT DEMONSTRATED"
    print(f"  RULE 1  boltz - descriptors {d:+.4f}  95% CI [{lo:+.4f}, {hi:+.4f}]  -> {verdict}")
    comb = oof(np.column_stack([D, np.nan_to_num(arr, nan=np.nanmedian(arr))]))
    a_comb = roc_auc_score(Y, comb); lo3, hi3 = boot(comb, desc)
    v3 = "SUPPORTED" if (a_comb - a_desc > MARGIN and lo3 > 0) else "NOT DEMONSTRATED"
    print(f"  RULE 3  desc+boltz {a_comb:.4f}  delta {a_comb-a_desc:+.4f}  "
          f"CI [{lo3:+.4f}, {hi3:+.4f}]  -> {v3}")
print("\nRule 2 (beat Vina 0.4530 by >0.04) is the weaker literature claim; read from the "
      "AUROC printed above.")
