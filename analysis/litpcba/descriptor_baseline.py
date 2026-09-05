"""Seven free descriptors vs LIT-PCBA, the field's standard 'unbiased' benchmark.

We measured docking losing to these descriptors on our own two panels. Those panels could be
confounded. LIT-PCBA was explicitly built to remove property bias -- and the AVE_unbiased
variant debiases further -- so if descriptors still separate actives here, the finding is not
an artefact of our panels. If they do not, our panels were confounded and we should say so.

Out-of-fold logistic regression, same seven descriptors as every other run in this project.
Inactives subsampled (seeded) for tractability: AUROC is unbiased under random subsampling of
one class, and the subsample size is reported so it can be checked.
"""
import sys, json, random, warnings, os
warnings.filterwarnings("ignore")
import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, Crippen
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler
RDLogger.DisableLog("rdApp.*")

ROOT = sys.argv[1] if len(sys.argv) > 1 else "analysis/litpcba/full"
MAX_INACT = int(os.environ.get("MAX_INACT", "25000"))
rng = random.Random(0)

def feats(smi):
    m = Chem.MolFromSmiles(smi)
    if m is None: return None
    return [Descriptors.MolWt(m), Crippen.MolLogP(m), Descriptors.NumHDonors(m),
            Descriptors.NumHAcceptors(m), Descriptors.NumRotatableBonds(m),
            Descriptors.TPSA(m), Chem.GetFormalCharge(m)]

def load(path, cap=None):
    rows = [l.split()[0] for l in open(path) if l.strip()]
    if cap and len(rows) > cap: rows = rng.sample(rows, cap)
    return rows

def find_sets(d):
    """LIT-PCBA ships two layouts: full_data has actives.smi/inactives.smi, AVE_unbiased has
    active_T/V.smi + inactive_T/V.smi (train/validation). Returning None for an unrecognised
    layout used to make this script skip every target and exit 0 with an empty table -- a
    silent failure that looked like a clean run."""
    if os.path.exists(f"{d}/actives.smi") and os.path.exists(f"{d}/inactives.smi"):
        return ("cv", f"{d}/actives.smi", f"{d}/inactives.smi", None, None)
    if all(os.path.exists(f"{d}/{x}.smi") for x in ("active_T","active_V","inactive_T","inactive_V")):
        return ("tv", f"{d}/active_T.smi", f"{d}/inactive_T.smi",
                      f"{d}/active_V.smi", f"{d}/inactive_V.smi")
    return None

out = {}
targets = sorted(d for d in os.listdir(ROOT) if os.path.isdir(os.path.join(ROOT, d)))
_usable = [t for t in targets if find_sets(f"{ROOT}/{t}")]
if not _usable:
    sys.exit(f"FATAL: no usable target under {ROOT}. Found {len(targets)} dirs but none carry a "
             f"recognised layout (actives.smi/inactives.smi, or active_T/V + inactive_T/V). "
             f"Exiting loudly rather than printing an empty table.")
print(f"  layout: {find_sets(f'{ROOT}/{_usable[0]}')[0]}   usable targets: {len(_usable)}/{len(targets)}")
print(f"{'target':<12}{'act':>6}{'inact':>8}{'AUROC':>9}{'EF@1%':>8}")
for t in targets:
    fs = find_sets(f"{ROOT}/{t}")
    if not fs: continue
    mode, ap, ip, avp, ivp = fs
    def build(a_path, i_path):
        X, Y = [], []
        for smi, y in [(s,1) for s in load(a_path)] + [(s,0) for s in load(i_path, MAX_INACT)]:
            f = feats(smi)
            if f: X.append(f); Y.append(y)
        return np.array(X), np.array(Y)
    X, Y = build(ap, ip)
    XV, YV = build(avp, ivp) if mode == "tv" else (None, None)
    if Y.sum() < 5 or len(Y) - Y.sum() < 50: 
        print(f"{t:<12}{int(Y.sum()):>6}{len(Y)-int(Y.sum()):>8}      skipped (too few)")
        continue
    if mode == "tv" and YV is not None and YV.sum() >= 3:
        # AVE ships its own train/validation split -- honour it rather than re-folding
        sc = StandardScaler().fit(X)
        p = LogisticRegression(max_iter=3000, class_weight="balanced").fit(
            sc.transform(X), Y).predict_proba(sc.transform(XV))[:, 1]
        Y = YV
    else:
        n_splits = min(5, int(Y.sum()))
        p = np.zeros(len(Y))
        for tr, te in StratifiedKFold(n_splits, shuffle=True, random_state=0).split(X, Y):
            sc = StandardScaler().fit(X[tr])
            p[te] = LogisticRegression(max_iter=3000, class_weight="balanced").fit(
                sc.transform(X[tr]), Y[tr]).predict_proba(sc.transform(X[te]))[:, 1]
    auc = roc_auc_score(Y, p)
    k = max(1, int(0.01 * len(Y)))
    top = np.argsort(-p)[:k]
    ef = (Y[top].sum() / k) / Y.mean()
    out[t] = {"actives": int(Y.sum()), "inactives": int(len(Y) - Y.sum()),
              "auroc": round(float(auc), 4), "ef1pct": round(float(ef), 2)}
    print(f"{t:<12}{int(Y.sum()):>6}{len(Y)-int(Y.sum()):>8}{auc:>9.4f}{ef:>8.2f}")
if out:
    a = [v["auroc"] for v in out.values()]
    print(f"\n  {len(out)} targets   median AUROC {np.median(a):.4f}   "
          f"mean {np.mean(a):.4f}   range {min(a):.4f}-{max(a):.4f}")
    print(f"  above 0.7: {sum(x>0.7 for x in a)}/{len(a)}   above 0.6: {sum(x>0.6 for x in a)}/{len(a)}")
json.dump(out, open(f"analysis/litpcba/descriptor_baseline_{os.path.basename(ROOT)}.json","w"), indent=1)
