#!/usr/bin/env python3
"""Does the descriptor + Boltz-2 combination transfer between targets?

Rules fixed in PREREGISTRATION_TRANSFER.md (fc6816b51) before any transfer estimate existed.

Fit on one target in full, evaluate cold on the other. No fold splitting on the evaluated side,
so there is no fold noise in the primary quantity -- which matters given that fold assignment
has already misled two results in this project.
"""
import json, glob
import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, rdMolDescriptors
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score
RDLogger.DisableLog("rdApp.*")

K = "analysis/retrospective_benchmark/kaggle"
SEED, NBOOT, FLOOR = 20260907, 10000, 0.0201

def desc(s):
    m = Chem.MolFromSmiles(s) if s else None
    if m is None: return None
    return [Descriptors.MolWt(m), Descriptors.MolLogP(m), rdMolDescriptors.CalcNumHBD(m),
            rdMolDescriptors.CalcNumHBA(m), rdMolDescriptors.CalcNumRotatableBonds(m),
            rdMolDescriptors.CalcTPSA(m), Chem.GetFormalCharge(m)]

def panel(name):
    if name == "Mpro":
        lab = {}
        for m in glob.glob(f"{K}/bundle_CHEMBL4523582_PROTONATED/shard*/manifest.json"):
            for c in json.load(open(m))["compounds"]: lab[c["name"]] = c["label"]
        smi = json.load(open(f"{K}/mpro_smiles.json"))
        bz = json.load(open("analysis/boltz2/results/boltz2_scores.json"))
    else:
        lab = {c["name"]: c["label"]
               for c in json.load(open("analysis/boltz2/bundle_fxa/manifest.json"))["compounds"]}
        smi = json.load(open(f"{K}/fxa_smiles.json"))
        bz = {}
        for f in sorted(glob.glob("analysis/boltz2/results_fxa/boltz2_scores_s*.json")):
            bz.update(json.load(open(f)))
    X, b, y = [], [], []
    for n, s in smi.items():
        if n not in lab or n not in bz: continue
        v = bz[n].get("affinity_probability_binary")
        d = desc(s)
        if v is None or d is None: continue
        X.append(d); b.append(float(v)); y.append(lab[n])
    return np.array(X, float), np.array(b, float), np.array(y, int)

def fit_apply(Xtr, ytr, Xte, recipient_scaler):
    sc = StandardScaler().fit(Xte if recipient_scaler else Xtr)
    m = LogisticRegression(max_iter=5000).fit(sc.transform(Xtr), ytr)
    return m.predict_proba(sc.transform(Xte))[:, 1]

def boot(pa, pb, y, n=NBOOT, seed=SEED):
    rng = np.random.default_rng(seed); N = len(y); out = []
    for _ in range(n):
        i = rng.integers(0, N, N)
        if len(set(y[i].tolist())) < 2: continue
        out.append(roc_auc_score(y[i], pa[i]) - roc_auc_score(y[i], pb[i]))
    o = np.sort(out)
    return float(np.mean(o)), float(o[int(.025*len(o))]), float(o[int(.975*len(o))])

def in_target(X, b, y):
    """The published in-target out-of-fold comparison, for reference."""
    A = np.zeros(len(y)); B = np.zeros(len(y))
    XB = np.hstack([X, b.reshape(-1, 1)])
    for tr, te in StratifiedKFold(5, shuffle=True, random_state=SEED).split(X, y):
        scA = StandardScaler().fit(X[tr]); scB = StandardScaler().fit(XB[tr])
        A[te] = LogisticRegression(max_iter=5000).fit(scA.transform(X[tr]), y[tr]) \
                 .predict_proba(scA.transform(X[te]))[:, 1]
        B[te] = LogisticRegression(max_iter=5000).fit(scB.transform(XB[tr]), y[tr]) \
                 .predict_proba(scB.transform(XB[te]))[:, 1]
    return roc_auc_score(y, A), roc_auc_score(y, B), boot(B, A, y)

P = {n: panel(n) for n in ("Mpro", "FXa")}
for n, (X, b, y) in P.items():
    print(f"{n}: n={len(y)} actives={int(y.sum())} ({100*y.mean():.1f}%)")

print("\n=== IN-TARGET (published design, out-of-fold on own labels) ===")
IT = {}
for n, (X, b, y) in P.items():
    a, bb, (d, lo, hi) = in_target(X, b, y); IT[n] = (a, bb, d, lo, hi)
    print(f"  {n:5s} descriptors {a:.4f}  +Boltz-2 {bb:.4f}   delta {d:+.4f} [{lo:+.4f}, {hi:+.4f}]")

print("\n=== TRANSFER (fit on donor in full, applied cold to recipient) ===")
res = {}
for donor, recip in (("Mpro", "FXa"), ("FXa", "Mpro")):
    Xd, bd, yd = P[donor]; Xr, br, yr = P[recip]
    XdB = np.hstack([Xd, bd.reshape(-1, 1)]); XrB = np.hstack([Xr, br.reshape(-1, 1)])
    for regime, rs in (("fully cold", False), ("recipient-standardised", True)):
        pA = fit_apply(Xd, yd, Xr, rs); pB = fit_apply(XdB, yd, XrB, rs)
        aA, aB = roc_auc_score(yr, pA), roc_auc_score(yr, pB)
        d, lo, hi = boot(pB, pA, yr)
        ok = d > FLOOR and lo > 0
        res[f"{donor}->{recip} [{regime}]"] = dict(auroc_A=aA, auroc_B=aB, delta=d, ci=[lo, hi],
                                                   clears=bool(ok))
        print(f"  {donor}->{recip:5s} {regime:22s} A {aA:.4f}  B {aB:.4f}   "
              f"delta {d:+.4f} [{lo:+.4f}, {hi:+.4f}]  {'clears' if ok else 'does not clear'}")

print("\n=== CONTROLS ===")
for n, (X, b, y) in P.items():
    XB = np.hstack([X, b.reshape(-1, 1)])
    p = fit_apply(XB, y, XB, False)
    a = roc_auc_score(y, p)
    print(f"  sanity: {n} model applied to its own training target {a:.4f} "
          f"vs in-target OOF {IT[n][1]:.4f} -> {'PASS' if a >= IT[n][1] - 0.02 else 'FAIL'}")
rng = np.random.default_rng(SEED)
Xd, bd, yd = P["Mpro"]; Xr, br, yr = P["FXa"]
yperm = rng.permutation(yr)
XdB = np.hstack([Xd, bd.reshape(-1, 1)]); XrB = np.hstack([Xr, br.reshape(-1, 1)])
dp, _, _ = boot(fit_apply(XdB, yd, XrB, False), fit_apply(Xd, yd, Xr, False), yperm)
print(f"  permuted recipient labels: delta {dp:+.4f} -> {'PASS' if abs(dp) <= 0.05 else 'FAIL'}")

both = all(v["clears"] for k, v in res.items() if "fully cold" in k)
both_rs = all(v["clears"] for k, v in res.items() if "recipient-standardised" in k)
verdict = "TRANSFERS" if (both or both_rs) else "DOES NOT TRANSFER"
print(f"\n=== VERDICT: {verdict} ===")
print(f"  clears in both directions, fully cold: {both}")
print(f"  clears in both directions, recipient-standardised: {both_rs}")
json.dump({"in_target": {k: list(v) for k, v in IT.items()}, "transfer": res,
           "verdict": verdict, "floor": FLOOR},
          open("analysis/boltz2/TRANSFER.json", "w"), indent=1)
print("\nwrote analysis/boltz2/TRANSFER.json")
