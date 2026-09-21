#!/usr/bin/env python3
"""Chemical-series structure of the Mpro panel, and what it does to the reported intervals.

The JCIM editor named three remedies. This addresses the third — "splits controlling
chemical-series effects" — in the form that actually threatens the manuscript's conclusions.

The obvious reading of that point is that series memorisation inflates absolute AUROC. The
sharper one is about the INTERVALS. `analysis/seed_floor/RESULT.md` records that the
manuscript's intervals "bootstrap compounds". If compounds fall into scaffold series they are
not independent draws, so resampling them independently understates the sampling variance and
every reported CI is too narrow. The three surviving effects are judged by whether their CI
lower bound clears the floor, so an interval that is too narrow can manufacture a survivor.

EXPECTATION, STATED BEFORE RUNNING: a cluster bootstrap that resamples SERIES rather than
compounds should give intervals at least as wide as the compound bootstrap, and the ratio grows
with how concentrated the panel is into few large series. If the widening is small the
manuscript's intervals stand; if it is large, the surviving effects must be re-judged.

    python analysis/seed_floor/series_structure.py
"""
from __future__ import annotations

import glob
import json
from collections import Counter
from pathlib import Path

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from rdkit.DataStructs import BulkTanimotoSimilarity
from rdkit.ML.Cluster import Butina
from scipy import stats

RDLogger.DisableLog("rdApp.*")
HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
CUTOFF = 0.35          # Tanimoto distance; 0.65 similarity, the usual Butina series threshold
NBOOT = 4000
SEED = 20260921


def auroc(y, s):
    y = np.asarray(y, int); s = np.asarray(s, float)
    pos, neg = s[y == 1], s[y == 0]
    if not len(pos) or not len(neg):
        return float("nan")
    r = stats.rankdata(np.concatenate([pos, neg]))
    return float((r[:len(pos)].sum() - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg)))


def load_arm(exh, seed=1):
    """Merge EVERY complete shard of one seed. The first version kept only the first file,
    which silently reduced the exh=32 arm to one shard of 87 compounds and made the whole
    comparison look far narrower than it is."""
    scores, labels, shards = {}, {}, set()
    for f in sorted(glob.glob(str(HERE / f"results/seed_exh{exh}_s{seed}_*.json"))):
        j = json.loads(Path(f).read_text())
        if not j.get("complete"):
            continue
        scores.update(j["scores"]); labels.update(j["labels"])
        shards.add((j.get("shard", 0), j.get("n_shards", 1)))
    ns = max((n for _, n in shards), default=1)
    got = {sh for sh, n in shards if n == ns}
    if got != set(range(ns)):
        raise SystemExit(f"exh={exh} seed={seed}: incomplete, have shards {sorted(got)} of {ns}")
    return scores, labels


def main() -> int:
    sm = json.loads((REPO / "analysis/boltz2/bundle/mpro_smiles.json").read_text())
    s4, lab = load_arm(4)
    names = sorted(n for n in s4 if n in sm and n in lab)
    print(f"panel: {len(names)} compounds with SMILES and a score")

    fps, keep = [], []
    for n in names:
        m = Chem.MolFromSmiles(sm[n])
        if m is None:
            continue
        fps.append(AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=2048))
        keep.append(n)
    print(f"  {len(keep)} parsed for fingerprinting")

    # Butina: standard scaffold-series clustering.
    dists = []
    for i in range(1, len(fps)):
        dists.extend(1 - np.array(BulkTanimotoSimilarity(fps[i], fps[:i])))
    clusters = Butina.ClusterData(dists, len(fps), CUTOFF, isDistData=True)
    cid = {}
    for ci, members in enumerate(clusters):
        for mi in members:
            cid[keep[mi]] = ci
    sizes = Counter(cid.values())
    big = sorted(sizes.values(), reverse=True)
    y = np.array([lab[n] for n in keep], int)
    print(f"\nseries at Tanimoto {1 - CUTOFF:.2f}: {len(clusters)} clusters over {len(keep)} compounds")
    print(f"  singletons {sum(1 for v in sizes.values() if v == 1)}, "
          f"largest {big[0]}, top-5 sizes {big[:5]}")
    print(f"  effective n (Kish, by series): {sum(big)**2 / sum(v*v for v in big):.1f} "
          f"against {len(keep)} compounds")

    # Compound bootstrap (what the manuscript does) vs cluster bootstrap (resampling series).
    s32, lab32 = load_arm(32)
    common = [n for n in keep if n in s32]
    y = np.array([lab[n] for n in common], int)
    a4 = np.array([s4[n] for n in common], float)
    a32 = np.array([s32[n] for n in common], float)
    obs = auroc(y, a32) - auroc(y, a4)
    groups = np.array([cid[n] for n in common])
    uniq = np.unique(groups)
    idx_by_g = {g: np.where(groups == g)[0] for g in uniq}
    rng = np.random.default_rng(SEED)

    comp, clus = [], []
    for _ in range(NBOOT):
        i = rng.integers(0, len(common), len(common))
        comp.append(auroc(y[i], a32[i]) - auroc(y[i], a4[i]))
        gs = rng.choice(uniq, size=len(uniq), replace=True)
        j = np.concatenate([idx_by_g[g] for g in gs])
        clus.append(auroc(y[j], a32[j]) - auroc(y[j], a4[j]))
    comp = np.array(comp); clus = np.array(clus)
    cl_ = np.nanpercentile(comp, [2.5, 97.5]); ck = np.nanpercentile(clus, [2.5, 97.5])
    print(f"\nexh32 - exh4 on {len(common)} common compounds, observed {obs:+.5f}")
    print(f"  compound bootstrap 95% CI : [{cl_[0]:+.5f}, {cl_[1]:+.5f}]  width {cl_[1]-cl_[0]:.5f}")
    print(f"  CLUSTER  bootstrap 95% CI : [{ck[0]:+.5f}, {ck[1]:+.5f}]  width {ck[1]-ck[0]:.5f}")
    print(f"  widening factor           : {(ck[1]-ck[0]) / (cl_[1]-cl_[0]):.2f}x")
    # The discriminating case. A paired DIFFERENCE sees the same series in both arms, so the
    # clustering largely cancels. An ABSOLUTE AUROC has no such protection, and that is where a
    # Kish effective n of ~137 against 749 compounds should bite. Measuring both separates
    # "the intervals are fine" from "the intervals are fine for the quantity we report".
    ca, ka = [], []
    for _ in range(NBOOT):
        i = rng.integers(0, len(common), len(common))
        ca.append(auroc(y[i], a4[i]))
        gs = rng.choice(uniq, size=len(uniq), replace=True)
        j = np.concatenate([idx_by_g[g] for g in gs])
        ka.append(auroc(y[j], a4[j]))
    ca = np.nanpercentile(ca, [2.5, 97.5]); ka = np.nanpercentile(ka, [2.5, 97.5])
    print(f"\nABSOLUTE AUROC of the exh=4 arm, {auroc(y, a4):.4f}")
    print(f"  compound bootstrap 95% CI : [{ca[0]:.4f}, {ca[1]:.4f}]  width {ca[1]-ca[0]:.4f}")
    print(f"  CLUSTER  bootstrap 95% CI : [{ka[0]:.4f}, {ka[1]:.4f}]  width {ka[1]-ka[0]:.4f}")
    print(f"  widening factor           : {(ka[1]-ka[0]) / (ca[1]-ca[0]):.2f}x")
    print()
    print("READING:")
    print("  Paired DIFFERENCES -- the manuscript's five interventions -- are barely affected,")
    print("  because both arms carry the same series structure and it cancels.")
    print("  ABSOLUTE AUROCs are affected, and any claim quoting one needs the cluster interval.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
