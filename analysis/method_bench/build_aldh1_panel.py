#!/usr/bin/env python3
"""Build the ALDH1 panel — the clean third target. Selection recorded in THIRD_TARGET.md.

Re-checks admissibility ON THE BUILT PANEL and refuses to write if it fails. The published
LIT-PCBA baseline is computed on 25,000 imbalanced inactives; balance changes AUROC, so the
number that matters is the one for the panel we actually use.

    python analysis/method_bench/build_aldh1_panel.py
"""
from __future__ import annotations

import json
import random
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from harness import MAX_DESCRIPTOR_BASELINE, target_admissible  # noqa: E402

SRC = HERE.parents[1] / "analysis" / "litpcba" / "full" / "ALDH1"
SEED = 20260919          # fixed before selection
N_ACT, N_INACT = 300, 450


def feats(smi):
    from rdkit import Chem, RDLogger
    from rdkit.Chem import Crippen, Descriptors, rdMolDescriptors
    RDLogger.DisableLog("rdApp.*")
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return None
    return [Descriptors.MolWt(m), Crippen.MolLogP(m), rdMolDescriptors.CalcTPSA(m),
            rdMolDescriptors.CalcNumHBD(m), rdMolDescriptors.CalcNumHBA(m),
            rdMolDescriptors.CalcNumRotatableBonds(m),
            rdMolDescriptors.CalcNumAromaticRings(m)]


def take(path, n):
    rows = [l.split()[0] for l in open(path) if l.strip()]
    random.Random(SEED).shuffle(rows)
    out, skipped = [], 0
    for s in rows:
        f = feats(s)
        if f is None:
            skipped += 1
            continue
        out.append((s, f))
        if len(out) >= n:
            break
    return out, skipped


def main() -> int:
    A, skip_a = take(SRC / "actives.smi", N_ACT)
    I, skip_i = take(SRC / "inactives.smi", N_INACT)
    print(f"actives   {len(A)} (skipped {skip_a} unparseable)")
    print(f"inactives {len(I)} (skipped {skip_i} unparseable)")

    smis = [s for s, _ in A + I]
    D = np.array([f for _, f in A + I])
    y = np.array([1] * len(A) + [0] * len(I))
    ok, base = target_admissible(y, D)
    print(f"descriptor baseline on THIS panel: {base:.4f} "
          f"(limit {MAX_DESCRIPTOR_BASELINE}) -> {'ADMISSIBLE' if ok else 'INADMISSIBLE'}")
    if not ok:
        print("REFUSING to write: a panel properties already separate cannot test whether a "
              "method exceeds properties. Amend THIRD_TARGET.md and choose another target.")
        return 1

    out = HERE / "panels"
    out.mkdir(exist_ok=True)
    p = out / "ALDH1_panel.json"
    p.write_text(json.dumps({
        "target": "ALDH1", "source": "LIT-PCBA full", "seed": SEED,
        "n": len(y), "n_active": int(y.sum()), "pct_active": round(float(y.mean()), 4),
        "descriptor_baseline": base, "admissible": ok,
        "unparseable_skipped": {"actives": skip_a, "inactives": skip_i},
        "ids": [f"ALDH1_{i:04d}" for i in range(len(y))],
        "smiles": smis, "y": y.tolist(), "D": D.tolist(),
    }, indent=1))
    print(f"-> {p}  ({len(y)} compounds, {int(y.sum())} active)")
    print(f"Boltz-2 cost at 101 s/ligand: {len(y)*101/3600:.1f} GPU hours "
          f"-> {int(np.ceil(len(y)*101/3600/7))} shards of ~7h")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
