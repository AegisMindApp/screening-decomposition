#!/usr/bin/env python3
"""Build harness input from score files already on disk.

Produces scores/<method>__<target>.json holding {"y","score","D","ids"} row-aligned.
Descriptors are the seven the manuscript uses, computed from SMILES with RDKit.

Nothing here re-runs a calculation; it only reshapes what exists. Targets or methods with no
per-compound data are reported as missing rather than silently skipped -- a method evaluated on
fewer targets than it appears to cover is exactly the error the harness exists to prevent.

    python analysis/method_bench/extract.py
"""
from __future__ import annotations

import glob
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
RB = REPO / "analysis" / "retrospective_benchmark"
B2 = REPO / "analysis" / "boltz2"
OUT = HERE / "scores"


def strip(n):
    return n[len("bench_"):] if n.startswith("bench_") else n


def load_map(*paths, key=None):
    m = {}
    for p in paths:
        for f in glob.glob(str(p)):
            for k, v in json.load(open(f)).items():
                m[strip(k)] = v[key] if key and isinstance(v, dict) else v
    return m


def panel(path):
    """name -> (smiles, label) from a panel file."""
    d = json.load(open(path))
    cs = d["compounds"] if isinstance(d, dict) else d
    return {c["name"]: (c["smiles"], int(c["label"])) for c in cs}


def mpro_panel():
    """Mpro labels come from the seed-floor files; SMILES from the pdbqt REMARK line."""
    lab = json.load(open(REPO / "analysis/seed_floor/results/seed_exh4_s1_shard0of1.json"))["labels"]
    smi = {}
    for f in glob.glob(str(RB / "kaggle/bundle_CHEMBL4523582_7VU6*/shard*/ligands/*.pdbqt")):
        name = strip(Path(f).stem)
        if name in smi:
            continue
        for line in open(f):
            if line.startswith("REMARK SMILES") and "IDX" not in line:
                smi[name] = line.split("REMARK SMILES", 1)[1].strip()
                break
    return {n: (smi[n], int(lab[n])) for n in lab if n in smi}


def descriptors(smiles):
    from rdkit import Chem, RDLogger
    from rdkit.Chem import Crippen, Descriptors, rdMolDescriptors
    RDLogger.DisableLog("rdApp.*")
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return None
    return [Descriptors.MolWt(m), Crippen.MolLogP(m), rdMolDescriptors.CalcTPSA(m),
            rdMolDescriptors.CalcNumHBD(m), rdMolDescriptors.CalcNumHBA(m),
            rdMolDescriptors.CalcNumRotatableBonds(m),
            rdMolDescriptors.CalcNumAromaticRings(m)]


def build(method, target, scores, pan):
    ids, y, s, D = [], [], [], []
    for n, sc in scores.items():
        if n not in pan:
            continue
        d = descriptors(pan[n][0])
        if d is None:
            continue
        ids.append(n); y.append(pan[n][1]); s.append(float(sc)); D.append(d)
    if len(ids) < 50:
        print(f"  {method:10} {target:6} SKIP -- only {len(ids)} compounds with score+SMILES")
        return None
    OUT.mkdir(exist_ok=True)
    p = OUT / f"{method}__{target}.json"
    p.write_text(json.dumps({"ids": ids, "y": y, "score": s, "D": D}))
    print(f"  {method:10} {target:6} {len(ids):5d} compounds, {sum(y):4d} active -> {p.name}")
    return p


def main() -> int:
    pans = {"Mpro": mpro_panel(),
            "FXa": panel(RB / "fxa_panel.json"),
            "PDL1": panel(RB / "pdl1_panel.json")}
    for t, p in pans.items():
        print(f"panel {t}: {len(p)} compounds with SMILES+label")

    print("\nVina (more negative = better, so negate for 'higher is better'):")
    vina = {
        "Mpro": {k: -v for k, v in load_map(RB / "vina_cache_CHEMBL4523582_7VU6.json").items()},
        "FXa": {k: -v for k, v in load_map(RB / "kaggle/vina_cache_CHEMBL244_shard*.json").items()},
        "PDL1": {k: -v for k, v in load_map(RB / "pdl1_vina_merged.json").items()},
    }
    for t, sc in vina.items():
        build("vina", t, sc, pans[t])

    print("\nBoltz-2 (affinity_probability_binary, higher = better):")
    b2 = {
        "Mpro": load_map(B2 / "results/boltz2_scores.json", key="affinity_probability_binary"),
        "FXa": load_map(B2 / "results_fxa/boltz2_scores_s*.json",
                        key="affinity_probability_binary"),
    }
    for t, sc in b2.items():
        build("boltz2", t, sc, pans[t])
    print("  boltz2     PDL1   MISSING -- never run; this is what blocks a transferability "
          "verdict for Boltz-2")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
