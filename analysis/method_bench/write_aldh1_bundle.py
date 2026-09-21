#!/usr/bin/env python3
"""Write the ALDH1 Boltz-2 bundle in the exact layout boltz2_worker.py reads.

The layout is copied from bundle_fxa so the third target runs through the identical code path:
sequence.txt, <target>_smiles.json, manifest.json (compounds with name+label), pocket.json.

The one trap, stated in the FXa pocket.json itself: pocket entries are **1-BASED SEQUENCE
INDICES into sequence.txt, not PDB residue numbers** -- Boltz-2 receives a sequence and never
sees PDB numbering. ALDH1's structure starts at residue 8 and may have gaps, so the mapping is
by rank in the ordered residue list, not a fixed offset.

    python analysis/method_bench/write_aldh1_bundle.py
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
SRC = HERE.parents[1] / "analysis" / "litpcba" / "full" / "ALDH1"
OUT = HERE / "bundle_aldh1"
PDB, CUT = "4wp7", 8.0
T31 = {"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E","GLY":"G",
       "HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P","SER":"S",
       "THR":"T","TRP":"W","TYR":"Y","VAL":"V","MSE":"M","HID":"H","HIE":"H","HIP":"H",
       "CYX":"C","ASH":"D","GLH":"E"}


def atoms(p):
    m = re.search(r"@<TRIPOS>ATOM(.*?)(@<TRIPOS>|\Z)", Path(p).read_text(), re.S)
    out = []
    for line in m.group(1).strip().splitlines():
        f = line.split()
        if len(f) >= 8:
            out.append((np.array([float(f[2]), float(f[3]), float(f[4])]), f[7]))
    return out


def split_res(tag):
    m = re.match(r"([A-Za-z]+)(\d+)", tag)
    return (m.group(1).upper(), int(m.group(2))) if m else (None, None)


def main() -> int:
    prot, lig = atoms(SRC / f"{PDB}_protein.mol2"), atoms(SRC / f"{PDB}_ligand.mol2")
    seen = {}
    for xyz, tag in prot:
        name, num = split_res(tag)
        one = T31.get(name)
        if one and num not in seen:
            seen[num] = one
    order = sorted(seen)                       # PDB numbers, ascending
    seq = "".join(seen[n] for n in order)
    # PDB number -> 1-BASED index into seq. By rank, not offset: gaps would break an offset.
    pdb_to_idx = {n: i + 1 for i, n in enumerate(order)}

    centre = np.mean([x for x, _ in lig], axis=0)
    near = {split_res(t)[1] for x, t in prot if np.linalg.norm(centre - x) <= CUT}
    pocket_pdb = sorted(n for n in near if n in pdb_to_idx)
    pocket = [["A", pdb_to_idx[n]] for n in pocket_pdb]
    annotated = [["A", pdb_to_idx[n], f"{seen[n]}{n}",
                  round(float(min(np.linalg.norm(centre - x) for x, t in prot
                                  if split_res(t)[1] == n)), 2)] for n in pocket_pdb]

    panel = json.loads((HERE / "panels" / "ALDH1_panel.json").read_text())
    ids, smis, ys = panel["ids"], panel["smiles"], panel["y"]

    OUT.mkdir(exist_ok=True)
    (OUT / "sequence.txt").write_text(seq + "\n")
    (OUT / "aldh1_smiles.json").write_text(json.dumps(dict(zip(ids, smis)), indent=1))
    (OUT / "manifest.json").write_text(json.dumps(
        {"compounds": [{"name": i, "label": int(y)} for i, y in zip(ids, ys)]}, indent=1))
    (OUT / "pocket.json").write_text(json.dumps({
        "chains": ["A"], "pocket": pocket, "pocket_annotated": annotated,
        "note": ("1-BASED SEQUENCE INDICES into sequence.txt, NOT PDB residue numbers. "
                 "Boltz-2 receives a sequence and never sees PDB numbering. ALDH1's structure "
                 "starts at residue 8, so indices are by RANK in the ordered residue list."),
        "box_centre": [round(float(x), 2) for x in centre]}, indent=1))

    print(f"sequence {len(seq)} aa (PDB {order[0]}-{order[-1]})")
    print(f"pocket   {len(pocket)} residues, e.g. {annotated[:3]}")
    print(f"panel    {len(ids)} compounds, {sum(ys)} active")
    # The trap, checked rather than trusted: index 1 must be the first residue of the sequence.
    first = order[0]
    assert pdb_to_idx[first] == 1 and seq[0] == seen[first], "sequence index mapping is wrong"
    mx = max(i for _, i in pocket)
    assert mx <= len(seq), f"pocket index {mx} exceeds sequence length {len(seq)}"
    print(f"index mapping verified: PDB {first} -> index 1; max pocket index {mx} <= {len(seq)}")
    print(f"-> {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
