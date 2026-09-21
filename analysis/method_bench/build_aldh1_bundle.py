#!/usr/bin/env python3
"""Extract the ALDH1 sequence and binding pocket for Boltz-2, from the LIT-PCBA structure.

Boltz-2 co-folds from SEQUENCE plus ligand SMILES, so no receptor preparation is needed -- but
the pocket constraint used on Mpro arm 2 and Factor Xa is reproduced here so the third target
is run the same way as the other two. A method compared across targets under different
conditions is not being compared.

    python analysis/method_bench/build_aldh1_bundle.py
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
SRC = HERE.parents[1] / "analysis" / "litpcba" / "full" / "ALDH1"
PDB = "4wp7"                     # first structure in the LIT-PCBA set for this target
POCKET_CUTOFF = 8.0              # A, matching Mpro arm 2 and Factor Xa

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
    "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
    "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M", "HID": "H", "HIE": "H", "HIP": "H", "CYX": "C", "ASH": "D", "GLH": "E",
}


def atoms(path):
    txt = Path(path).read_text()
    m = re.search(r"@<TRIPOS>ATOM(.*?)(@<TRIPOS>|\Z)", txt, re.S)
    out = []
    for line in m.group(1).strip().splitlines():
        f = line.split()
        if len(f) >= 8:
            out.append({"xyz": np.array([float(f[2]), float(f[3]), float(f[4])]),
                        "res": f[7]})
    return out


def split_res(tag):
    m = re.match(r"([A-Za-z]+)(\d+)", tag)
    return (m.group(1).upper(), int(m.group(2))) if m else (None, None)


def main() -> int:
    prot = atoms(SRC / f"{PDB}_protein.mol2")
    lig = atoms(SRC / f"{PDB}_ligand.mol2")
    print(f"{PDB}: {len(prot)} protein atoms, {len(lig)} ligand atoms")

    # sequence, in residue-number order, skipping anything not a standard amino acid
    seen, seq = {}, []
    for a in prot:
        name, num = split_res(a["res"])
        if num is None or num in seen:
            continue
        one = THREE_TO_ONE.get(name)
        if one:
            seen[num] = one
    for num in sorted(seen):
        seq.append(seen[num])
    seq = "".join(seq)
    unknown = {split_res(a["res"])[0] for a in prot} - set(THREE_TO_ONE) - {None}
    print(f"sequence: {len(seq)} residues, span {min(seen)}-{max(seen)}")
    if unknown:
        print(f"  non-standard residues skipped: {sorted(unknown)}")

    # Pocket from the ligand CENTROID, not from every ligand atom.
    #
    # Mpro arm 2 and Factor Xa both derived their pocket from a single point -- the Vina box
    # centre -- at 8 A, giving 13 residues each. Measuring from all 63 ligand atoms instead
    # sweeps up 66 residues, five times as many, which is not a pocket constraint and not the
    # condition the other two targets ran under. A method compared across targets under
    # different conditions is not being compared.
    L = np.array([a["xyz"] for a in lig])
    centre = L.mean(axis=0)
    pocket = set()
    for a in prot:
        name, num = split_res(a["res"])
        if num in seen and np.linalg.norm(centre - a["xyz"]) <= POCKET_CUTOFF:
            pocket.add(num)
    pocket = sorted(pocket)
    print(f"pocket: {len(pocket)} residues within {POCKET_CUTOFF} A of the ligand centroid "
          f"(Mpro and FXa both gave 13 under the same rule)")
    print(f"  {[f'{seen[n]}{n}' for n in pocket[:10]]}...")

    if len(seq) < 300 or not pocket:
        print("REFUSING to write: sequence or pocket extraction failed")
        return 1

    out = HERE / "panels" / "ALDH1_bundle.json"
    out.write_text(json.dumps({
        "target": "ALDH1", "pdb": PDB, "chain": "A",
        "sequence": seq, "n_residues": len(seq),
        "residue_span": [int(min(seen)), int(max(seen))],
        "pocket_residues": pocket, "pocket_cutoff_A": POCKET_CUTOFF,
        "pocket_origin": "ligand centroid (matches the Vina box centre used for Mpro/FXa)",
        "box_centre": [round(float(x), 3) for x in centre],
        "pocket_labels": [f"{seen[n]}{n}" for n in pocket],
        "note": ("Pocket conditioning reproduces the Mpro arm-2 and Factor Xa configuration so "
                 "the third target is run identically to the other two."),
    }, indent=1))
    print(f"-> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
