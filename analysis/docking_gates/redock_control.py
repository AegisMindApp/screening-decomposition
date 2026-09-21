"""Redocking positive control — the pre-flight gate that CTSS would have failed.

WHY THIS EXISTS

The CTSS/3N4C screen produced hypotheses (sertraline, aripiprazole analogs) from docking scores
that turned out to rest on nothing. The box was correct, centred on the co-crystal ligand's own
centroid to 0.00 A. The ligand was correct. And Vina still placed EF3 -- 3N4C's own crystallised
inhibitor -- 7.2 A from where the crystal puts it, then scored that wrong pose better than 61 of
72 decoys.

None of the existing pre-flight gates catch this. Target-identity passes (right protein, right
chain). Box-sanity passes (box is on the ligand). Counter-screen passes (we ran one). The assay
is broken at a level none of them look at: whether the protocol can reproduce a known answer.

THE SUBTLETY THAT MATTERS

A naive check -- "is a good pose somewhere in the output?" -- would have PASSED CTSS. At
exhaustiveness 128 Vina generates near-native poses twice (2.22 A and 2.34 A) and ranks them
6th and 9th of 9. Sampling finds the binding mode; the scoring function rejects it. Since a
prospective screen only ever sees the top pose, a protocol that buries the right answer at rank 6
is worthless even though it "found" it.

So this gate fails a target on EITHER condition:
  1. no pose within `pass_rmsd` of the crystal pose        -> sampling cannot reach it
  2. such a pose exists but is not rank 1                  -> scoring rejects it  [the CTSS case]

RMSD is computed with RDKit CalcRMS, which enumerates graph automorphisms so a symmetric ring
flip is not counted as displacement, and which does NOT superimpose -- redocking asks whether the
pose is in the right place, so aligning the two molecules first would answer the wrong question. An earlier hand-rolled element-matching RMSD was used during
the CTSS investigation; it is a LOWER BOUND (it may pair chemically distinct carbons) and is not
good enough to gate on.

The ligand is rebuilt from SMILES and prepared through the same RDKit/Meeko path production uses,
so the gate tests the pipeline's own prep chain rather than a privileged one. If prep is what is
broken, this gate catches that too.

CLI:
    python redock_control.py --receptor path/to/receptor.pdbqt --pdb 3N4C \
        --center 14.004 9.016 21.712 --size 24 24 24 --exhaustiveness 16
"""
from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
import tempfile
import urllib.request
from pathlib import Path

import numpy as np

# Residues that are crystallisation additives, not the ligand of interest. Without this the
# "largest HETATM" heuristic happily picks a sulfate.
_JUNK = {
    "HOH", "WAT", "DOD", "SO4", "PO4", "GOL", "EDO", "PEG", "DMS", "ACT", "MPD", "TRS",
    "IMD", "FMT", "CL", "NA", "MG", "CA", "ZN", "K", "MN", "FE", "NI", "CD", "CU", "IOD",
    "BR", "NO3", "ACY", "BME", "CIT", "TLA", "SIN", "PG4", "1PE", "P6G", "MES", "EPE",
    # Detergents, lipids and PEG oligomers. These are crystallisation additives, but they are
    # often LARGER than the drug, so "largest non-junk HETATM" hands the gate a surfactant: 3L4D's
    # N8E (C18H38O6, 24 heavy atoms) outranks fluconazole (22) and would have been redocked as if
    # it were the co-crystal ligand. A blocklist is the only workable shape for junk -- there is no
    # finite allowlist of "real ligands" -- and it fails SAFE: a wrong entry here makes the gate
    # refuse to find a ligand and issue no verdict, never makes it issue a passing one.
    "N8E", "C8E", "LDA", "LMT", "DMU", "BOG", "BNG", "OCT", "HEX", "MYR", "PLM", "OLA", "OLB",
    "OLC", "SDS", "TRT", "PEE", "PCW", "PE4", "PE5", "PE8", "PGE", "2PE", "7PE", "12P", "15P",
    "XPE", "33O", "P33", "UNL", "UNX", "MRD", "BU3", "SUC", "TAR", "MLI", "FLC", "NHE",
}

# Cofactors: part of the BINDING SITE, not the thing being redocked. They must survive receptor
# prep (an apo pocket is a different pocket -- deleting 3L4D's heme moved the fluconazole redock
# to 9.02 A) and must never be selected as the ligand (the "largest non-junk HETATM" rule picked
# HEM's 43 atoms over fluconazole's 22, and RDKit then choked on the Fe-coordinated nitrogens).
#
# An ALLOWLIST, deliberately. A blocklist of "things that aren't ligands" cannot be verified; this
# set can be read against the validation targets. Names follow screen_cyp51.py's vocabulary rather
# than inventing a second one that can drift from it.
#
# Nucleotides (ADP/ATP/GDP/NAD...) are NOT here even though they are cofactors in context: in a
# kinase or a dehydrogenase they are frequently the ligand you want to redock, and mis-classifying
# the ligand as scenery is the failure this set exists to prevent. Same reason BTN is absent --
# biotin is 1STP's ligand, and listing it would silently delete a validation target.
_COFACTOR = {
    "HEM", "HEC", "HEA", "HEB", "HAS", "SRM", "VER",   # haem variants
    "FAD", "FMN", "FES", "SF4", "F3S",                 # flavins, Fe-S clusters
    "PLP", "SAM", "SAH", "COA", "MQ7", "UQ1",          # PLP, methyl donors, CoA, quinones
}
# Standalone ions. Already in _JUNK (so already excluded from ligand choice); named separately
# because prepare_receptor has to KEEP them, which is the opposite question.
_METAL_ION = {"MG", "ZN", "CA", "MN", "FE", "FE2", "NI", "CU", "CU1", "CO", "CD", "K", "NA"}

_SCORE_RE = re.compile(r"^\s+(\d+)\s+(-?\d+\.\d+)", re.M)


# --------------------------------------------------------------------------- structure input

def fetch_pdb(pdb_id: str, workdir: Path) -> Path:
    dest = workdir / f"{pdb_id.upper()}.pdb"
    if not dest.exists():
        url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
        with urllib.request.urlopen(url, timeout=60) as r:
            dest.write_bytes(r.read())
    return dest


def fetch_smiles(resname: str) -> str | None:
    """Canonical SMILES from the RCSB chemical component dictionary.

    Taken from the CCD rather than inferred from coordinates: bond orders and protonation guessed
    from a PDB block are exactly the kind of silent error this gate exists to catch.
    """
    url = f"https://data.rcsb.org/rest/v1/core/chemcomp/{resname.upper()}"
    try:
        with urllib.request.urlopen(url, timeout=60) as r:
            d = json.load(r)
    except Exception:                                                    # noqa: BLE001
        return None
    for entry in d.get("rcsb_chem_comp_descriptor", {}), d.get("pdbx_chem_comp_descriptor", []):
        if isinstance(entry, dict):
            for k in ("smiles_stereo", "smiles"):
                if entry.get(k):
                    return entry[k]
        elif isinstance(entry, list):
            best = None
            for e in entry:
                if e.get("type", "").upper() == "SMILES_CANONICAL":
                    if e.get("program", "").upper().startswith("OPENEYE"):
                        return e["descriptor"]
                    best = best or e.get("descriptor")
            if best:
                return best
    return None


def covalent_to_protein(pdb_path: Path, lig_atoms, cutoff=1.95):
    """Is the co-crystal ligand covalently bonded to the protein?

    Checks ligand heavy atoms against protein SG/OG/OG1/NE2/SD -- the nucleophiles that form
    covalent adducts. This matters because a NON-covalent redocking control cannot reproduce a
    covalent pose, so a bare "sampling failed" verdict would misdiagnose the protocol. Cathepsins
    are cysteine proteases and most of their inhibitors are covalent, which is how this was found.
    """
    import numpy as _np
    # ANY protein heavy atom, not just named nucleophiles. Restricting to SG/OG/... missed 1HNE,
    # where the selected HETATM is a capping fragment of a peptidyl inhibitor bonded to an
    # adjacent modified residue rather than to a catalytic side chain -- 1.31 A, plainly covalent,
    # and silently redocked as if it were a free 8-atom ligand.
    nuc = []
    for l in pdb_path.read_text().splitlines():
        if l.startswith("ATOM"):
            el = (l[76:78].strip() or l[12:16].strip()[0]).upper()
            if el == "H":
                continue
            nuc.append((l[17:20].strip() + l[22:27].strip() + l[21],
                        _np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])))
    best = (1e9, None)
    for _, p in lig_atoms:
        for tag, q in nuc:
            d = float(_np.linalg.norm(p - q))
            if d < best[0]:
                best = (d, tag)
    return (best[0] <= cutoff), round(best[0], 2), best[1]


def metal_coordination(pdb_path: Path, lig_atoms, cutoff=2.6):
    """Does the co-crystal ligand coordinate a cofactor metal? Reported, NOT gating.

    `covalent_to_protein` scans ATOM records only, so a HETATM metal is invisible to it: in 3L4D
    fluconazole's N5 sits 1.88 A from the haem Fe -- inside that function's own 1.95 A cutoff --
    yet the nearest thing it can see is PHE289A at 3.08 A, and it reports "not covalent".

    Deliberately NOT wired to a NOT_APPLICABLE verdict. Coordination is not a covalent adduct, and
    Vina models Fe-N as vdW rather than refusing it, so whether the pose is reproducible is an
    empirical question. Auto-declining every metalloenzyme would be a disqualification dressed up
    as a measurement -- the same shape as the finding this gate exists to prevent. It is recorded
    so that a FAIL on a metalloenzyme can be read for what it is instead of blamed on sampling.
    """
    import numpy as _np
    metals = []
    for l in pdb_path.read_text().splitlines():
        if not l.startswith("HETATM"):
            continue
        rn = l[17:20].strip().upper()
        if rn not in _COFACTOR and rn not in _METAL_ION:
            continue
        el = (l[76:78].strip().upper() if len(l) >= 78 else "")
        if el not in AD4_TYPES or el in ("C", "N", "O", "S", "H", "F", "P"):
            continue  # AD4_TYPES' remaining keys are exactly the metals
        metals.append((f"{rn}{l[22:27].strip()}{l[21]}.{l[12:16].strip()}",
                       _np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])))
    if not metals:
        return False, None, None
    best = (1e9, None)
    for _, p in lig_atoms:
        for tag, q in metals:
            d = float(_np.linalg.norm(p - q))
            if d < best[0]:
                best = (d, tag)
    return (best[0] <= cutoff), round(best[0], 2), best[1]


def find_ligand(pdb_path: Path, resname: str | None = None, chain: str | None = None):
    """Pick the co-crystal ligand: the largest non-junk HETATM residue, or the one named.

    Returns (resname, chain, [(element, xyz), ...], raw PDB lines).
    A structure with several copies (3N4C is a dimer) yields ONE chain's copy -- merging them
    produces a meaningless reference, which cost an hour during the CTSS investigation.
    """
    groups: dict[tuple, list] = {}
    for line in pdb_path.read_text().splitlines():
        if not line.startswith("HETATM"):
            continue
        rn, ch = line[17:20].strip(), line[21]
        # A cofactor is scenery, never the probe. Without this the heuristic below picks HEM (43
        # heavy atoms) over fluconazole (22) in 3L4D and the gate returns ERROR on every P450.
        # An explicit `resname=` still wins -- asking for HEM by name is a deliberate act.
        if rn in _JUNK or (rn in _COFACTOR and not resname):
            continue
        if resname and rn.upper() != resname.upper():
            continue
        if chain and ch != chain:
            continue
        groups.setdefault((rn, ch, line[22:27].strip()), []).append(line)

    # Drop modified amino acids. A residue carrying a full backbone (N, CA, C, O) is part of the
    # polymer, not a bound ligand -- 1TRN's PTR is phosphotyrosine and 1HNE's ALV is an alanine
    # variant inside a peptidyl inhibitor. Redocking one of these measures nothing: the "ligand"
    # is covalently in the chain, so no free-ligand pose exists to reproduce. Left unfiltered,
    # this silently manufactured 8-12 A "failures" on perfectly ordinary receptors.
    #
    # Names alone cannot decide this. A FREE ligand carrying an amino-acid fragment has the same
    # four names: methotrexate's glutamate tail is literally `C O N CA`, so MTX -- the co-crystal
    # ligand of BOTH 1E7W (PTR1) and 1U72 (DHFR) -- was silently discarded and the gate returned
    # "no ligand" on a target it should have ruled on. Glutathione and GCG escape only because
    # their equivalent atoms happen to be numbered (CA1/C1). Nor does any name-level rescue work:
    # 1TRN's PTR carries P/O1P/O2P/O3P and 1HNE's ALV is a bare 5-atom alanine, so neither
    # "non-standard atom names" nor "terminal carboxylate" separates free ligand from polymer.
    #
    # So the heuristic is kept for auto-selection and skipped when the caller NAMES the ligand --
    # the same rule the cofactor filter above already follows. Naming MTX is a deliberate act;
    # guessing at it is not. This cannot loosen a verdict: an explicitly named residue that really
    # is in the chain still meets covalent_to_protein below and returns NOT_APPLICABLE, which is
    # how 1HNE's peptidyl inhibitor is caught. Verified inert for the validation suite -- all six
    # manifest entries auto-select (validate_gate.py calls find_ligand(p) with no resname), and
    # both polymer cases were already resolved by selection moving on, 1HNE to MSU and 1TRN to ISP.
    BACKBONE = {"N", "CA", "C", "O"}
    modified = []
    if not resname:
        for k, ls in list(groups.items()):
            names = {l[12:16].strip().upper() for l in ls}
            if BACKBONE <= names:
                modified.append(k[0])
                del groups[k]
    if not groups:
        return None
    # Largest by heavy-atom count; ties break toward the earliest chain for determinism.
    def heavy(ls):
        return sum(1 for l in ls if ((l[76:78].strip() or l[12:16].strip()[0]).upper() != "H"))
    key = max(groups, key=lambda k: (heavy(groups[k]), -ord(k[1])))
    lines = groups[key]
    atoms = []
    for l in lines:
        el = (l[76:78].strip() or l[12:16].strip()[0]).upper()
        if el == "H":
            continue
        atoms.append((el, np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])))
    return key[0], key[1], atoms, lines


_ENCLOSE_MIN_FRAC = 0.5


def crystal_enclosure(pdb_path, lig_atoms, within=5.0):
    """The same enclosure fraction, measured against the RAW crystal: every ATOM record, all chains,
    no preparation of any kind.

    This exists to stop RECEPTOR_SITE_MISMATCH asserting a cause it has not checked. The verdict's
    message used to say the receptor was "probably missing the chain the ligand binds to" -- true for
    3THW, where ADP is in chain A and the receptor was chain B only, but simply wrong for 4ADW, whose
    GCG sits in the wide trypanothione site with 65% of its atoms in solvent. Measured on the
    untouched crystal with both chains present, 4ADW's enclosure is 35% -- so there is no chain to
    add and no prep defect to fix, and telling someone to "prepare the correct chain" sends them
    after a bug that does not exist.

    A low fraction here means the site is genuinely open, which is a fact about the target rather
    than about our preparation of it -- and a shallow, solvent-exposed site is exactly where a
    docking score is least interpretable, since buried contact is what the function scores.
    """
    import numpy as _np
    P = [[float(l[30:38]), float(l[38:46]), float(l[46:54])]
         for l in Path(pdb_path).read_text().splitlines() if l.startswith("ATOM")]
    if not P:
        return None
    L = _np.asarray([a[1] for a in lig_atoms])
    d = _np.linalg.norm(_np.asarray(P)[None, :, :] - L[:, None, :], axis=2)
    return round(float((d.min(axis=1) < within).mean()), 2)


def receptor_encloses(receptor_pdbqt, lig_atoms, within=5.0, min_frac=_ENCLOSE_MIN_FRAC):
    """Does the PREPARED receptor actually contain this ligand's binding site?

    3THW's ADP sits in chain A while the MSH3 receptor is chain B only, so redocking it measured
    nothing but the distance between two chains -- and reported a 12 A "sampling failure" that
    said nothing about the protocol. That is a receptor/ligand mismatch and deserves its own
    verdict, since the fix is to prepare the right chain, not to abandon the target.
    """
    import numpy as _np
    R = []
    for l in Path(receptor_pdbqt).read_text().splitlines():
        if l.startswith(("ATOM", "HETATM")):
            R.append([float(l[30:38]), float(l[38:46]), float(l[46:54])])
    if not R:
        return False, 0.0, None
    R = _np.asarray(R)
    L = _np.asarray([a[1] for a in lig_atoms])
    d = _np.linalg.norm(R[None, :, :] - L[:, None, :], axis=2)
    frac = float((d.min(axis=1) < within).mean())
    return frac >= min_frac, round(frac, 2), round(float(d.min()), 2)


# --------------------------------------------------------------------------- ligand prep

def ligand_pdbqt_from_smiles(smiles: str, out_path: Path) -> bool:
    """Same RDKit embed + MMFF + Meeko path production uses (run_acrb_efflux_screen.py)."""
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False
    mol = Chem.AddHs(mol)
    if AllChem.EmbedMolecule(mol, randomSeed=0xF00D) != 0:
        return False
    try:
        AllChem.MMFFOptimizeMolecule(mol, maxIters=2000)
    except Exception:                                                    # noqa: BLE001
        pass
    from meeko import MoleculePreparation, PDBQTWriterLegacy
    setups = MoleculePreparation().prepare(mol)
    if not setups:
        return False
    res = PDBQTWriterLegacy.write_string(setups[0])
    out_path.write_text(res[0] if isinstance(res, tuple) else res)
    return True


# --------------------------------------------------------------------------- RMSD

def _pdbqt_poses_as_pdb(path: Path) -> list[str]:
    """Split a Vina output PDBQT into per-pose PDB blocks (heavy atoms, AutoDock types mapped)."""
    out, cur = [], None
    for l in path.read_text().splitlines():
        if l.startswith("MODEL"):
            cur = []
        elif l.startswith("ENDMDL"):
            if cur is not None:
                out.append("\n".join(cur) + "\nEND\n")
            cur = None
        elif cur is not None and l.startswith(("ATOM", "HETATM")):
            el = l[77:79].strip().upper()
            if el.startswith("H"):
                continue
            el = {"A": "C", "NA": "N", "OA": "O", "SA": "S", "NS": "N", "OS": "O"}.get(el, el)
            cur.append(l[:76] + f"{el:>2s}" + " " * 2)
    return out


def best_rms(template_smiles: str, probe_pdb_block: str, ref_pdb_block: str) -> float | None:
    """Symmetry-corrected, non-superimposed RMSD via RDKit CalcRMS.

    Both poses get their bond orders from the same CCD template, so the graphs match and the
    automorphism search is meaningful.
    """
    from rdkit import Chem, RDLogger
    from rdkit.Chem import AllChem, rdMolAlign
    RDLogger.DisableLog("rdApp.warning")     # symmetric matches are expected here, not a problem

    tpl = Chem.MolFromSmiles(template_smiles)
    if tpl is None:
        return None
    # RemoveAllHs, not RemoveHs: CCD SMILES carry explicit stereo-defining hydrogens (EF3's imine
    # is written [H]/N=C/...), which RemoveHs keeps. The template then has 34 atoms against the
    # PDB block's 33 and every substructure match fails with "No matching found" -- which surfaces
    # as "RMSD could not be computed", i.e. a silent gate bypass rather than a loud error.
    tpl = Chem.RemoveAllHs(tpl)
    try:
        probe = Chem.MolFromPDBBlock(probe_pdb_block, removeHs=True, sanitize=False)
        ref = Chem.MolFromPDBBlock(ref_pdb_block, removeHs=True, sanitize=False)
        if probe is None or ref is None:
            return None
        probe = AllChem.AssignBondOrdersFromTemplate(tpl, probe)
        ref = AllChem.AssignBondOrdersFromTemplate(tpl, ref)
        # GetBestRMS superimposes; we want in-place displacement, so map without alignment.
        return rdMolAlign.CalcRMS(probe, ref)
    except Exception:                                                    # noqa: BLE001
        return None


def _fallback_rms(probe_pdb_block: str, ref_pdb_block: str) -> float | None:
    """Optimal same-element assignment. A LOWER BOUND on true RMSD -- it may pair chemically
    distinct carbons -- so it is only trustworthy when it already exceeds the pass threshold.
    Used when the CCD template will not match (charged groups such as nitro are a common cause)."""
    from scipy.optimize import linear_sum_assignment

    def parse(b):
        A = []
        for l in b.splitlines():
            if l.startswith(("ATOM", "HETATM")):
                el = (l[76:78].strip() or l[12:16].strip()[0]).upper()
                if el == "H":
                    continue
                A.append((el, np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])))
        return A
    X, Y = parse(ref_pdb_block), parse(probe_pdb_block)
    if not X or not Y:
        return None
    C = np.full((len(X), len(Y)), 1e9)
    for i, (e, p) in enumerate(X):
        for j, (f, q) in enumerate(Y):
            if e == f:
                C[i, j] = float(np.sum((p - q) ** 2))
    r, c = linear_sum_assignment(C)
    d = [C[i, j] for i, j in zip(r, c) if C[i, j] < 1e8]
    return float(np.sqrt(np.mean(d))) if d else None


# --------------------------------------------------------------------------- the gate

def redock_control(receptor_pdbqt, pdb_id, box_center, box_size=(24, 24, 24), *,
                   ligand_resname=None, chain=None, exhaustiveness=16, num_modes=9,
                   vina_bin="vina", pass_rmsd=2.0, box_tolerance=6.0, workdir=None, seed=42):
    """Dock a structure's own co-crystal ligand back into the production box.

    Returns a verdict dict. PASS requires a pose within `pass_rmsd` of the crystal position AND
    that pose ranking first -- see the module docstring for why the second clause is not optional.
    """
    tmp = Path(workdir or tempfile.mkdtemp(prefix="redock_"))
    tmp.mkdir(parents=True, exist_ok=True)
    out: dict = {"pdb_id": pdb_id.upper(), "receptor": str(receptor_pdbqt),
                 "box_center": list(map(float, box_center)), "box_size": list(box_size),
                 "exhaustiveness": exhaustiveness, "pass_rmsd": pass_rmsd,
                 "verdict": "ERROR", "reasons": []}

    pdb = fetch_pdb(pdb_id, tmp)
    found = find_ligand(pdb, ligand_resname, chain)
    if not found:
        out["reasons"].append("no co-crystal ligand found — cannot run a redocking control")
        return out
    resn, ch, atoms, raw = found
    cen = np.mean([a[1] for a in atoms], axis=0)
    offset = float(np.linalg.norm(cen - np.asarray(box_center, dtype=float)))
    out["ligand"] = {"resname": resn, "chain": ch, "n_heavy": len(atoms),
                     "crystal_centroid": [round(float(v), 3) for v in cen]}
    out["box_offset_from_ligand"] = round(offset, 2)

    encl, frac, closest = receptor_encloses(receptor_pdbqt, atoms)
    out["receptor_encloses_site"] = {"fraction_within_5A": frac, "closest_atom": closest}
    if not encl:
        out["verdict"] = "RECEPTOR_SITE_MISMATCH"
        # Which of the two causes it is, measured rather than assumed -- see crystal_enclosure().
        xfrac = crystal_enclosure(pdb, atoms)
        out["receptor_encloses_site"]["crystal_fraction_within_5A"] = xfrac
        if xfrac is not None and xfrac < _ENCLOSE_MIN_FRAC:
            out["site_mismatch_cause"] = "solvent_exposed_site"
            out["reasons"].append(
                f"only {frac:.0%} of {resn}'s atoms lie within 5 A of the prepared receptor "
                f"(closest {closest} A) — but the untouched crystal, all chains included, encloses "
                f"it no better ({xfrac:.0%}), so this is NOT a preparation defect and there is no "
                f"missing chain to add. {resn} is largely solvent-exposed in its own site. No "
                f"redocking control can be constructed here: the pose is held by too little buried "
                f"contact for reproducing it to test anything, which is also why a docking score "
                f"on this site is the least interpretable kind.")
        else:
            out["site_mismatch_cause"] = "receptor_missing_chain"
            out["reasons"].append(
                f"the prepared receptor does not enclose {resn}'s site — only {frac:.0%} of ligand "
                f"atoms lie within 5 A of it (closest {closest} A), while the raw crystal encloses "
                f"it {'%.0f%%' % (xfrac * 100) if xfrac is not None else 'fully'}. The receptor is "
                f"missing the chain the ligand binds to. Prepare the correct chain; this is not a "
                f"verdict on the docking protocol.")
        return out

    is_cov, cov_d, cov_res = covalent_to_protein(pdb, atoms)
    out["covalent"] = {"is_covalent": is_cov, "min_dist_to_nucleophile": cov_d, "residue": cov_res}

    # Diagnostic only -- see metal_coordination() for why this does not decide a verdict.
    co_m, co_d, co_tag = metal_coordination(pdb, atoms)
    if co_d is not None:
        out["metal_coordination"] = {"coordinates_metal": co_m, "min_dist": co_d, "atom": co_tag}
        if co_m:
            out["reasons"].append(
                f"ligand coordinates a cofactor metal ({co_d} A to {co_tag}). Vina scores this as "
                f"vdW, so a poor RMSD here may be a force-field limit rather than a sampling one; "
                f"this does not by itself decide the verdict.")
    if is_cov:
        out["verdict"] = "NOT_APPLICABLE"
        out["reasons"].append(
            f"the co-crystal ligand is COVALENTLY bound ({cov_d} A to {cov_res}). A non-covalent "
            f"redocking control cannot reproduce a covalent pose, so this gate cannot validate "
            f"the protocol here — and a non-covalent screen against this site is itself suspect.")
        return out

    smiles = fetch_smiles(resn)
    if not smiles:
        out["reasons"].append(f"no SMILES in the CCD for {resn}")
        return out
    out["ligand"]["smiles"] = smiles

    lig = tmp / f"{resn}_probe.pdbqt"
    if not ligand_pdbqt_from_smiles(smiles, lig):
        out["reasons"].append(f"ligand prep failed for {resn} — the pipeline's own prep path")
        return out

    pose_out = tmp / f"{resn}_redock_out.pdbqt"
    cx, cy, cz = box_center
    sx, sy, sz = box_size
    cmd = [str(vina_bin), "--receptor", str(receptor_pdbqt), "--ligand", str(lig),
           "--out", str(pose_out),
           "--center_x", str(cx), "--center_y", str(cy), "--center_z", str(cz),
           "--size_x", str(sx), "--size_y", str(sy), "--size_z", str(sz),
           "--exhaustiveness", str(exhaustiveness), "--num_modes", str(num_modes),
           "--seed", str(seed)]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    scores = [float(s) for _, s in _SCORE_RE.findall(proc.stdout)]
    if not scores or not pose_out.exists():
        out["reasons"].append(f"vina produced no poses (rc={proc.returncode})")
        return out

    heavy_raw = [l for l in raw
                 if ((l[76:78].strip() or l[12:16].strip()[0]).upper() != "H")]
    ref_block = "\n".join(heavy_raw) + "\nEND\n"
    poses = []
    for i, blk in enumerate(_pdbqt_poses_as_pdb(pose_out), 1):
        r = best_rms(smiles, blk, ref_block)
        poses.append({"rank": i, "score": scores[i - 1] if i <= len(scores) else None,
                      "rmsd": None if r is None else round(r, 2)})
    got = [p for p in poses if p["rmsd"] is not None]
    if not got:
        # Fall back to element-matched assignment rather than returning ERROR. An unevaluated gate
        # is a bypass, and a lower-bound RMSD still answers "is the top pose anywhere near the
        # crystal answer?" -- if the LOWER bound already exceeds the threshold, the true value
        # does too, so a FAIL from here is sound. A PASS is flagged as unverified.
        out["rmsd_method"] = "element-matched (LOWER BOUND — CCD template did not match)"
        for i, blk in enumerate(_pdbqt_poses_as_pdb(pose_out), 1):
            r = _fallback_rms(blk, ref_block)
            poses[i - 1]["rmsd"] = None if r is None else round(r, 2)
        got = [p for p in poses if p["rmsd"] is not None]
        if not got:
            out["poses"] = poses
            out["reasons"].append("RMSD could not be computed for any pose, by either method")
            return out

    best = min(got, key=lambda p: p["rmsd"])
    out["poses"] = poses
    out["best_rmsd"] = best["rmsd"]
    out["best_rank"] = best["rank"]
    out["top_pose_rmsd"] = poses[0]["rmsd"]

    if offset > box_tolerance:
        out["reasons"].append(
            f"box centre is {offset:.1f} A from the co-crystal ligand centroid "
            f"(tolerance {box_tolerance}) — box-sanity failure, the box may not be on the site")

    # JUDGE THE TOP POSE. A prospective screen consumes rank 1 and nothing else, so that is the
    # thing that has to be right.
    #
    # An earlier version of this gate demanded that the near-native pose BE rank 1. Calibration
    # killed that rule: on 3PTB (trypsin/benzamidine, the easiest redock there is) the near-native
    # pose lands at rank 1, 2 and 2 across seeds 7, 1 and 99 -- while the TOP pose is itself
    # near-native every time (0.41-0.45 A) and the score gap is 0.00-0.03 kcal/mol. Rank was
    # measuring which of two indistinguishable correct poses won an RNG tie, so the gate failed
    # the textbook case two times in three. A gate that rejects benzamidine rejects everything.
    top = poses[0]
    if top["rmsd"] is None:
        out["reasons"].append("top pose RMSD could not be computed — target not evaluated")
    elif top["rmsd"] > pass_rmsd:
        # Both are failures; which one says whether the protocol is worth repairing.
        if best["rmsd"] <= pass_rmsd:
            gap = (best["score"] - top["score"]) if (best["score"] and top["score"]) else None
            out["failure_mode"] = "scoring"
            out["reasons"].append(
                f"top pose is {top['rmsd']} A from the crystal pose. A near-native pose "
                f"({best['rmsd']} A) IS generated at rank {best['rank']} but scores "
                f"{gap:+.2f} kcal/mol worse, so the scoring function prefers a wrong pose. More "
                f"sampling cannot fix this — it optimises harder against the same objective.")
        else:
            out["failure_mode"] = "sampling"
            out["reasons"].append(
                f"top pose is {top['rmsd']} A and no pose reaches {pass_rmsd} A "
                f"(best {best['rmsd']} A) — sampling never visits the known binding mode.")

    out["verdict"] = "FAIL" if out["reasons"] else "PASS"
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--receptor", required=True)
    ap.add_argument("--pdb", required=True)
    ap.add_argument("--center", nargs=3, type=float, required=True)
    ap.add_argument("--size", nargs=3, type=float, default=[24, 24, 24])
    ap.add_argument("--ligand-resname", default=None)
    ap.add_argument("--chain", default=None)
    ap.add_argument("--exhaustiveness", type=int, default=16)
    ap.add_argument("--vina", default=os.environ.get("VINA_BIN", "vina"))
    ap.add_argument("--pass-rmsd", type=float, default=2.0)
    ap.add_argument("--workdir", default=None)
    ap.add_argument("--json", action="store_true", help="emit the verdict dict as JSON")
    a = ap.parse_args()

    v = redock_control(a.receptor, a.pdb, a.center, a.size, ligand_resname=a.ligand_resname,
                       chain=a.chain, exhaustiveness=a.exhaustiveness, vina_bin=a.vina,
                       pass_rmsd=a.pass_rmsd, workdir=a.workdir)
    if a.json:
        print(json.dumps(v, indent=2))
    else:
        lig = v.get("ligand", {})
        print(f"{v['pdb_id']}  ligand {lig.get('resname','?')} chain {lig.get('chain','?')} "
              f"({lig.get('n_heavy','?')} heavy atoms)")
        print(f"box offset from ligand centroid: {v.get('box_offset_from_ligand','?')} A")
        for p in v.get("poses", []):
            print(f"   rank {p['rank']}: score {p['score']}  RMSD {p['rmsd']} A")
        print(f"\nVERDICT: {v['verdict']}")
        for r in v["reasons"]:
            print(f"  - {r}")
    return 0 if v["verdict"] == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())


# --------------------------------------------------------------------------- receptor prep

AD4_TYPES = {"C": "C", "N": "NA", "O": "OA", "S": "SA", "H": "H", "F": "F", "CL": "Cl",
             "BR": "Br", "P": "P", "FE": "Fe", "ZN": "Zn", "MG": "Mg", "CA": "Ca", "MN": "Mn"}


def prepare_receptor(pdb_path: Path, out_path: Path, chains_remove=None,
                     keep_cofactors: bool = True) -> Path:
    """Minimal PDB -> PDBQT receptor: polymer atoms plus any cofactor in the _COFACTOR allowlist.

    Production (`run_acrb_efflux_screen.pdb_to_pdbqt_receptor`) keeps only ATOM records, and this
    function used to copy that convention on the argument that the gate should test the receptor
    the pipeline actually docks against. That argument fails when the convention is the defect:
    dropping 3L4D's haem leaves an apo cavity, and fluconazole -- which binds by coordinating the
    haem iron at 1.88 A -- redocks 9.02 A away. The gate was certifying a receptor nobody should
    dock into, and would have rejected a correctly built one.

    `keep_cofactors=False` restores the old behaviour, for measuring what the cofactor is worth.
    NOTE the production prep is still ATOM-only, so for a cofactor target this gate now tests a
    BETTER receptor than production builds; that gap is recorded, not closed, by this change.
    """
    drop = {c.upper() for c in (chains_remove or [])}
    lines = []
    kept_cof: dict[str, int] = {}
    bad_types: dict[str, set] = {}
    for line in pdb_path.read_text().splitlines():
        rn = line[17:20].strip().upper()
        is_atom = line.startswith("ATOM")
        is_cof = (line.startswith("HETATM") and keep_cofactors
                  and (rn in _COFACTOR or rn in _METAL_ION))
        if not (is_atom or is_cof):
            continue
        if line[21].upper() in drop or rn in ("HOH", "WAT", "DOD"):
            continue
        el = (line[76:78].strip().upper() if len(line) >= 78 else "")
        if not el:
            # Fallback from the atom name. One character for polymer atoms -- " CA " in a residue
            # is alpha-carbon, and reading it as calcium would retype the whole backbone. Two-
            # character elements are only trusted for a standalone ion, where the residue name IS
            # the element; otherwise an FE atom silently types as F (fluorine) and the haem docks
            # as four nitrogens and a halogen, which looks exactly like it worked.
            nm = line[12:16].strip().lstrip("0123456789").upper()
            el = nm if (is_cof and rn in _METAL_ION and nm == rn) else nm[:1]
        if is_cof:
            kept_cof[rn] = kept_cof.get(rn, 0) + 1
            # Fail closed on an element we cannot type. Falling back to el[:1] on a cofactor is how
            # Fe becomes F: the receptor still docks, the run still completes, and the result is
            # indistinguishable from a good one. Checked only for cofactors -- the ATOM branch is
            # unchanged and its verdicts are the control for this whole change.
            if el not in AD4_TYPES:
                bad_types.setdefault(rn, set()).add(el or line[12:16].strip())
        lines.append(f"{line[:66].ljust(66)}    {0.000:6.3f} {AD4_TYPES.get(el, el[:1]):<2}")
    if bad_types:
        raise ValueError(
            f"{pdb_path.name}: untypeable atoms in retained cofactor(s) "
            f"{ {k: sorted(v) for k, v in bad_types.items()} } — add them to AD4_TYPES rather than "
            f"letting them fall back to a one-character element guess.")
    out_path.write_text("\n".join(lines) + "\n")
    if kept_cof:
        types: dict[str, int] = {}
        for l in lines:
            if l[17:20].strip().upper() in kept_cof:
                types[l[77:].strip()] = types.get(l[77:].strip(), 0) + 1
        out_path.with_suffix(".cofactors.json").write_text(
            json.dumps({"residues": kept_cof, "ad4_types": types}, indent=2))
    return out_path
