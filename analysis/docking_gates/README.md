# Redocking gate

Whether a docking protocol can put a co-crystallised ligand back where the crystal says it is,
on the receptor and at the exhaustiveness a screen actually uses. A screen against a receptor
that fails this is not a weak screen; it is an unvalidated one.

`VALIDATION.md` is the record. Three things in it are the point:

**The first criterion was wrong, and the positive control is what caught it.** The gate originally
required the near-native pose to *be rank 1*. That fails trypsin/benzamidine (3PTB) — the easiest
redocking case in structural biology — two times in three, because rank was recording which of two
indistinguishable correct poses won an RNG tie at 0.03 kcal/mol against a Vina RMSE of ~2.8. A gate
that rejects benzamidine rejects everything, and would have been switched off inside a week. The
criterion is now top-pose RMSD ≤ 2.0 Å, which is what a prospective screen consumes.

**It has to be able to both pass and fail.** Across seeds 1/7/42/99 at exhaustiveness 16: 3PTB
PASSes at 0.43/0.41/0.43/0.45 Å, 1STP (streptavidin/biotin) FAILs at 6.00/5.99/5.98/5.99 Å, 3AI8
(cathepsin B) FAILs at 2.61/2.61/2.60/2.60 Å. No verdict flips on seed and the failing RMSDs
reproduce to 0.01–0.02 Å, so these are properties of the protocol rather than sampling noise. A
suite where everything fails cannot be distinguished from a broken gate.

**There is a third verdict, and it is load-bearing.** Cathepsin S (3N4C) returns NOT_APPLICABLE,
not FAIL: its native ligand is covalent to the catalytic Cys25 (measured 1.77 Å), and a
non-covalent scoring function cannot represent that binding mode at all. A FAIL there would read
as "this protocol is bad here" when the truth is "this question cannot be asked this way". The
same applies to 3RXX, whose ligand is 1.59 Å from Ser70.

## Files

| file | what it is |
|---|---|
| `VALIDATION.md` | the calibration and validation record |
| `gate_validation.json` | the seed-robustness data behind that record |
| `gate_manifest.json` | per-receptor verdicts, with the md5 of the file each was computed on |
| `redock_control.py` | the control itself |
| `build_gate_manifest.py`, `gate_manifest.py`, `validate_gate.py` | manifest construction and lookup |

The receptor PDBQT files the manifest refers to are not in this repository; paths in
`gate_validation.json` are repo-relative to the private tree the gate was run in. `receptor_md5`
in the manifest is what makes a verdict traceable to a specific file — gating one receptor and
then bundling another is two different operations until someone compares hashes, which is a
failure this project has recorded.
