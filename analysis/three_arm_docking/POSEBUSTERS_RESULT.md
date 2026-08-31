# PoseBusters: DiffDock's chemistry is fine. It puts the ligand too close to the protein.

**31 August 2026.** `run_posebusters.py`, 60 pocket-constrained poses, full-protein receptor,
results in `posebusters_pocket.csv`.

## I predicted the wrong failure mode

When 39% of pocket-constrained poses scored as clashes, I wrote that a diffusion model "is not
minimising a force field, so **bond lengths, angles and contacts** need not satisfy one." The
first two are wrong.

| check family | pass rate |
|---|---|
| bond lengths | **100%** |
| bond angles | **100%** |
| aromatic ring flatness | **100%** |
| non-aromatic ring non-flatness | **100%** |
| double bond flatness | **100%** |
| internal energy | **100%** |
| internal steric clash | 95% |
| sanitization / connectivity / no radicals | 100% |
| **minimum distance to protein** | **47%** ← the only real failure |
| volume overlap with protein | 97% |

**PB-VALID (every check passes): 28/60 = 47%.**

DiffDock's **internal** chemistry is essentially perfect — every bond length, every angle, every
ring. What fails is **intermolecular**: 53% of poses sit closer to the protein than van der Waals
contact allows.

| | all checks pass |
|---|---|
| intramolecular (ligand geometry alone) | **57/60** |
| intermolecular (ligand vs protein) | **28/60** |

## Why this matters for the next step

It sharpens the minimisation prediction rather than weakening it. Local minimisation in the
receptor field relieves **exactly** this failure — it pushes a ligand out of a steric overlap —
and there is nothing left for it to fix internally, because the internal geometry already passes.
So the minimise-and-rescore pass, already queued, is now the right experiment for a precise
reason instead of a plausible one.

It also predicts the shape of the result: if minimisation works, PB-validity should rise toward
the intramolecular ceiling of ~95%, and the Vina clash rate should fall from 39% toward zero.
If it does not, the poses are in positions that cannot be relieved locally, and arm 2 stays
unreportable.

## Is the truncated receptor to blame?

Checked, not assumed. Poses were generated against a receptor truncated to residues with any atom
within 18 Å of the box centre, then scored against the full protein — so a re-added atom could in
principle clash. It cannot account for this: pose centroids sit a median 3.6 Å from centre and at
most 10.4 Å, ligand radii add roughly 5 Å, so ligand atoms reach ~16 Å, while every removed
residue has all atoms beyond 18 Å. The ≥2 Å margin is far outside the contact distance PoseBusters
flags. The clashes are with residues that were present during generation.

## Against the published benchmark

PoseBusters' paper reports deep-learning docking producing **up to 95% invalid poses**. Ours are
**53%** invalid, and only intermolecularly. Same phenomenon, milder, and now attributed to a
specific check rather than inferred from a positive Vina score.

That the improvised proxy (39% clashing) and the proper test (53% failing minimum-distance) agree
in direction and rough magnitude is a useful cross-check on both.

## Correction recorded

`ARM2_TWO_FAILURES.md` attributes the residual to "bond lengths, angles and contacts". Only
*contacts* is right. The other two are measured at 100% and that sentence overstated the case
against generative pose geometry.
