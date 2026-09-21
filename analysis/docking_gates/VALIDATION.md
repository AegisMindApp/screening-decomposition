# Redocking gate — calibration and validation record

## Why the first criterion was wrong

The gate originally required the near-native pose to **be rank 1**. That rule fails
trypsin/benzamidine (3PTB), the easiest redocking case in structural biology, two times in three:

| seed | near-native RMSD | its rank | **top-pose RMSD** | score gap |
|---|---|---|---|---|
| 1 | 0.42 A | 2 | **0.43 A** | +0.030 kcal/mol |
| 7 | 0.41 A | 1 | **0.41 A** | +0.000 |
| 99 | 0.43 A | 2 | **0.45 A** | +0.029 |

The top pose is *itself* near-native every time. Rank was recording which of two indistinguishable
correct poses won an RNG tie, at a score separation of 0.03 kcal/mol against a Vina RMSE of ~2.8.

**A gate that rejects benzamidine rejects everything**, so it would have been switched off within a
week and the pipeline would have gone back to ungated docking.

Criterion is now **top-pose RMSD <= 2.0 A**, which is what a prospective screen consumes: it reads
rank 1 and nothing else.

## Seed robustness of the revised criterion (exhaustiveness 16, seeds 1/7/42/99)

| target | verdict | top-pose RMSD across seeds | failure mode |
|---|---|---|---|
| 3PTB trypsin/benzamidine | **PASS** stable | 0.43, 0.41, 0.43, 0.45 | — |
| 1STP streptavidin/biotin | **FAIL** stable | 6.00, 5.99, 5.98, 5.99 | scoring |
| 3AI8 cathepsin B | **FAIL** stable | 2.61, 2.61, 2.60, 2.60 | scoring |
| 3N4C cathepsin S | **NOT_APPLICABLE** stable | — | covalent native ligand |

No verdict flips on seed, and the failing RMSDs are reproducible to 0.01-0.02 A — these are
properties of the protocol, not sampling noise.

**The gate can both pass and fail**, which is the only thing that makes a PASS informative. A suite
where everything fails cannot distinguish a working gate from a broken one.

## sampling vs scoring

Both are failures; the distinction says whether a protocol is worth repairing.

- **sampling** — no pose reaches the threshold. More exhaustiveness might help.
- **scoring** — a near-native pose *is* generated but something worse outranks it. More
  exhaustiveness cannot help; it optimises harder against the same objective. Observed directly on
  CTSS, where raising exhaustiveness 16 -> 128 *improved* the top score (-9.081 -> -9.130) while it
  stayed ~7 A wrong.

1STP is the cleanest example: biotin/streptavidin is among the tightest known non-covalent
complexes, and Vina puts its top pose 6 A from the crystal answer while generating a 1.36 A pose at
rank 6.

## Failure modes of the gate itself (all found by testing, all silently disabling)

1. `Chem.RemoveHs` keeps stereo-defining explicit hydrogens, so CCD templates came out 34 atoms
   against a 33-atom PDB block and **every** RMSD returned `None` — surfacing as "could not
   compute", i.e. a bypass rather than an error. Fixed with `RemoveAllHs`.
2. `MolFromPDBBlock(removeHs=True, sanitize=False)` does not remove hydrogens; a crystal ligand
   with riding H gave a 58-atom block against a 35-atom template.
3. An `ERROR` verdict left a target unevaluated. Now falls back to element-matched RMSD, which is a
   **lower bound** — so a FAIL from it is sound, and a PASS is explicitly flagged unverified.
4. My own CLI test passed `--center 0 0 0` and produced a "box-sanity failure" that was an artefact
   of the input, not the target.

The pattern is worth stating: **every bug in a gate defaults to passing.** A gate must be tested on
a case it should reject *and* one it should accept.

## Covalent ligands

Detected by proximity to protein nucleophiles (SG/OG/OG1/NE2/SD/OH, cutoff 1.95 A). Verified not to
over-fire: 3PTB 3.04 A and 1STP 2.58 A both read non-covalent; 3N4C reads 1.77 A to CYS25A.

Cathepsins are cysteine proteases and most of their inhibitors are covalent, so a **non-covalent**
redocking control cannot validate a protocol there — and a non-covalent *screen* against such a
site is itself questionable. That is the finding that retrospectively invalidated the CTSS
repurposing screen: it scored non-covalent CNS drugs against a covalent-inhibitor site.


## The manifest, and a false alarm I raised on the way to it

First build reported **three FAILs** (MSH3 12.19 A, elastase 8.44 A, trypsin 11.77 A) and I began
reporting that "every receptor in the pipeline fails". **That was wrong**, and the discriminator
that caught it was already in the suite: trypsin/benzamidine (3PTB) PASSes at 0.43 A while the
pipeline's *own* trypsin receptor (1TRN) failed at 11.77 A. Same protein. A protocol does not
work on one trypsin structure and miss by 11 A on another.

Two bugs in ligand auto-detection, both of which manufacture failures:

1. **Modified amino acids were selected as ligands.** 1TRN's PTR is phosphotyrosine and 1HNE's
   ALV sits inside a peptidyl inhibitor; both carry full N/CA/C/O backbones. Redocking a residue
   that is covalently in the chain measures nothing, because no free-ligand pose exists to
   reproduce. Now filtered on backbone-atom presence.
2. **Covalent detection only looked at named nucleophiles** (SG/OG/...). 1HNE's fragment is bonded
   to an adjacent modified residue at 1.31 A, not to a catalytic side chain, so it read as free.
   Now checks all protein heavy atoms. The cutoff separates cleanly with a wide margin:

   | | closest protein atom | |
   |---|---|---|
   | 3PTB benzamidine | 2.82 A | free |
   | 1STP biotin | 2.58 A | free |
   | 3THW ADP | 2.75 A | free |
   | 3N4C EF3 nitrile | 1.77 A | covalent |
   | 3RXX / 1TRN | 1.59 A | covalent |
   | 1HNE peptidyl | 1.31 A | covalent |

3. A receptor that does not enclose its ligand's site now gets its own verdict rather than a
   meaningless RMSD. 3THW's ADP sits in chain A while the MSH3 receptor is chain B only, so the
   "12.19 A sampling failure" was measuring the distance between two chains.

**Corrected manifest** — no receptor passes, but for reasons that are specific and actionable
rather than a blanket condemnation:

| receptor | verdict | why |
|---|---|---|
| 3RXX beta-lactamase | NOT_APPLICABLE | boronic acid, covalent to SER70 |
| trypsin 1TRN | NOT_APPLICABLE | organophosphate, covalent to SER195 |
| elastase 1HNE | NOT_APPLICABLE | peptidyl inhibitor, covalent |
| MSH3 3THW | RECEPTOR_SITE_MISMATCH | ADP is in chain A; the receptor is chain B only |

Three are serine-hydrolase structures solved with covalent inhibitors, so they cannot validate a
*non-covalent* protocol — and by the same token a non-covalent screen against them is
questionable. The MSH3 entry is **fixable**: prepare the chain the ligand actually binds. It is
also the independent confirmation of the known MSH3 defect, that the box sat off-pocket for
chain B, reached here without any knowledge of that history.


## Re-validation after the detection fixes (17 Aug 2026)

The suite above was run BEFORE the modified-residue filter and the all-heavy-atom covalent check
went in, so its results were produced by superseded code and could not be relied on. Re-run on
seeds 1/42/99:

| target | verdict | top-pose RMSD | unchanged? |
|---|---|---|---|
| 3PTB benzamidine | PASS stable | 0.43, 0.43, 0.45 | yes |
| 1STP biotin | FAIL (scoring) stable | 6.00, 5.98, 5.99 | yes |
| 3AI8 cathepsin B | FAIL (scoring) stable | 2.61, 2.60, 2.60 | yes |
| 3N4C cathepsin S | NOT_APPLICABLE stable | — | yes |

Every verdict is unchanged, which is the expected result: the fixes only alter which HETATM is
selected as the ligand, and these four were already selecting the right one. They mattered for the
pipeline receptors (1TRN, 1HNE, 3THW), where the wrong selection was fabricating 8-12 A failures.

The check that had to pass was 3PTB: the generalised covalent test inspects every protein heavy
atom, and benzamidine's closest contact is 2.82 A. Had the cutoff crept above that, the gate's only
positive control would have become NOT_APPLICABLE and the gate would have lost its ability to
return PASS at all — failing closed on everything, undetectably.
