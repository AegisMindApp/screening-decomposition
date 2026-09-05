# Arm 2: pocket conditioning does not help Boltz-2. It slightly hurts.

**5 September 2026.** 751 dispatched, 750 scored, 1 failure, 20.4 GPU-hours on a GCE L4.
Rules fixed in `PREREGISTRATION_ARM2.md` (`e96d4fea1`) before the run. Paired against arm 1 on
the 750 compounds both arms scored.

## Result

| | AUROC |
|---|---|
| seven free descriptors | 0.7654 |
| **arm 1 — blind** | **0.7913** |
| **arm 2 — pocket-conditioned** | **0.7783** |
| descriptors + arm 1 | 0.8588 |
| descriptors + arm 2 | 0.8435 |

| pre-registered rule | Δ | 95% CI | verdict |
|---|---|---|---|
| **2.** arm 2 − arm 1 — *the question this arm exists to answer* | **−0.0130** | [−0.0264, +0.0000] | **NOT DEMONSTRATED** |
| **1.** arm 2 − descriptors | +0.0129 | [−0.0370, +0.0614] | NOT DEMONSTRATED |
| **3.** desc + arm 2 − descriptors | +0.0781 | [+0.0473, +0.1089] | SUPPORTED |

**Arm 1 vs arm 2 score correlation: r = 0.952.** Supplying the binding site barely changed the
predictions at all.

The pocket was not a guess and not mis-specified: 14 residues within 8 Å of the Vina box
centre, containing the catalytic dyad **Cys145 (4.41 Å)** and **His41 (5.02 Å)**, supplied as
1-based sequence indices after the PDB-numbering bug was found and the run restarted. The log
confirms `pocket from bundle: 14 residues [('A', 162), ('A', 143), ('A', 39), ...]` — 143 is
Cys145, 39 is His41.

## This retracts a caveat I attached to arm 1

Arm 1's write-up recorded that it ran **blind** against a pocket-constrained Vina comparator,
called that a handicap, and described its +0.0934 as a **lower bound**. That caveat was
recorded honestly and in advance — and it is now shown to be unnecessary. Running blind cost
Boltz-2 nothing; if anything the blind arm scored slightly higher. **Arm 1's result stands on
its own terms and needs no allowance.**

## Where this leaves the decomposition

Pocket conditioning joins the set of interventions that do not survive this benchmark:

| intervention | ΔAUROC |
|---|---|
| **Changing what does the scoring** (descriptors + Boltz-2) | **+0.093** |
| Repairing the receptor | +0.045 |
| Pose ensemble vs top pose | +0.019 |
| Swapping the pose generator | +0.017 |
| **Supplying the binding site to Boltz-2** | **−0.013** |
| Eight times the search effort | −0.019 |

Five interventions in the "help the model find the site" family — pose source, search effort,
pose ensembles, receptor preparation, and now pocket conditioning — and only receptor repair
clears the noise floor, at +0.045. The one thing that consistently works is changing what does
the scoring.

## Cost

20.4 GPU-hours, ~AUD 27, inside the authorised 25.0 GPU-hour ceiling. Instance verified
TERMINATED. Cumulative Boltz-2 GCE spend ~AUD 54 against the AUD 100 project budget; Factor Xa
is running free on Kaggle.
