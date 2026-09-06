# Pre-registration: proper reconstruction of the pose-ensemble arm

**Written 7 September 2026, before any reconstructed estimate is read.**

## The problem this has to solve

`RESULT.md` published **+0.0194 [−0.0169, +0.0548]** and called the arm refuted. The script that
produced it was never committed (`git log --all` over `analysis/pose_ensemble/` returns only the
pre-registration and the write-up). A reimplementation from the prose description gives
**+0.0393 [+0.0054, +0.0733]** on the identical poses. The two cannot be reconciled by inspection
because one side no longer exists.

Matching the original is therefore not achievable and is **not** the goal. The goal is to
establish whether the arm's effect is stable enough to support *any* verdict.

## The hypothesis being tested

A diagnostic points at the cause. The descriptor baseline — same seven descriptors, same 751
compounds, a deterministic function of the fold assignment — is published at **0.7651** and my
reimplementation gives **0.7709**. That is a pure fold-assignment difference, and it is a fifth
of the gap between the two ensemble estimates.

**H1: the ensemble effect is unstable under cross-validation fold assignment, and +0.0194 and
+0.0393 are both ordinary draws from its sampling distribution.**

If H1 holds, neither the published refutation nor the reimplementation's positive reading is
readable, and the arm is unresolved for a reason that has nothing to do with the receptor.

## Method

The identical feature construction already committed in `analyse_protonated.py`, run across
**40 cross-validation fold seeds** on both pose sets. Reported as the distribution of
ΔAUROC(ensemble − top pose), not a point estimate. Seeds are `20260907 + i`, fixed here.

**No seed selection.** Every seed drawn is reported. Under no circumstances is a seed chosen for
reproducing a target value; that would be fitting the analysis to the answer.

## Pipeline controls — all must pass before any effect is read

Registered because the previous run's control validated the docking and not the feature pipeline,
which is where the risk was.

1. **Degenerate identity.** The ensemble model restricted to its single top-score feature must
   reproduce the top-pose arm to within 1e-9 AUROC, per seed. If the pipeline cannot reproduce a
   result it contains as a special case, nothing else it produces is readable.
2. **Descriptor baseline in range.** Across the 40 seeds the descriptor baseline must span the
   published 0.7651, confirming the fold machinery is the same kind of object as the original's.
3. **Permuted labels.** Shuffling labels must give mean AUROC in [0.45, 0.55].
4. **Permuted features.** Shuffling the ensemble feature rows must not produce a positive gain.

## Reading rules — fixed now

Let *lo*, *hi* be the 2.5th and 97.5th percentiles of ΔAUROC across the 40 seeds.

- **UNRESOLVED (H1 supported)** — the seed range spans the scoring floor 0.0201, i.e.
  *lo* < 0.0201 < *hi*. The estimate is not stable enough to place the arm relative to the floor,
  and both prior numbers are draws from it. The manuscript reports the arm as unresolved.
- **RESOLVABLE, above floor** — *lo* > 0.0201. The arm clears regardless of fold assignment and
  the published refutation was wrong.
- **RESOLVABLE, below floor** — *hi* < 0.0201. The published refutation stands and the
  reimplementation's positive reading was a fold artefact.

Reported for both receptors. The receptor effect is the difference between the two distributions
and is only interpretable if both are resolvable.

## Pre-committed disclosure

Whatever this returns is reported and the manuscript is corrected to match, including if the
outcome is that a published result of ours is unreproducible and unresolvable. The reconstruction
script is committed **with** its results, so this cannot recur.
