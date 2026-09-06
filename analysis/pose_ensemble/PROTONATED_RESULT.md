# The pose-ensemble arm clears the floor on the repaired receptor — but the published null is not reproducible

**7 September 2026.** `analyse_protonated.py`, results in `PROTONATED_RESULT.json`. Rules fixed
in `PREREGISTRATION_PROTONATED.md` (`6bac07fd9`) before the re-dock started. 751/751 docked on
three Kaggle CPU shards, **zero failures**, all shards reporting `complete`.

**Read the second section before quoting the first.**

## Pre-registered controls — both pass

| control | result | |
|---|---|---|
| Repaired top-mode raw AUROC reproduces published 0.4530 (±0.010) | **0.4490** | **PASS** |
| Coverage ≥ 90% | **100.0%** (751/751) | **PASS** |
| Cross-platform replicate: 40 compounds docked here *and* on Kaggle | r = **0.9842**, mean \|Δ\| 0.101 kcal/mol | consistent |

## Primary result

| | repaired | defective |
|---|---|---|
| raw Vina top-pose AUROC | 0.4490 | 0.4079 |
| top pose (logistic) | 0.5468 | 0.5918 |
| **ensemble (8 features)** | **0.6043** | **0.6313** |
| permuted control | 0.5208 | 0.5085 |
| descriptors | 0.7709 | 0.7709 |
| descriptors + ensemble | 0.7663 | 0.7780 |

| rule | repaired | defective |
|---|---|---|
| **1.** ensemble − top pose | **+0.0576 [+0.0112, +0.1040]** | +0.0393 [+0.0054, +0.0733] |
| permuted control | −0.0251 [−0.0889, +0.0369] | −0.0828 [−0.1463, −0.0210] |
| **3.** desc+ens − descriptors | −0.0047 [−0.0173, +0.0073] | +0.0070 [−0.0060, +0.0200] |

By the rule as written — Δ > 0.0201 *and* CI excluding zero — rule 1 on the repaired receptor
reads **SUPPORTED**. Rule 3 fails on both: the ensemble still adds nothing over descriptors.

## Why that verdict cannot be attributed to the receptor

**The published defective-receptor number is +0.0194 [−0.0169, +0.0548]. Running identical code
on the identical poses, I get +0.0393 [+0.0054, +0.0733].** Twice the effect, and an interval
that excludes zero where the published one spans it.

**The original analysis script was never committed.** `git log --all` over
`analysis/pose_ensemble/` returns only `PREREGISTRATION.md` and `RESULT.md`. The published
+0.0194 was produced by code that no longer exists anywhere in the repository, so the discrepancy
cannot be diagnosed. The eight features here are reconstructed from `RESULT.md`'s prose
description — Boltzmann-weighted ΔG at 298 K, top score, mean, SD, gap to second pose, count
within 1 kcal/mol, RMSD spread of those, pose count — which is a specification, not an
implementation.

So there are two candidate explanations for +0.0576 and they cannot be separated here:

1. Repairing the receptor genuinely lifts the arm above the floor. Under this implementation the
   repair effect is +0.0576 − 0.0393 = **+0.018**, which is *below* the scoring floor.
2. The reconstruction differs from the original method, and the arm was above the floor all
   along on both receptors.

The second is at least as likely as the first, because the reconstruction already disagrees with
the original by +0.020 on data where the receptor is held constant.

## What my pre-registration got wrong

The positive control validated the **docking** (top-mode raw AUROC) and passed cleanly. It did
not validate the **feature pipeline**, which is where the risk actually was. A control that
passes on the part that was never in doubt is the failure mode this project has hit repeatedly:
it looks like verification and is not.

The control that should have been registered: reproduce the published +0.0194 on the defective
poses before reading anything from the repaired ones. Had it been, this run would have halted at
that step instead of producing a verdict.

## Consequence for the manuscript

§3.2 currently reports the pose-ensemble arm as **+0.019, not resolvable**, citing a
pre-registered refutation. That row rests on unreproducible code and must not stand as written.
The honest replacement is that the arm is **unresolved**: the published null cannot be
reproduced, a documented reimplementation places the effect above the scoring floor on both
receptors, and the two cannot be reconciled without the original script.

This moves the paper's headline from "four of six fall inside the floor" to **three of six, one
unresolved** — a weaker claim than currently published, and it must be corrected before
submission rather than after.

## What is solid

- The docking itself: 751/751, zero failures, cross-platform r = 0.9842.
- Rule 3 replicates on both receptors — pose-ensemble features add nothing over seven free
  descriptors, which was the more practically important of the two published findings.
- The permuted control behaves correctly on both arms: permuting the features destroys the gain.
