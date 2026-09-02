# Pose-ensemble scoring: +0.0194 AUROC, inside the noise floor — REFUTED

**2 September 2026.** Reading rules fixed in `PREREGISTRATION.md`, committed as `448d954e5`
**before** this number existed. 751 Mpro compounds, 5 Vina modes each, exhaustiveness 4,
original (donor-defective) receptor. All arms share identical compounds and identical
StratifiedKFold(5) folds, out-of-fold predictions only.

Ensemble features: Boltzmann-weighted ΔG at 298 K, mean and SD of pose scores, gap to the
second pose, count of poses within 1 kcal/mol of the top, and RMSD spread of those
near-degenerate poses (well width).

| arm | AUROC |
|---|---|
| top pose only (current practice) | 0.5895 |
| **ensemble (8 features)** | **0.6089** |
| permuted control | 0.5658 |
| descriptors | 0.7651 |
| descriptors + ensemble | 0.7596 |

| pre-registered rule | result | 95% CI | verdict |
|---|---|---|---|
| **1.** ensemble − top pose > 0.04 | **+0.0194** | [−0.0169, +0.0548] | **FAIL** — below floor, CI spans 0 |
| **3.** permuted control shows no gain | −0.0237 | [−0.0511, +0.0019] | **PASS** — control clean, no leak |
| **2.** desc+ens − descriptors > 0.04 | −0.0056 | [−0.0175, +0.0060] | **FAIL** — adds nothing |

## Verdict: refuted

The point estimate is **positive and in the predicted direction** — the pose distribution does
appear to carry a little information the top pose discards. But +0.0194 sits below the ~0.04
measurement floor and its interval spans zero, which is precisely the case Rule 1 was written
in advance to call refuted rather than promising. Adding ensemble features to descriptors
makes the model slightly *worse*.

The permuted-feature control behaved correctly (no gain, slightly negative), so the pipeline
is not leaking and both readings stand.

## One number that needs explaining

"Top pose only" scores 0.5895 here, against 0.4451 for raw Vina ranking on the same receptor.
The difference is that this arm passes the score through logistic regression, which is free to
**fit the sign** — and Vina's relationship to activity on Mpro is inverted. This is the
"Vina as a fitted feature" effect seen earlier (0.5692). It makes the comparison in this table
fitted-vs-fitted, which is the fair test for this hypothesis, but it means 0.5895 is not
comparable to the raw-ranking numbers elsewhere in this repo.

## Consequence

This was the last direction from the external review that our existing data could not already
refute. With it refuted, all four are closed:

| direction | outcome |
|---|---|
| Residual docking | Reverses sign between targets, both significant |
| Disagreement-as-signal | 0.5321, not the predicted 0.60+; inside the noise floor |
| Few-shot target adaptation | Pre-empted — descriptors already are a target-specific fit, 3D adds nothing |
| **Pose-ensemble scoring** | **+0.0194, inside the noise floor** |

The docking-methodology line is closed. Not because the ideas were bad — the pose-ensemble
point estimate was positive — but because every candidate improvement lands inside a
measurement floor of ~0.04 AUROC that this benchmark cannot see past.
