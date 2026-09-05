# The resolution floor was measuring two things at once

**6 September 2026.** Found while bootstrapping an interval for the floor before submission.
`floor_interval.py`, results in `FLOOR_INTERVAL.json`. The positive control failed, which is why
this file exists.

## What MARGINAL_VALUE.md §3 says

> gnina's `--score_only` re-evaluation of **the identical Vina poses** reproduces the cached Vina
> score at r = 0.9369 … AUROC from cache 0.4182, AUROC from gnina 0.3793 — 0.039 apart.

## What it actually compared

The numbers reproduce exactly — but only when the cached side is the **exhaustiveness-32** run,
while gnina rescored the **exhaustiveness-4** poses. Those are two independent stochastic
searches. The poses are not identical.

| | cached side | gnina side | n | r | gap |
|---|---|---|---|---|---|
| **as published** | exh **32** cache | `--score_only` on exh **4** poses | 745 | 0.9372 | **+0.0393** |
| **identical poses** | exh **4** cache | `--score_only` on the same exh 4 poses | 743 | 0.9700 | **+0.0201** |

Control targets from MARGINAL_VALUE.md §3 — r 0.9369, mean −0.043 kcal/mol, 0/745 above
2 kcal/mol, AUROC 0.4182 / 0.3793 — are reproduced by the first row to within 4×10⁻⁴ and by the
second row not at all.

## Both numbers are real floors. They gate different comparisons.

| | value | 95% CI | what it bounds |
|---|---|---|---|
| **Protocol floor** | **0.0393** | [+0.0235, +0.0557] | run-to-run reproducibility of the whole docking protocol — search included. The right bar for comparing two *protocols*, which is what most published improvements are |
| **Scoring floor** | **0.0201** | [+0.0112, +0.0286] | numerical and implementation noise with placement held fixed. The right bar for comparing two *scoring functions* on shared poses |

Paired bootstrap over compounds, 10,000 resamples, unprotonated receptor.

## Consequences

The single 0.04 bar was applied to all six interventions regardless of whether each held poses
fixed. Re-gated against the appropriate floor, **every qualitative reading survives** and one
improves:

| intervention | Δ | appropriate floor | before | after |
|---|---|---|---|---|
| Scoring function (gnina, identical poses) | +0.140 | scoring, 0.020 | yes | **yes** |
| Receptor repair (re-docked) | +0.045 | protocol, 0.039 | marginal | **marginal — clears, barely** |
| Pose ensemble (same run, different selection) | +0.019 | scoring, 0.020 | no | **no — sits on the floor** |
| Pose generator, DiffDock-L (re-placed) | +0.017 | protocol, 0.039 | no | **no** |
| Pocket conditioning, Boltz-2 (different model) | −0.013 | protocol, 0.039 | no | **no** |
| Eight-fold search effort (exh 4 → 32) | −0.019 | protocol, 0.039 | no | **no** |

The last row is now the sharpest line in the analysis: **the measured effect of 8× search effort
is smaller than the run-to-run reproducibility of the protocol it was measured in.** The exh4/exh32
pair *is* the protocol floor — so that intervention is, quite literally, indistinguishable from
running the same protocol twice.

## What this changes in the write-ups

`MARGINAL_VALUE.md` §3's phrase "the identical Vina poses" is wrong and is corrected there. The
0.039 value stands, with its meaning restated. Nothing that was called resolvable becomes
unresolvable, and nothing that was called unresolvable becomes resolvable.

## Why the control caught it

The floor had never been recomputed from the raw score files since it was first written down; it
was quoted forward. Reproducing it from `exh4_scores_merged.json` gave 0.0201, and the discrepancy
against the published 0.0393 was only explicable by testing which cache the published number came
from. A number that is quoted forward and never re-derived is indistinguishable from a correct one.
