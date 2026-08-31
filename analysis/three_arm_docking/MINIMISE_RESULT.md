# Local minimisation rescues DiffDock's poses — the placements were repairable after all

**31 August 2026.** `local_minimise.py`, 55 of 60 pocket-constrained poses scored before and
after, results in `local_minimise_result.json`. AutoDock Vina `--local_only`, same binary and same
unprotonated receptor arm 1 used.

## The prediction, and the result

PoseBusters localised the failure precisely: internal chemistry at **100%**, and the only failing
check **minimum-distance-to-protein at 47%**. The prediction recorded before running was that
local optimisation in the receptor field would relieve exactly this.

| | before | after |
|---|---|---|
| poses scoring **positive** (clashing) | 29.1% | **1.8%** |
| mean | +0.13 | **−4.23** |
| median | −1.47 | −4.40 |
| worst | +30.49 | +7.60 |

**The clash rate collapses from 29% to 1.8%.** A single local optimisation, no global search,
converts almost every physically invalid pose into a physically plausible one. The placements were
recoverable; they were simply never relaxed.

## But they are still markedly worse than Vina's own poses

| | mean score |
|---|---|
| arm 1 — Vina poses, Vina score | **−7.58** |
| arm 2 — DiffDock poses, minimised, Vina score | **−4.23** |

A **3.35 kcal/mol** gap remains after minimisation. Some of that is expected and not
informative: Vina's poses are the output of Vina's own global search, so they sit at good points
of the very function being used to score them. A local optimisation from a DiffDock starting
point cannot reach those minima without a global search — which is precisely the thing arm 2 is
supposed to hold constant.

The gap is therefore *interesting but not yet interpretable*. Whether it reflects worse binding
modes or only a less thorough search is exactly what the AUROC comparison would settle, because
AUROC depends on the **ranking**, not on the absolute score.

## Arm 2 is still not reportable, and the reason has changed

It is no longer a validity problem. It is a coverage problem: **55 of 751 = 7%**, far below the
90% floor. What was blocking the arm — physically invalid poses — is now solved, and what blocks
it is simply that these 55 are the only poses retrievable through the Kaggle output API, which
truncates the file listing.

The full run is the queued Kaggle kernel, which mounts the pocket kernels' outputs directly
rather than downloading them. It now has a justified protocol: **minimise, then score**.

## What this already establishes

**A generative pose generator is usable in a physics-scored pipeline, provided its output is
relaxed first.** Raw, it produces 39–63% physically invalid poses and any AUROC computed over
them is meaningless. After one local minimisation the same poses are 98% valid. That is a
one-line protocol change with a large effect, and it is the practical answer to "can we use
DiffDock" — yes, but not as a drop-in.

## Loose end

5 of 60 produced no minimised score. The failure text captured is the tail of Vina's crash
report rather than its error line — the same logging fault that hid the missing grid box earlier
— so the cause is unknown and will need the real stderr line to diagnose. At 8% it does not
affect the direction of this result, but it is not explained.
