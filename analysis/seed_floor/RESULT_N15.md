# The stage-1 floor at fifteen seeds: the bound tightens 21%, and the cross-target verdict is robust to it

24 September 2026. Five further seeds of the unchanged exhaustiveness-4 protocol, taking stage 1
from ten to fifteen. Re-derive: `analyse_seeds.py --exh 4 --dir analysis/seed_floor/results`
→ `exh4_floor_n15.json`.

## The panel did not change, and that was checked rather than assumed

The new seeds record `panel: 771` where the originals record `panel: None`, which looks like a
different benchmark. It is not: all five score **exactly the same 764 compounds** as seed 10
(symmetric difference zero), the field was simply added to the runner later. The intersection
across all fifteen seeds is **755 compounds — identical to the published ten-seed panel**, with
nothing lost.

That check matters here. This repository has already been bitten by a superseded 771-compound
run being merged into a 350-compound seed, and by a Factor Xa seed overwriting an Mpro one
because the filename carried no target. A count that looked wrong was the visible symptom both
times.

## The measurement

| | n = 10 (published) | **n = 15** |
|---|---|---|
| AUROC sd | 0.00497 | **0.00441** |
| χ² penalty | 1.645× | **1.460×** |
| 95% upper bound — the quotable floor | 0.00817 | **0.00644** |

| n | 6 | 8 | 10 | 12 | **15** |
|---|---|---|---|---|---|
| sd | 0.00327 | 0.00526 | 0.00497 | 0.00485 | **0.00441** |
| bound | 0.00682 | 0.00946 | 0.00817 | 0.00751 | **0.00644** |

The pattern is the same one the exhaustiveness-32 arm showed going from six seeds to ten, and it
is the paper's own point: **replication buys the bound, not the estimate**. The sd wanders by
±20% up to n≈10 and is settling by 15; the bound falls monotonically after n=8 because the χ²
penalty shrinks. Since the bound is what gets charged against a claim, a floor from few
replicates is loose in the direction that admits effects.

## The cross-target verdict does not depend on which n is used

`PREREGISTRATION_FXA.md` fixed its prediction from the ten-seed value: 0.00497 × √(755/886) =
**0.00459**, band ±25%, measured 0.00427, ratio 0.93 — TRANSFERS. That test stands exactly as
pre-registered and is **not** re-read here; retro-fitting a prediction to a later measurement
would destroy the thing that made it worth running.

As a robustness note only: the fifteen-seed value would have predicted 0.00407, against the same
measured 0.00427 — ratio **1.05**, comfortably inside the same band. So the transfer conclusion
survives either choice of *n*, which is worth knowing precisely because the prediction's input
has now moved.

## What this does not license

It does not license re-reading the within-target panel-size check. That comparison is **matched
on shared seeds** between the exhaustiveness-4 and exhaustiveness-32 arms, and the exh-32 arm has
ten. Substituting a fifteen-seed exh-4 sd into a ten-seed matched comparison would mix
estimators — the kind of substitution this paper is about. The matched figures stay as they are
until both arms have the same seeds.
