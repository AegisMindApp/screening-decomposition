# Stage 2 at ten seeds: the search-effort verdict flips from NEGLIGIBLE to INCONCLUSIVE

24 September 2026. Seeds 7-10 completed overnight, taking both arms to ten. Nothing about the
estimator or the reading rules changed; only the number of replicates did.

Re-derive: `analyse_seeds.py --exh 32 --dir analysis/seed_floor/results`,
`floor_exh4_vs_exh32.py`, `equivalence.py`, `positive_control.py`.

## The floor

| | n = 6 | **n = 10** | change |
|---|---|---|---|
| exh=32 AUROC sd, 343 compounds | 0.00641 | **0.00619** | −3% |
| 95% upper bound (the quotable floor) | 0.01338 | **0.01019** | **−24%** |
| χ² penalty on an sd | 2.09× | 1.65× | |

**Replication buys the bound, not the estimate.** The sd was already settled at six seeds; what
four more seeds bought was a 24% tighter upper bound, because the χ² penalty falls. Since it is
the bound that gets charged against a claim, a floor quoted from few replicates is loose in
exactly the direction that admits effects.

## The verdict that changed

| | n = 6 | **n = 10** |
|---|---|---|
| paired Δ AUROC, exh32 − exh4 | +0.004 | **+0.0082** |
| 95% CI | [−0.006, +0.014] | **[+0.0016, +0.0148]** |
| equivalence margin *d* | 0.01346 | **0.01025** |
| TOST p vs −d / +d | 0.0033 / 0.0328 | 0.0001 / **0.25** |
| verdict | NEGLIGIBLE | **INCONCLUSIVE** |

Two things moved and both went against the claim. The point estimate **doubled**, so +0.004 was
an underestimate rather than scatter around zero — the CI now excludes zero, making this a
detectable effect. And the margin **shrank 24%**, so the bar it had to stay under came down to
meet it. An equivalence verdict depends on the margin as much as on the effect; here more data
moved both.

The n = 6 write-up already called the verdict "close to its boundary … rather than a comfortable
null", with a 90% upper bound of +0.0121 against a margin of 0.01346. That caveat was right and
it was not enough: the verdict still had to be withdrawn.

## What it costs

**NEGLIGIBLE now has no instance among the five interventions.** The category the framework
exists to distinguish is exercised only by the planted controls. That is a weaker paper than one
reporting a real measured-negligible, and it is what the evidence says.

The control anticipated it. Tier 2 of `positive_control.py` — a planted effect carried on a
*different* seed's run, so real reseed noise is in the interval — returns a planted +0.004 as
**INCONCLUSIVE**, not NEGLIGIBLE. The control was telling us that an effect of that size is not
resolvable under reseeding before the real measurement said the same thing.

Re-run against the tightened margin, all five planted effects still recover with the verdicts
fixed in advance (+0.000, +0.004 negligible; +0.050, +0.140, −0.050 survive), so the rule itself
is unchanged.

## Unchanged

The two floors remain statistically indistinguishable: matched on 339 compounds and all ten
shared seeds, exh=4 gives 0.00844 and exh=32 gives 0.00618, ratio 0.73, Pitman-Morgan **p = 0.38**
(it was 0.62 at n = 6 — the p moved without resolving). The common assumption that less search is
more stochastic still has no support here.

The within-target panel-size check moves slightly: the 755-compound sd of 0.00497 predicts 0.00742
on 339 and the matched re-measurement over ten seeds gives **0.00844**, a ratio of **1.14** rather
than the 1.09 recorded at six seeds. Still within the ±25% band the cross-target test used.
