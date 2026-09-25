# Two objections closed from committed data: the metric, and the AUROC level

Answers `PREREGISTRATION_METRIC_SCOPE.md`, 24 September 2026. No new compute; both use the ten
exhaustiveness-4 and ten exhaustiveness-32 seeds already in the repository.
Re-derive: `python analysis/seed_floor/metric_scope.py` → `METRIC_SCOPE.json`.

## 1. The floor is strongly metric-specific, and the metric the field reports is the bad one

| exh | metric | mean | sd | 95% bound | CV = sd/mean |
|---|---|---|---|---|---|
| 4 | AUROC | 0.4175 | 0.00497 | 0.00817 | **0.0119** |
| 4 | **EF@1%** | 2.2768 | 0.23225 | **0.38210** | **0.1020** |
| 4 | EF@5% | 1.3143 | 0.10308 | 0.16958 | 0.0784 |
| 4 | BEDROC(α=20) | 0.4489 | 0.01671 | 0.02749 | 0.0372 |
| 32 | AUROC | 0.4531 | 0.00619 | 0.01019 | 0.0137 |
| 32 | **EF@1%** | 2.5751 | 0.54287 | **0.89313** | **0.2108** |
| 32 | EF@5% | 1.4905 | 0.11496 | 0.18913 | 0.0771 |
| 32 | BEDROC(α=20) | 0.4909 | 0.01509 | 0.02482 | 0.0307 |

**Pre-registered prediction: CV(EF@1%) > 3 × CV(AUROC). Measured 8.57× at exh 4 — CONFIRMED**,
and 15.4× at exh 32.

Read the EF@1% rows directly. At exhaustiveness 4 the enrichment is **2.28 with a resolution
floor of 0.38** — reseeding the identical protocol moves EF@1% by up to a sixth of its own value.
At exhaustiveness 32, on the smaller panel, the floor is **0.89 against a mean of 2.58**: better
than a third. An EF@1% improvement has to be very large before a benchmark of this size can tell
it from the same protocol run twice.

This matters more than any result in the manuscript's AUROC sections, because **EF@1% is what
screening papers report**. AUROC weights the whole ranking; screening cares about the top of it.
A field convention of quoting EF@1% from a single run is quoting a number whose own
irreproducibility is a sixth of its magnitude.

The mechanism is not subtle: EF@1% on 755 compounds is decided by the **top 8 compounds**, and
which 8 they are changes when the seed changes. BEDROC(α=20), designed to weight early
recognition smoothly rather than at a hard cutoff, sits between the two at CV 0.037 — about 3×
AUROC's and a third of EF@1%'s. If a paper must report early recognition, BEDROC is the more
resolvable choice, and that is a recommendation this measurement supports.

## 2. The floor is only weakly level-specific, so it does transport

Vina is **below chance on this panel** — mean AUROC 0.4175 at exhaustiveness 4. A referee whose
benchmark reaches 0.75 is entitled to ask whether a floor measured where the method has no signal
applies to them.

Everything was held fixed except the level: the same 755 compounds, the same labels, the same
class balance and the **same score multiset**, reassigned across compounds so the panel attains
a target AUROC, then perturbed by **resampling each compound's own measured seed deviations**
rather than an assumed Gaussian.

**Positive control.** On the real panel the simulation gives sd **0.00487** against the
**0.00497** measured across ten real re-docking runs — a ratio of **0.98**. The noise model
reproduces the thing it is standing in for, so the trend below is usable. Without this the
simulation would be an assumption dressed as a measurement.

| realised AUROC | induced sd | 95% bound | relative to peak |
|---|---|---|---|
| 0.2785 | 0.00473 | 0.00507 | 1.00 |
| 0.3821 | 0.00466 | 0.00500 | 0.98 |
| **0.5158** | **0.00474** | 0.00508 | **1.00 (peak)** |
| 0.5738 | 0.00447 | 0.00480 | 0.94 |
| 0.6870 | 0.00407 | 0.00437 | 0.86 |
| 0.7939 | 0.00423 | 0.00453 | 0.89 |
| 0.9021 | 0.00419 | 0.00449 | 0.88 |

**Pre-registered prediction: sd peaks near 0.5 and falls away — CONFIRMED.** The peak is at
realised AUROC 0.516, 0.016 from chance.

But the size of the effect is the useful part, and it is **small**: across the entire range from
0.28 to 0.90 the floor varies by about **12%**. So the objection is answered in the direction
that helps. A floor measured on a below-chance panel is conservative for a better-performing one,
and conservative by roughly a tenth, not by a factor. Panel size, at up to 1.5× in this work,
dominates the level term by an order of magnitude.

## What this does not establish

The level sweep is a **simulation**, validated against one real point. Its noise is measured and
its positive control passes, but it has not been checked by re-docking a panel that genuinely
achieves AUROC 0.8 — that would need a target where Vina works, which this panel is not.

Neither result says anything about whether the docking protocol is any *good*. It is not: 0.4175
is below chance and far below the 0.7136 that seven physicochemical descriptors reach on the same
compounds. A resolution floor is a statement about reproducibility, not accuracy, and the two are
independent — a perfectly reproducible method can be worthless, which is close to the case here.
