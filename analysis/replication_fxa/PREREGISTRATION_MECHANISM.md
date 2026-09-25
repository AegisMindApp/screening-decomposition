# Pre-registration: is gnina a constant-quality ranker, or does it track the target?

Written 25 September 2026, **before the analysis is run**, using only data already committed.
Tests the one hypothesis in `resolution_v4_draft.md` §3.6 that is labelled as unproven: that
gnina's **+0.140 on Mpro measured recovery from a broken baseline**, not a property of the
scoring function.

## The obvious test is circular, which is why the design is not the obvious one

The tempting analysis is to split the panels into sub-panels, regress the gnina **delta** on the
Vina baseline, and report a negative slope. That is guaranteed to succeed and means nothing:
for *any* score uncorrelated with Vina, delta = (its AUROC) − (Vina's AUROC), so the slope
against the Vina baseline is **−1 by construction**. A permuted-gnina control would "confirm" the
hypothesis perfectly. Selecting on the baseline and then regressing on it is the same error as
deriving a threshold from the data it judges.

## The non-circular formulation

Regress **gnina's absolute AUROC** on **Vina's absolute AUROC** across sub-panels. The slope
distinguishes the two explanations directly:

| slope | meaning |
|---|---|
| ≈ **0** | gnina ranks at a fixed quality regardless of whether the target is easy for Vina. Its apparent gain is whatever gap that leaves. **Supports** the broken-baseline reading. |
| ≈ **1** | gnina tracks target difficulty the way Vina does; the two succeed and fail together, and the Mpro gain is a property of the scoring function on that target. **Refutes** it. |

The delta-slope is just this slope minus one, so nothing is lost — but the quantity being tested
is no longer guaranteed by construction.

## Sub-panels

Disjoint groups formed from **Butina scaffold clusters** (Morgan radius 2, 2048 bits, Tanimoto
0.65 — the clustering already used in §3.7), aggregated into bins of at least 80 compounds so
each sub-panel's AUROC is estimable. Clusters are assigned to bins by chemistry, never by score
or label. Both targets, analysed separately and pooled.

## The control that decides whether the test can fire at all

Sub-panel AUROCs are noisy, and **noise in the x-variable biases a slope toward zero** — that is,
toward the hypothesis. A slope near 0 is therefore only meaningful if the method can recover a
slope near 1 when one truly exists.

**Positive control: Vina at a different seed, regressed on Vina.** Two runs of the same protocol
track target difficulty perfectly, so the true slope is 1. If this control returns a slope
materially below 1, regression dilution is eating the signal and **the test is void** — the gnina
result is not read at all.

A permuted-gnina negative control is also reported: its slope must be ≈ 0 with its absolute AUROC
≈ 0.5 everywhere, confirming that a genuinely uninformative score looks like the "constant
quality" case.

## Reading rules, fixed now

Let *b* be the gnina-on-Vina slope with a bootstrap 95% CI over sub-panels, and *b_ctrl* the
positive control's slope.

- **VOID** — *b_ctrl* < 0.70. The design cannot detect tracking; nothing is read.
- **SUPPORTS the broken-baseline reading** — *b_ctrl* ≥ 0.70 **and** *b*'s CI excludes 0.5 from
  below (i.e. gnina tracks the target much less than a re-run of Vina does).
- **REFUTES it** — *b*'s CI excludes 0.5 from above.
- **INCONCLUSIVE** — otherwise. The hypothesis is then **cut from the manuscript**, not softened:
  an explanation that cannot be tested on the data in hand does not belong in the paper as
  anything more than an absence.

0.5 is the midpoint between the two competing predictions and is chosen before seeing any slope.

## What no outcome here establishes

Not that gnina is a bad scoring function — it is being run `--score_only` on Vina's poses, which
is not what it is designed for. Not anything about a third target. And a SUPPORTS verdict explains
the sign reversal without proving it: the mechanism would be consistent with the data, not
demonstrated to be its cause.
