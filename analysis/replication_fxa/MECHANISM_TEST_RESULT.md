# Testing the quotable explanation: INCONCLUSIVE, so it is cut

25 September 2026. Answers `PREREGISTRATION_MECHANISM.md`, written before the analysis ran.
Re-derive: `python analysis/replication_fxa/mechanism_test.py` → `MECHANISM_TEST.json`.

## What was tested

`resolution_v4_draft.md` §3.6 carried one explicitly unproven idea: that gnina's **+0.140 on
Mpro measured recovery from a broken baseline** rather than a property of the scoring function.
It was the paper's most quotable sentence and its least supported.

The obvious test — regress the gnina *delta* on the Vina baseline across sub-panels — is
worthless: for any score uncorrelated with Vina the slope is **−1 by construction**, so a
permuted-gnina control would "confirm" the hypothesis perfectly. The pre-registration therefore
regressed gnina's **absolute** AUROC on Vina's, where slope 0 means constant-quality ranking
(supports) and slope 1 means tracking target difficulty (refutes).

## Result

19 scaffold-binned sub-panels of 80–95 compounds across both targets, Vina AUROC spanning
**0.241 to 0.854**.

| regression | slope | 95% CI | reading |
|---|---|---|---|
| **positive control** — Vina (seed 2) on Vina | **+0.958** | [+0.918, +1.006] | **PASS**: the design detects tracking |
| negative control — permuted gnina on Vina | +0.272 | [+0.034, +0.600] | mean AUROC 0.502, as expected |
| **gnina on Vina** | **+0.577** | **[+0.173, +0.956]** | **straddles 0.5** |

**INCONCLUSIVE.** gnina sits between the two competing predictions and 19 sub-panels cannot say
which. Per the pre-registered rule the hypothesis is **cut from the manuscript**, not softened.

## Why this is a real verdict and not a broken instrument

The positive control is the part that makes the null readable. Sub-panel AUROCs are noisy, and
noise in the x-variable biases a slope **toward zero** — that is, toward the hypothesis. Had the
control come back at 0.6, a gnina slope of 0.577 would have meant nothing. It came back at
**0.958**, so the design recovers a true slope of 1 essentially intact, and gnina's 0.577 is a
measurement rather than an artefact of dilution.

The negative control also behaved: a permuted gnina score ranks at chance (0.502) on every
sub-panel and its slope is near zero, confirming that a genuinely uninformative score looks like
the "constant quality" case.

So the honest position is that gnina partially tracks target difficulty — nominally less than a
Vina re-run does — and the data do not resolve whether it tracks enough to call the Mpro gain a
property of the scoring function. **We do not now switch to comparing gnina's slope against the
control's rather than against 0.5.** The 0.5 threshold was fixed in advance for exactly this
situation, and moving it after seeing a slope of 0.577 would be the error this project keeps
catching in its own work.

## What is cut, and what stays

**Cut:** the broken-baseline explanation, in all its forms.

**Stays**, because all of it is measured:

- gnina rescoring gains **+0.1397** on Mpro and loses **0.0583** on Factor Xa, both with intervals
  clear of their targets' floors.
- On Mpro, Vina ranks below chance (0.399) and gnina reaches 0.539; on Factor Xa, Vina reaches
  0.687 and gnina 0.629. Those are the absolute numbers, and a reader may draw their own
  inference from them.
- The reversal is therefore reported as a **fact without an established mechanism**.

That is a smaller claim than the paper had this morning and a better-supported one. An
explanation that cannot be tested on the data in hand does not belong in a paper as anything more
than an absence — and the cost of finding out was one afternoon of analysis on data already
committed.
