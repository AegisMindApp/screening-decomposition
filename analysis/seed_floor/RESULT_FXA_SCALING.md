# The floor's scaling constant TRANSFERS across targets — and the test that says so barely discriminates

Answers `analysis/seed_floor/PREREGISTRATION_FXA.md`, completed 23 September 2026. Ten seeds,
Factor Xa (ChEMBL244), 886 compounds, AutoDock Vina exhaustiveness 4, varying only `--seed`.
Zero failures; the cross-seed intersection is the full 886, well clear of the 700 floor for
INCONCLUSIVE.

Re-derive: `python analysis/seed_floor/analyse_seeds.py --exh 4 --dir analysis/seed_floor/results_fxa`
(output `fxa_floor_n10.json`).

## Verdict: TRANSFERS

| | value |
|---|---|
| Predicted sd, fixed before the run | 0.00497 × √(755/886) = **0.00459** |
| Pre-registered band (±25%) | [0.00344, 0.00573] |
| **Measured sd (n = 10)** | **0.004272** |
| Error against prediction | **6.9%** |
| 95% one-sided upper bound (the quotable floor) | 0.00703 |

6.9% against the 9% the within-target check achieved on Mpro. The 1/n scaling is a property of
the AUROC statistic, not of Mpro, so a floor can be stated for a panel of size *n* without
running anything on that panel.

## Two things that would be easy to leave out

**The verdict turned on the tenth seed.** Recomputing the sd as seeds accumulate — holding the
panel fixed at the n = 10 intersection, so this is not identical to the live numbers reported at
the time, which used their own smaller intersections:

| n | 4 | 5 | 6 | 7 | 8 | 9 | **10** |
|---|---|---|---|---|---|---|---|
| sd | 0.003906 | 0.003445 | 0.003720 | 0.003397 | 0.003441 | 0.003436 | **0.004272** |

At n = 9 the value sat at 0.003436 against a lower band edge of 0.003440 — **outside by
4×10⁻⁶**, which is a DOES NOT TRANSFER verdict on a margin of one part in a thousand. Seed 10
returned the lowest AUROC of the ten (0.66322) and moved the sd 24%. The pre-registration fixed
n = 10 in advance, so the reading is the pre-registered one and not a choice made after seeing
the trajectory. But nobody should take an n = 10 sd as a stable quantity.

**The ±25% band is tighter than the estimator's own precision.** The two-sided 95% CI on an sd
from 10 replicates is [0.00294, 0.00780] — a range of [0.69×, 1.83×] the point estimate, against
a band of [0.75×, 1.25×]. The pre-registration justified ±25% by citing a χ² penalty of "about
1.24×"; the actual one-sided upper penalty here is **1.65×**, and the two-sided interval is wider
still. So the document made, in its own justification, the error it named and was trying to
avoid.

This does **not** void the verdict, because the consequence is asymmetric and runs in the
favourable direction:

- **Passing a band narrower than your measurement precision is strong evidence.** The measured
  value had ample room to land outside by sampling noise alone, and did not.
- **Failing it would have been weak evidence**, and largely uninterpretable — a DOES NOT TRANSFER
  could have come purely from the sd estimator's own scatter, as the n = 9 knife-edge shows.

Stated here rather than inferred later, per the standing rule that the pass/null asymmetry belongs
in the verdict function and not in the discussion.

## What this licenses, and what it does not

**Licensed.** The manuscript may state the floor as a function of panel size for AUROC on a
docking benchmark, rather than as an observation about Mpro. Two targets, 755 and 886 compounds,
both within 9% of the 1/n prediction.

**Not licensed.** Nothing here says Factor Xa's *interventions* behave like Mpro's, and the
pre-registration says so explicitly. A second panel makes the floor general; it does not make the
intervention results general. Two targets is also two, not many — the constant is now measured on
a pair of serine-protease-and-cysteine-protease panels of similar size, and a panel an order of
magnitude smaller or larger remains an extrapolation.

---

## Addendum, 25 September 2026 — the knife-edge settles, at n = 12

The verdict above turned on the tenth seed: at n = 9 the sd sat **4×10⁻⁶ outside** the lower band
edge, which would have read DOES NOT TRANSFER on a margin of one part in a thousand. Two further
seeds were run to find out whether the published reading was luck. Same 886-compound panel,
verified identical across all twelve (zero symmetric difference, nothing lost from the
intersection).

| n | 8 | 9 | **10 (published)** | 11 | **12** |
|---|---|---|---|---|---|
| sd | 0.003441 | 0.003436 | **0.004272** | 0.004477 | **0.004401** |
| two-sided 95% CI | [0.00228, 0.00700] | [0.00232, 0.00658] | [0.00294, 0.00780] | [0.00313, 0.00786] | **[0.00312, 0.00747]** |

**The pre-registered test is not re-read.** It fixed n = 10 and a prediction of 0.00459 before any
Factor Xa docking ran; it measured 0.00427, ratio 0.93, TRANSFERS. That stands as published.
Re-reading it against a later, larger sample would retro-fit the thing that made it worth running.

**As robustness only:** at n = 12 the sd is **0.00440**, ratio **0.96** to the same prediction —
still inside the ±25% band, and now **21% of the prediction clear of the lower edge** rather than
sitting on it. Against the newer fifteen-seed Mpro prediction (0.00407) the ratio is 1.08, also
inside. The conclusion holds under every combination of *n* tried on either target.

What settled is the *stability*, not the verdict. The n = 9 reading was a genuine knife-edge and
the caveat above was right to flag it; three more seeds show the estimator had simply not
converged, and the low values at n = 8–9 were the outliers rather than n = 10 being a lucky draw.

The quotable floor also tightens: **0.00682 at n = 12** against 0.00703 at n = 10, as the χ²
penalty falls from 1.645× to 1.551×.
