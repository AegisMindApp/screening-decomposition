# Stage 2: the eight-fold search-effort effect does not survive seed replication

20 Sep 2026. Estimator fixed in `paired_exhaustiveness.py` before any stage-2 data existed.
Both the t-test and Wilcoxon are reported because the script computes both unconditionally and
choosing between them after seeing which is friendlier is the circularity this whole exercise
exists to remove.

**Every figure below names its panel.** A floor is a property of (method, protocol, panel size),
and the JCIM rejection turned on a threshold imported from a different measurement context.

## The effect, paired within seed

Exhaustiveness 32 minus exhaustiveness 4, same seed, **339 compounds common to all arms**:

| seed | exh=4 | exh=32 | delta |
|---|---|---|---|
| 1 | 0.4562 | 0.4504 | −0.0058 |
| 2 | 0.4440 | 0.4588 | +0.0148 |
| 3 | 0.4358 | 0.4408 | +0.0050 |
| 4 | 0.4469 | 0.4499 | +0.0029 |
| 5 | 0.4553 | 0.4480 | −0.0073 |
| 6 | 0.4405 | 0.4558 | +0.0153 |

**paired delta = +0.00417, sd 0.00968, n=6, same sign 4/6**
**95% CI [−0.00598, +0.01433] — includes zero. t-test p=0.339, Wilcoxon p=0.5625.**

The manuscript reports **+0.0185** from a single unpaired exh=32 run minus a single unpaired
exh=4 run. That value lies **outside the upper bound of this CI**, so the discrepancy is not
merely that the original estimate was noisy: re-measured properly it is roughly a third of the
size and indistinguishable from nothing.

*Caveat on that comparison.* The manuscript's +0.0185 was computed on its own panel, not on
these 339 compounds. Comparing a point estimate across panels is the very substitution this
document warns about, so treat "outside the CI" as indicative of magnitude, not as a test.
What is measured here without qualification is that **on a fixed panel with seeds held in
common, the effect is +0.004 with a CI spanning zero.**

Pairing was worth doing, and that was measured rather than assumed: the paired sd of the
difference is 0.00968 against an unpaired expectation of sqrt(2) x sd = 0.01148.

## The exh=32 floor, measured at last

`method_ladder.json` carried `vina_exh32` as a `MEASURING` placeholder, which meant the ladder
had nothing to charge the manuscript's actual protocol against.

| | value | panel |
|---|---|---|
| AUROC sd across 6 seeds | **0.00641** | 343 compounds |
| 95% upper bound (chi-squared, one-sided) | **0.01338** | 343 compounds |
| per-ligand score sd, median | 0.0346 eV | |
| rank span, median | 10 | |

For contrast, exh=4 on **755** compounds gave sd 0.00497 and a 0.00817 bound.

## Does the exh=4 floor bound the exh=32 floor?

This is what `docs/papers/RESUBMISSION_GATE.md` Q6 blocks on. Matched on the same 339 compounds
and the same 6 seeds (`floor_exh4_vs_exh32.py`):

| | sd | 95% upper bound |
|---|---|---|
| exh=4 | 0.00812 | 0.01696 |
| exh=32 | 0.00626 | 0.01307 |

Ratio 0.77. Pitman–Morgan **paired** variance test (paired, not F, because the arms share seeds
and their AUROCs are correlated): r = −0.256, **p = 0.624**.

So the substitution is **conservative in direction but not established**. The honest resolution
is not to argue about whether it is safe — it is to use the exh=32 floor, which now exists.

## Floors are panel-size specific, and that is measured

Stage 1 put the exh=4 sd at 0.00497 on 755 compounds. If AUROC variance scales as 1/n, that
predicts 0.00742 on a 339-compound panel. The matched re-measurement gives **0.00812** — a
ratio of **1.094**, from a prediction using data the check does not otherwise touch.

Consequence: the 0.00817 bound measured on 755 compounds is about **1.5x too small** to charge
against a claim computed on 339, and using it there would wave through effects that are inside
the noise. Re-scale by sqrt(n_measured / n_claim), or measure the floor on the claim's own panel.

## An apparent paradox that is not one

Per-ligand reproducibility improves sharply with search effort — score sd **0.096 → 0.035 eV**,
median rank span **70 → 10** — while the AUROC sd barely moves. Both are true because AUROC sd
is set by panel size, not by per-ligand precision. An earlier reading of mine, that exh=32 was
*noisier*, came from comparing a 343-compound sd against a 755-compound one and does not
survive matching.

## What this licenses for the resubmission

- The exhaustiveness claim as stated does not survive. It should be reported as **+0.004
  [−0.006, +0.014], n=6 seeds, 339 compounds**, not as +0.0185.
- Q6 of the gate can now be answered with a measured exh=32 floor instead of a substituted one.
- Nothing here rescues the original number, and nothing here depends on rescuing it.
