# What is docking actually worth? +0.011 AUROC over seven free descriptors

**2 September 2026.** Run to answer whether building a derivative docking algorithm is a
sensible use of our time. Five-fold CV, out-of-fold predictions, same compounds and folds
for every arm.

## The measurement

| target | Vina alone | 7 descriptors | descriptors + Vina | **docking adds** |
|---|---|---|---|---|
| **Mpro** (751, 257 act) | 0.4264 | 0.7565 | 0.7564 | **−0.0002** |
| **Factor Xa** (886, 405 act) | 0.6645 | 0.6969 | 0.7083 | **+0.0114** |

The descriptors are molecular weight, clogP, H-bond donors, H-bond acceptors, rotatable
bonds, TPSA and formal charge. They cost microseconds. Docking costs ~20 minutes per
compound at exhaustiveness 32.

**On Mpro, docking contributes nothing** — not a little, nothing to four decimal places.
**On Factor Xa, the target where docking works and passed every gate, it contributes
+0.011 AUROC.**

One more number worth having: Vina *as a fitted feature* on Mpro reaches 0.5692 against
0.4264 raw. So the score does carry information — the sign and scale are simply wrong — and
that information turns out to be **entirely redundant with the descriptors**. Which is what
the AcrB/CTSS result already implied: 87% of a docking score is not target-specific.

## So: could we build a better docking algorithm than what is published?

**No, and the binding constraint is not algorithm design.**

**1. The prize is +0.011 AUROC.** That is the measured marginal value of the entire docking
computation over descriptors that are free, on our best-behaved target. A better docking
algorithm competes to improve *that*, not to improve 0.7083.

**2. We cannot tell whether we succeeded.** Published methods are evaluated on LIT-PCBA (15
targets), DUD-E (102), PoseBusters. We have **one** target properly gated (Factor Xa), one
below chance (Mpro), one property-confounded (PD-L1). A method that improved our Factor Xa
number by 0.02 would be indistinguishable from noise, and we would have no way to know it
generalised.

**3. No wet lab.** Every advance in this space this year — Mythos's binders at ~50% hit rate,
ArtiDock's pose accuracy — was validated experimentally. A retrospective AUROC on one target
is not a competitive claim.

**4. Eight GPU hours a week.** Against groups with hundreds of funded full-time researchers.

## What is actually within reach, and is what this week produced

Not an algorithm — a **decomposition**. On one benchmark, with the reading rule committed
before each number existed, changing exactly one thing at a time:

| intervention | ΔAUROC |
|---|---|
| scoring function (identical poses) | **+0.140** |
| receptor protonation | +0.045 Mpro / +0.013 FXa / ~0 AcrB |
| 8× search effort | **−0.019** (top-10 unchanged) |
| **docking over free descriptors** | **−0.000 Mpro / +0.011 FXa** |

That is rare, and not because it is hard. It is because holding everything else constant is
slow and produces unflattering numbers. Most published comparisons change the method, the
benchmark and the preparation together.

The honest positioning: **we are not competitive at building docking methods, and we are
unusually well set up to measure them.** The receptor-donor defect, the size-selected dropout
that would have biased arm 2, the exhaustiveness null — each was found by a control that
could have failed, not by a better algorithm.

## The one caveat that keeps this from being a general claim

Two targets. Factor Xa's +0.011 is a single estimate with no confidence interval computed
here, and Mpro is a target where our own protocol scores below chance. This says what docking
is worth *to us, on these benchmarks*. It does not say what docking is worth in general, and a
group with a properly powered benchmark suite could reasonably measure something different.
