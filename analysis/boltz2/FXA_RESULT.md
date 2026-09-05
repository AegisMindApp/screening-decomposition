# Factor Xa: Boltz-2 replicates the Mpro result on an independent target

**5 September 2026.** 886 dispatched, 885 scored, 1 failure. 27.8 GPU-hours across six free
Kaggle P100 sessions — no paid compute. Rules fixed in `PREREGISTRATION_TARGETS.md`
(`68bd9b94a`) before any shard ran: descriptor baseline 0.7022, bar 0.7422.

## Controls — both pass

| control | result | |
|---|---|---|
| **Positive** — reproduce the pre-registered baseline | **0.7041** vs 0.7022 ± 0.005 | **PASS** |
| **Permutation** — shuffled labels | 0.5165, inside [0.45, 0.55] | **PASS** |
| **Dropout** — 1 inactive of 886 | bounded 0.7212–0.7232 vs 0.7227 | **immaterial** |

The positive control is the one that matters: the baseline was written down *before* the first
Kaggle session and reproduced to within 0.002 afterwards. A pipeline that reproduces a number
we already knew is one whose unknown numbers are worth reading. The dropout rule is degenerate
at n = 1 for the same reason as Mpro — a single failure gives 0% or 100% by construction — so
the impact was bounded instead: imputing the dropped inactive at best and worst possible scores
moves AUROC by **0.002**.

## Result

| | AUROC |
|---|---|
| Vina (repaired receptor) | 0.6775 |
| seven free descriptors | 0.7041 |
| **Boltz-2** | **0.7227** |
| **descriptors + Boltz-2** | **0.7791** |

| rule | Factor Xa | Mpro (arm 1) |
|---|---|---|
| Boltz-2 − descriptors | +0.0186 [−0.030, +0.066] **not demonstrated** | +0.0259 [−0.023, +0.074] **not demonstrated** |
| **desc + Boltz-2 − descriptors** | **+0.0751 [+0.051, +0.099] SUPPORTED** | **+0.0934 [+0.062, +0.125] SUPPORTED** |

## What replicated, and what did not

**The combination result replicated.** Two unrelated proteins — a viral protease and a
coagulation factor — two independently assembled compound panels, different active rates (45.8%
vs 34.3%), six separate Kaggle sessions on different hardware to the Mpro run, and the
combination gain lands at +0.075 against +0.093 with overlapping intervals.

**The standalone null also replicated.** Boltz-2 alone fails the bar on both targets. Its
point estimate wandered +0.036 → +0.045 → +0.019 as the panel filled from 33% to 67% to 100% —
motion around a null, which is exactly what an interval spanning zero predicts and a useful
reminder not to read a drifting partial estimate as a trend.

So the claim this work supports is narrow and specific: **Boltz-2 carries information that
seven free descriptors do not, but is not a substitute for them.** On the evidence here the
useful configuration is both together, not either alone.

## Cost

27.8 GPU-hours, free. Sharded six ways for Kaggle's 12-hour session cap, interleaved by sorted
name so every shard held the panel's label mix (43.9%–48.9% against 45.7%) — which also made
the partial results legitimately analysable rather than order-biased.
