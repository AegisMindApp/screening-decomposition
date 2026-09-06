# The increment transfers; the model does not. And the only configuration usable on a novel target is the one we reported as failing.

**7 September 2026.** `transfer_test.py`, results in `TRANSFER.json`. Rules fixed in
`PREREGISTRATION_TRANSFER.md` (`fc6816b51`) before any transfer estimate existed. All three
controls pass.

## Verdict against the pre-registered rule: TRANSFERS

Δ_transfer > 0.0201 with a CI excluding zero, in both directions and under both scaling regimes.

| direction | regime | descriptors alone | + Boltz-2 | Δ_transfer |
|---|---|---|---|---|
| Mpro → FXa | fully cold | 0.4011 | 0.4890 | +0.0877 [+0.0562, +0.1188] |
| Mpro → FXa | recipient-standardised | 0.4013 | 0.4908 | +0.0893 [+0.0575, +0.1209] |
| FXa → Mpro | fully cold | 0.4553 | 0.7036 | +0.2481 [+0.2096, +0.2867] |
| FXa → Mpro | recipient-standardised | 0.4536 | 0.7011 | +0.2474 [+0.2091, +0.2856] |

Controls: each donor model scores its own target at or above its in-target out-of-fold value
(0.8698 vs 0.8647; 0.7873 vs 0.7803) — **PASS**. Permuted recipient labels give Δ = −0.0137,
inside [−0.05, +0.05] — **PASS**.

## Reading the rule as written is not the whole story

The rule asked whether the *increment* survives transfer. It does. But look at the absolute
numbers: **the cold descriptor model is below chance in both directions** (0.4011, 0.4553), and
the cold combination reaches only **0.4890** on Factor Xa — chance. A model at 0.489 is not
usable no matter how much better it is than one at 0.401.

So the honest statement is narrower than the verdict: **the Boltz-2 increment transfers; the
combined model does not.** In one direction (FXa → Mpro) transfer produces a usable 0.7036; in
the other it does not.

## What the failure of descriptor transfer means

Seven descriptors fit on one target and applied to another rank actives **worse than random**.
That is not weak transfer, it is inverted transfer: the property signature separating actives
from inactives on Mpro is close to the opposite of the one on Factor Xa. The descriptor
component is target-specific *bias*, not chemistry that generalises, which is precisely the
reading §3.3 gives it as a bias control. This is the strongest evidence in the project for that
framing, and it arrived from a test designed to ask something else.

## The finding that matters, sorted by what a screener actually has

| configuration | Mpro | FXa | labels required |
|---|---|---|---|
| **Boltz-2 alone, rank by score** | **0.7913** | **0.7227** | **none** |
| descriptors, in-target out-of-fold | 0.7673 | 0.6994 | the target's own actives |
| descriptors + Boltz-2, in-target out-of-fold | 0.8647 | 0.7803 | the target's own actives |
| descriptors, cold from the other target | 0.4553 | 0.4011 | another target's actives |
| descriptors + Boltz-2, cold from the other target | 0.7036 | 0.4890 | another target's actives |

**On a novel target with no known actives, the only usable row is the first.** Boltz-2 alone needs
no labels, no receptor preparation, no docking box, and reaches 0.72–0.79 on both targets.

That configuration is the one §3.4 reports as **not demonstrated** — it fails to beat the
descriptor baseline by the pre-registered +0.04 margin (+0.026 Mpro, +0.019 FXa, both intervals
spanning zero). Both statements are true and they are not in tension: the in-target comparison
asks whether Boltz-2 beats a baseline fitted on the answers, and the prospective question asks
what you can run when there are no answers to fit to.

## Consequence for the manuscript

§3.3 states that its descriptor comparison is supervised and therefore a bias control rather than
a deployable method. **§3.4 makes the same kind of claim without the same caveat** — the +0.093
and +0.075 combination gains are also fitted on each target's own labels. That asymmetry must be
fixed, and this result supplies both the caveat and its resolution:

- the combination gain is real but requires the target's labels, so it is not a screening
  protocol;
- the Boltz-2 score alone requires nothing and is the deployable configuration;
- the descriptor baseline, transferred, is actively harmful.

## Note on the in-target numbers

This script recomputes the in-target deltas as +0.0976 (Mpro) and +0.0809 (FXa) against the
published +0.093 and +0.075. The difference is cross-validation fold assignment, the same effect
documented in `../pose_ensemble/RECONSTRUCTION.md`. The transfer quantities above are not subject
to it: the model is fit on the donor in full and evaluated on a disjoint target, so there is no
fold assignment on the evaluated side.
