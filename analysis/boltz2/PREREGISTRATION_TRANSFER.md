# Pre-registration: does the descriptor + Boltz-2 combination transfer between targets?

**Written 7 September 2026, before any transfer estimate is computed.**

## The gap this closes

§3.4 reports that descriptors + Boltz-2 beats descriptors alone by **+0.093** (Mpro) and
**+0.075** (Factor Xa). Both are **out-of-fold on each target's own actives and inactives** —
the combination weights are fit using the labels of the very target being screened.

§3.3 already states this limitation for the descriptor baseline and frames that comparison as a
bias control rather than a deployable method. **§3.4 does not carry the same caveat although it
has the same property.** On a novel target with no known actives, there are no labels to fit the
combination with.

This tests whether the gain survives when the weights come from a *different* target.

## Design

Two directions, each fully specified:

- **Mpro → Factor Xa**: fit on all 750 Mpro compounds, apply cold to all 885 Factor Xa compounds.
- **Factor Xa → Mpro**: the reverse.

Two models per direction, identically fitted:

- **A** = descriptors only (7 features)
- **B** = descriptors + Boltz-2 `affinity_probability_binary` (8 features)

Logistic regression, as in §3.4. No fold splitting is needed or used: the model is fit on the
donor target in full and evaluated on a disjoint recipient target, so there is no fold noise on
the evaluated side. Confidence intervals are a bootstrap over recipient compounds, 10,000
resamples.

**Two scaling regimes, both reported:**

1. **Fully cold** — the donor's feature scaler is applied to the recipient. Nothing about the
   recipient is used.
2. **Recipient-standardised** — features are standardised using the recipient's own feature
   distribution. This uses recipient *data* but no recipient *labels*, so it remains prospective
   and is the realistic deployment setting.

## Reading rules — fixed now

Primary quantity: **Δ_transfer = AUROC(B) − AUROC(A)** on the recipient target.

- **TRANSFERS** — Δ_transfer > **0.0201** (the scoring floor) with a bootstrap 95% CI excluding
  zero, in **both** directions. The combination is a protocol that could be applied to a novel
  target, and §3.4 stands as a method claim.
- **DOES NOT TRANSFER** — otherwise. The +0.093 and +0.075 are per-target refits, §3.4 must carry
  the same supervision caveat as §3.3, and the contribution is a bias-control measurement rather
  than a screening method.

Reported alongside: the recipient's own in-target out-of-fold Δ (+0.093 / +0.075), so the loss
from transferring is visible; and the absolute AUROCs, so it is clear whether the cold model is
usable at all rather than merely better than another cold model.

## Controls

1. **Sanity** — the donor-fitted model applied back to the *donor* must reproduce that target's
   in-target AUROC to within 0.02. A model that cannot score its own training target is broken.
2. **Label permutation** — with recipient labels shuffled, Δ_transfer must fall in [−0.05, +0.05].

## Pre-committed disclosure

Reported either way. A negative result here weakens a claim we have already published in a
preprint and would require correcting §3.4 before journal submission; that is the point of
running it now rather than after review.
