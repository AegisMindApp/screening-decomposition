# Pre-registration: does a third scoring approach clear the docking residual ceiling?

**Written 8 September 2026, before any FlashBind score exists.**

## The question, taken verbatim from the paper

§3.4 of `docs/papers/screening_decomposition_preprint.md`: *"on these two targets the ceiling
appears to be a property of the scoring approach rather than of the task. A third scoring
approach clearing it would be the test that generalises this."*

Six docking scores residualise to **0.494–0.573**. Boltz-2 residualises to **0.6567
[0.6439, 0.6648]**. One point outside a band is not a generalisation.

FlashBind (bioRxiv 10.64898/2025.12.22.695983, MIT, `AIDD-Lab/FlashBind`) is the sharpest
available third case because it is *architecturally between the two families*: a fast docking
model for placement, a learned EGNN head for scoring. It claims early enrichment competitive
with Boltz-2 at ~50× lower cost. Surfaced by `analysis/recency_monitor` on its first production
scan.

## Reading rules — fixed now

Primary quantity: **residual AUROC** on Mpro — FlashBind's score with the seven descriptors
regressed out, out-of-fold gradient-boosted, mean over 12 fold seeds with the full range, exactly
as in `analysis/docking_value/RESIDUAL_BAND.md`.

- **CEILING IS THE SCORING APPROACH (claim generalises)** — residual > **0.573**, the top of the
  docking band, with the 12-seed range clear of it. Two non-docking approaches beat the ceiling
  and it is a property of the approach.
- **ENRICHMENT AND ORTHOGONAL SIGNAL COME APART** — residual within **[0.494, 0.573]** while raw
  AUROC is competitive with Boltz-2's 0.7913. A method can rank actives well and still carry no
  more target-specific signal than docking. This would *weaken* our §3.4 claim and is the more
  interesting outcome.
- **BELOW THE BAND** — residual < 0.494. Reported as such.

Secondary, reported regardless: raw AUROC; marginal value over descriptors alone
(paired bootstrap, 10,000 resamples) against Boltz-2's +0.093 [+0.062, +0.125]; and the same on
Factor Xa if Mpro completes.

## Controls — the arm is void without these

1. **Positive control.** The descriptor baseline on the same compounds must reproduce **0.7654**
   to ±0.005. If the panel is not the panel, no comparison is readable.
2. **Coverage floor.** ≥ 90% of 751 compounds scored, with failures reported by reason.
3. **Permutation.** Shuffled labels must give AUROC in [0.45, 0.55].

## Leakage — checked BEFORE the residual is read

FlashBind's binary task ships against MF-PCBA (PubChem-derived). Our panel is ChEMBL-derived
SARS-CoV-2 Mpro. Mpro is a heavily used benchmark target and these sets may intersect.

Before any number is interpreted, establish whether Mpro, or our specific compounds, appear in
FlashBind's training data. If the training set cannot be established, **the result is reported
with leakage unresolved and is not used to support the generalisation claim** — a model trained
on the target it is being tested on would clear the ceiling for the wrong reason.

## Abort conditions

- ESM3 representations cannot be generated (gated weights) → report blocked, do not substitute a
  different encoder and call it FlashBind.
- Checkpoints unavailable or licence-incompatible → report blocked.

## Pre-committed disclosure

Reported either way, including the outcome that weakens §3.4. The manuscript is not yet submitted
to a journal; if this lands before submission it goes in, and if it contradicts the current text
the text changes.
