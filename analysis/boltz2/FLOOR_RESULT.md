# Boltz-2's own reproducibility floor: ~0.005, not the 0.0201 we borrowed

21 Sep 2026. Reading rules fixed in [`PREREGISTRATION_FLOOR.md`](PREREGISTRATION_FLOOR.md)
before the run.

## Stage A: Boltz-2 is stochastic

48 ligands = **24 Factor Xa compounds duplicated** as arms A and B — identical SMILES, same
protein, same MSA, same instance, `BOLTZ_CHUNK=8` so the copies were scored in different
`boltz predict` invocations. 48/48 scored, 0 failures.

| | |
|---|---|
| pairs identical to 5 dp | **0 / 24** |
| max abs difference, binding probability | **0.0905** |
| mean abs difference | 0.0097 |
| **per-compound score sd, single run** | **0.01542** |

So unlike gnina — 745/745 bitwise identical — Boltz-2 does **not** have a zero floor. The
pre-registered rule therefore required the floor to be measured rather than set to 0.

The noise is **heteroscedastic**: |difference| correlates with distance from the 0/1 bounds at
r = +0.80 (p < 0.001). Mid-range predictions are the noisy ones, which is what a bounded
probability output should do.

## The floor, on the quantity the verdict actually uses

The harness margin is **residual** AUROC minus a measured null, so the floor must be the sd of
the *residual* AUROC under re-running — not of the raw AUROC. Different quantity, different
noise. The measured per-compound perturbations were resampled onto the **real 885-compound
Factor Xa panel** (200 draws) and the residual AUROC recomputed each time.

| noise model | floor sd | 95% upper bound |
|---|---|---|
| uniform resampling of observed differences | 0.00499 | 0.00545 |
| level-aware (scaled by distance from bound) | **0.00543** | **0.00592** |

The 24-compound sample is representative of the panel on the variable that drives the noise —
mean distance-from-bound 0.1214 against the panel's 0.1221, a ratio of 1.006 — so neither model
is materially biased by the sample. The level-aware figure is the one to quote, being the larger.

## What this does to the verdict

**The borrowed Vina floor of 0.0201 was about 4x too large.**

| target | margin | vs measured floor 0.00592 |
|---|---|---|
| ALDH1 | +0.0600 | clears, 10x over |
| Mpro | +0.0571 | clears, 10x over |
| **Factor Xa** | **+0.0042** | **does not clear** |

**`NOT_BEYOND_PROPERTIES` stands.** But the character of the result changes. Against the
borrowed floor, Factor Xa failed by 5x and looked decisive. Against Boltz-2's own floor it fails
by a whisker — 0.0042 against 0.0054–0.0059. The verdict is the same; the confidence behind it
is much lower than the original numbers implied.

Also worth saying: the borrowed floor was wrong in the *conservative* direction for this
verdict, i.e. it made a negative look stronger than the evidence supported. A borrowed threshold
is not safe merely because it is large.

## Limits of this measurement

- **24 pairs is a small base** for the per-compound noise estimate, and the decision now turns
  on the third decimal place. Quadrupling to 100 pairs costs ~5.8 GPU-hours (~$5) and would
  tighten the floor substantially — better value than the 10 x 100 replicate design originally
  pre-registered, because the bottleneck is the noise estimate rather than the propagation.
- **Propagation assumes independent per-compound noise.** Correlated run-level drift is largely
  irrelevant to AUROC, which is rank-based and invariant to monotone shifts, so the independent
  component is the one that matters — but structured differences (a chemotype systematically
  scored differently on one run) would not be captured.
- The floor is specific to **(Boltz-2 2.2.1, L4, this protocol, the 885-compound FXa panel)**.
  Per `method_ladder.json::_floors_are_panel_specific`, re-scale by sqrt(n) before charging it
  against a different panel.

## Cost

One g2-standard-8 on-demand instance, 45 minutes of prediction, ~$1.40 — against the ~$15
originally budgeted for stage B, because the cheap determinism test came first.
