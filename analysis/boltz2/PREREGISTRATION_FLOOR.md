# Pre-registration: measuring Boltz-2's own reproducibility floor

Written 21 Sep 2026, **before any run**. The reading rules below are fixed now so that they
cannot be chosen after seeing which way the number falls.

## Why

`analysis/method_bench/RESULT_3TARGET.md` returns **NOT_BEYOND_PROPERTIES** for Boltz-2 on three
admissible targets, and the verdict turns on Factor Xa's margin of **+0.0042** against a floor
of **0.0201**.

That floor is not Boltz-2's. Traced to source, `0.0201` is `FLOOR_SCORING` from
`analysis/pose_ensemble/analyse_protonated.py` — the **Vina scoring floor**, for scoring-function
changes on Vina poses. `PREREGISTRATION_TRANSFER.md` adopted it for Boltz-2 as "the scoring
floor". So the headline verdict rests on a threshold borrowed from a different method on a
different panel: the same substitution `docs/papers/RESUBMISSION_GATE.md` blocks on.

## Design — stage A (determinism), the cheap decisive test

Precedent: `analysis/gnina_determinism/RESULT.md` settled gnina's floor at 0 by showing 745/745
bitwise-identical scores, rather than by expensive replication.

**Construction.** One bundle of 48 ligands built from **24 Factor Xa compounds, each appearing
twice** under different IDs (`FXAREP_A_*`, `FXAREP_B_*`). Same protein, same MSA, same instance.
`BOLTZ_CHUNK=8`, so the A copies and B copies are scored in **different `boltz predict`
invocations** — the thing a re-run actually varies. Duplicating within one bundle rather than
running twice avoids paying setup and MSA twice, and controls everything except the invocation.

**Reading rules, fixed now:**

- **All 24 pairs identical to 5 decimal places** → Boltz-2 contributes **no run-to-run noise
  under this protocol**; floor = **0.0, MEASURED** (not assumed). Stage B is not run.
- **Any pair differs** → Boltz-2 is stochastic; stage B measures the floor properly.
- **Fewer than 20 of 24 pairs recovered** → inconclusive, re-run. A partial result must not be
  read as agreement.

## What a zero floor would and would not license

It would **not** make +0.0042 meaningful. A zero floor says a re-run reproduces the number; it
says nothing about whether the number is distinguishable from chance on this panel. So if stage
A returns determinism, the margin is additionally charged against a **bootstrap 95% CI over
compounds** (10,000 resamples, seed 20260921), and the verdict is reported against that.
Clearing a floor of zero is close to a tautology and will not be presented as a result.

The verdict that follows mechanically from floor = 0, given the measured margins and cold
transfers already in `method_bench.json`, is **CONTRIBUTES, NOT DEPLOYABLE** — every increment is
positive but only 3 of 6 cold absolutes reach 0.60. Stating that now, before the run, so it
cannot look like a conclusion chosen to be interesting.

## Design — stage B (only if stochastic)

10 replicates of a fixed 100-compound Factor Xa subset. Floor reported as the χ² one-sided 95%
upper bound on the AUROC sd, then **scaled to Factor Xa's 885-compound panel by √(n/885)** —
a relationship verified to within 3–9% in `analysis/seed_floor/floor_exh4_vs_exh32.py`, not
assumed. Estimated ~17.5 GPU-hours, ~$15.

## Abort conditions

Stage A aborts if the instance fails to produce affinity output for the first 8 ligands, or if
cost exceeds 3 GPU-hours (`BOLTZ_CEILING_HOURS=3.0`).
