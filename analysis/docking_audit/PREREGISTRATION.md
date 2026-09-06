# Pre-registration: retroactive audit of every historical docking claim

**Written 6 September 2026, before any panel other than Mpro and Factor Xa has been scored
with these diagnostics.** Mpro and Factor Xa results already exist and are used as positive
controls, not as findings.

## What this is

Not an attempt to rescue anything. This applies three diagnostics developed this week to every
panel this project has docked, to establish which historical docking claims were readable at
all. The expected outcome is that the number of standing claims **falls**.

## Panels in scope

Every target with labels, SMILES and at least one score vector on disk:

| target | n | active % | audited before? |
|---|---|---|---|
| CHEMBL4523582 Mpro | 751 | 34.2 | yes — **positive control** |
| CHEMBL244 Factor Xa | 886 | 45.7 | yes — **positive control** |
| CHEMBL3580522 PD-L1 | 418 | 79.9 | **no** |
| CHEMBL612545 | 2249 | 3.3 | **no** |

## Diagnostics, all pre-specified

1. **Descriptor baseline** — seven free descriptors, out-of-fold logistic regression, 5-fold.
   Bias control in the sense of Sieg et al.
2. **Raw docking AUROC** on the same compounds.
3. **Residual AUROC** — docking score with the descriptors regressed out, out-of-fold
   gradient-boosted, averaged over 12 fold seeds with the range reported.
4. **Marginal value** — AUROC(descriptors + docking) − AUROC(descriptors), paired bootstrap
   over compounds, 10,000 resamples.

## Reading rules — fixed now

Applied per (target, score vector). Floors from `../docking_value/FLOOR_INTERVAL.json`:
scoring **0.0201**, protocol **0.0393**.

- **PROPERTY-CONFOUNDED** — descriptor baseline ≥ 0.70. The panel is substantially separable
  from ligand properties alone, so any docking conclusion drawn on it is not evidence about
  binding. Reported regardless of what docking did.
- **NO MARGINAL VALUE** — the marginal-value CI includes zero, or the point estimate is below
  the protocol floor. Docking added nothing resolvable over free descriptors.
- **NO TARGET-SPECIFIC SIGNAL** — residual AUROC 95% range includes 0.50. The part of the
  docking score orthogonal to ligand properties does not rank actives.
- **SURVIVES** — none of the above fires.

A panel with a class imbalance beyond 75/25 is additionally flagged **IMBALANCED**: AUROC is
still interpretable but enrichment and the descriptor baseline both inflate, and no conclusion
is drawn from it without that caveat.

## Gate

The audit is void unless both positive controls reproduce to within **±0.005**:

| control | expected |
|---|---|
| Mpro, Vina repaired, raw AUROC | 0.4530 |
| Mpro, descriptor baseline | 0.7654 |
| Factor Xa, Vina repaired, raw AUROC | 0.6775 |
| Factor Xa, descriptor baseline | 0.7041 |

If a control fails, no new row is readable and the discrepancy is investigated first. A sweep
that silently produces plausible numbers on a broken join is the failure mode this gate exists
to prevent.

## Pre-committed disclosure

Every panel that enters the sweep is reported with its verdict, including panels that survive
and panels whose failure is inconvenient. Results feed back to the discovery review queue and
to solver.press where a public claim rests on an audited panel.

---

## Amendment 1 — 6 September 2026, after the first run, before any re-run

**The gate failed as written and the failure was mine, not the pipeline's.**

The Factor Xa descriptor control was specified as **0.7041**. That value belongs to the
**886**-compound panel used for the Boltz-2 Factor Xa arm (`analysis/boltz2/FXA_RESULT.md`).
The sweep joins `CHEMBL244_PROT` labels against `fxa_protonated_all.json`, which yields the
**854**-compound panel — the one used in `analysis/docking_value/MARGINAL_VALUE.md`, whose
published descriptor baseline is **0.7129**.

The sweep computed **0.7107** against that panel: a discrepancy of 0.0022, inside the ±0.005
tolerance. It also reproduced the panel's raw Vina AUROC exactly (0.6775) and its marginal value
closely (+0.0172 against a published +0.0150). The pipeline was correct; the expected value was
copied from the wrong panel.

**Amendment:** the Factor Xa descriptor control becomes **0.7129 ± 0.005**, the value published
for the 854-compound panel the sweep actually uses. No other control, threshold or reading rule
changes.

This is recorded rather than edited in place because a control that is relaxed after seeing the
data is worthless. What justifies the change is that the target value was misattributed to the
wrong compound set, verifiable against `MARGINAL_VALUE.md:15` — not that the observed number was
inconvenient. The first run's outputs are superseded by the re-run and no result below the gate
was read while the gate was failing.

**Also fixed:** the first run crashed writing `SWEEP.json` (`numpy.bool_` is not JSON
serialisable). The table had already printed, so the results above were visible before the
crash; they are recomputed rather than transcribed.
