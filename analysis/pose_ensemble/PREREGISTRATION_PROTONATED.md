# Pre-registration: pose-ensemble scoring on the REPAIRED receptor

**Written 6 September 2026, before the re-dock starts and before any outcome exists.**

## Why this run

`RESULT.md` measured pose-ensemble scoring at +0.0194 AUROC on the **donor-defective**
receptor and reported it as refuted. The manuscript currently carries that null as
*conditional on a receptor we later showed was broken* (§3.2, Supplementary Table S2). This
re-runs the arm on the repaired receptor so the null is unconditional, or is overturned.

The same docking run also produces the pose set needed to re-run the gnina scoring-function
arm on the repaired receptor. That arm is pre-registered separately, after this one lands.

## Protocol — identical to the recorded protonated run

Read from `bundle_CHEMBL4523582_PROTONATED/shard*/manifest.json`, not retyped:

| | |
|---|---|
| receptor md5 | `28657daa2db1856b36d66bf70651209b` (protonated; pre-repair was `d1bfc9ef…`) |
| exhaustiveness | 4 |
| num_modes | 5 |
| box centre | `[9.05, 8.9, -1.51]` |
| box size | `22 × 22 × 22` Å |
| vina md5 | `0c5d02550bdb3e661a77ee0badcfbe78` |

Only the receptor differs from `RESULT.md`'s run. Same compounds, same labels, same binary,
same box, same seed policy.

## Gate: the run is void unless the positive control passes

The top-scoring mode of this re-dock must reproduce the already-published protonated Vina
result, **AUROC 0.4530** (`protonated_scores_merged.json`), to within **±0.010**. Vina is
stochastic and this is a fresh run, so exact equality is not expected; a larger discrepancy
means this is not the protocol we think it is and no ensemble number may be read.

Receptor md5 is asserted equal to `28657daa…` before docking. Coverage floor: **≥ 90%** of 751
compounds scored, else void.

## Reading rules — fixed now

Ensemble features are built from the 5 modes exactly as in `RESULT.md`: per-compound best
affinity, mean affinity, spread, and RMSD-from-top, combined by out-of-fold logistic regression
under the same 5-fold splits used for the descriptor baseline.

1. **Primary.** ΔAUROC(ensemble − top pose) with a paired-bootstrap 95% CI over compounds,
   10,000 resamples. This arm **holds placement fixed** — both sides read the same docking run —
   so it is gated against the **scoring floor, 0.0201 [0.0112, 0.0286]**
   (`analysis/docking_value/FLOOR_INTERVAL.json`), not the protocol floor.
   - Δ **> 0.0201** *and* CI excluding zero → **supported**: pose ensembles carry resolvable
     information the top pose discards, and the original null was an artefact of the defective
     receptor.
   - otherwise → **refuted**, and the manuscript's null becomes unconditional rather than
     conditional.
2. **Comparison to the defective-receptor run.** ΔAUROC(this arm − 0.0194) reported with its
   own interval, so the effect of the repair on this arm is visible rather than inferred.
3. **Usefulness, secondary.** Ensemble + seven descriptors vs descriptors alone (0.7654). Reported
   regardless of rule 1.

## What no outcome licenses

One target, one receptor, one library, exhaustiveness 4. A positive result here would not
transfer to Factor Xa or to the exhaustiveness-32 protocol without running them. A negative
result does not show pose ensembles are useless in general — only that on this benchmark, at
this search effort, they do not clear the measured scoring floor.

## Pre-committed disclosure

Whatever this returns, it is reported. If it overturns the published null, the manuscript's
§3.2 and abstract change and the change is recorded as a correction rather than a silent edit.
