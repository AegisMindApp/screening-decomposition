# Pose ensemble on Factor Xa: INCONCLUSIVE, and that does not count as support

24 September 2026. Rules fixed in `PREREGISTRATION.md` (amendment) before the run.
Re-derive: `python analysis/replication_fxa/pose_ensemble_fxa.py` → `POSE_ENSEMBLE_FXA.json`.
Same 8 features, same RT, same 5-fold OOF logistic, same 10,000-resample paired bootstrap as
`analysis/pose_ensemble/analyse_protonated.py`. All 886 compounds carried poses, labels and
descriptors — no dropout.

## Both controls pass before the result is read

| control | Mpro | **Factor Xa** |
|---|---|---|
| degenerate identity — ensemble restricted to its top-score feature must equal the top-pose arm | 0.00e+00 | **0.00e+00 PASS** |
| permuted features vs top pose — must be clearly negative | negative | **−0.2400 [−0.2922, −0.1873]** |

The permuted control is worth reading, not just ticking. Shuffling the ensemble features across
compounds costs 0.24 AUROC, so the features carry substantial real signal on this target. The
null result below is not the features being inert.

## Result

| | Mpro (repaired receptor) | **Factor Xa** |
|---|---|---|
| raw Vina top-pose AUROC | 0.4490 | 0.6668 |
| top pose (logistic) | — | 0.6663 |
| ensemble, 8 features | — | 0.6528 |
| **ensemble − top pose** | **+0.0552 [+0.0434, +0.0680]** | **−0.0135 [−0.0315, +0.0043]** |
| floor charged | Mpro's | Factor Xa exh-4, **0.00703** at n=886 |

**INCONCLUSIVE.** The interval straddles both zero and the floor.

Per the pre-registration, written before any of this ran: *"A replication that returns
INCONCLUSIVE is reported as prominently as one that succeeds… it does **not** count as support."*
So this is not a second target agreeing with Mpro, and it is not a refutation either. The point
estimate is negative where Mpro's was positive, but the interval is too wide to make anything of
that, and reading a sign out of a straddling interval is exactly what the three-outcome framework
exists to prevent.

## Alongside the other completed arm

| intervention | Mpro | Factor Xa | verdict |
|---|---|---|---|
| scoring function (gnina CNN) | +0.1397 [+0.0913, +0.1867] | **−0.0583 [−0.0930, −0.0235]** | DOES NOT REPLICATE |
| pose ensemble | +0.0552 [+0.0434, +0.0680] | **−0.0135 [−0.0315, +0.0043]** | INCONCLUSIVE |

Both point estimates are negative on the second target. One clears its floor and one does not, so
only the first is a finding; but it is worth recording that neither replication came back
positive, because a reader will notice the pattern whether or not we name it.

A restraint worth stating: two arms is not a pattern yet, and the temptation to read one into
them is strong precisely because they agree in direction. The receptor-preparation and
search-effort arms are still running and will report against the same rules.

## Scope

One seed of Factor Xa poses (seed 1), whose 886 Vina scores are bit-identical to the committed
Kaggle seed-1 run. Both arms of the comparison read the same pose files, so pose-to-pose noise is
shared and cancels — the same design property that made the Mpro comparison paired.
