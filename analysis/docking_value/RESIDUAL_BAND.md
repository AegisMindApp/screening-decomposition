# The docking residual band, recomputed under one procedure

**6 September 2026.** `residual_band.py`, results in `RESIDUAL_BAND.json`. Written because the
band quoted in the manuscript mixed two residualisation procedures and could not be defended as
a like-for-like interval.

## The problem

`MARGINAL_VALUE.md` residualised with out-of-fold **linear** regression; `SIZE_BIAS.md` used
out-of-fold **gradient boosting**. The published band 0.507–0.555 drew from both, and the
Boltz-2 figure it was compared against (0.6418) was gradient-boosted. Different estimators
answer different questions — a linear residual leaves nonlinear descriptor dependence in the
remainder and inflates it — so the comparison was not clean.

## All seven scores, one procedure

Residual AUROC = AUROC of `score − E[score | seven descriptors]`, predictions out-of-fold,
5-fold. Because the fold assignment moves the estimate, every value is the mean over **12 fold
seeds** with the full range given.

| score | n | raw AUROC | residual (GB), 12-seed mean | range |
|---|---|---|---|---|
| Vina, repaired receptor (Mpro) | 751 | 0.4530 | 0.5187 | [0.5131, 0.5259] |
| Vina, raw receptor (Mpro) | 751 | 0.4079 | 0.5092 | [0.4980, 0.5222] |
| gnina Vina term (Mpro) | 743 | 0.3791 | 0.4941 | [0.4887, 0.5017] |
| gnina CNNaffinity (Mpro) | 743 | 0.5389 | 0.4994 | [0.4908, 0.5091] |
| gnina CNNscore (Mpro) | 743 | 0.5041 | 0.5426 | [0.5320, 0.5592] |
| Vina, repaired receptor (Factor Xa) | 854 | 0.6775 | 0.5733 | [0.5634, 0.5810] |
| **Boltz-2 `prob_binary` (Mpro)** | 750 | 0.7913 | **0.6567** | **[0.6439, 0.6648]** |

**Docking band: 0.494 – 0.573.** Boltz-2 sits at 0.657, outside it.

## The separation survives the worst case

The claim is not that two point estimates differ. Across all 12 seeds the **highest** docking
residual observed is **0.5810** and the **lowest** Boltz-2 residual is **0.6439**. The intervals
do not touch under any fold assignment tried.

It also survives the choice of estimator. Under linear residualisation the docking band is
0.418–0.614 and Boltz-2 is 0.748 — wider on both sides, same conclusion.

## What did not reproduce, stated plainly

The two published gradient-boosted values are close but **not** reproducible to within seed
noise:

| | published | recomputed (12-seed range) |
|---|---|---|
| Vina, repaired (Mpro) | 0.5382 | [0.5131, 0.5259] |
| Boltz-2 | 0.6418 | [0.6439, 0.6648] |

Both published values fall just outside the recomputed ranges, in opposite directions. The
original run did not record its fold seed or regressor hyperparameters, so the discrepancy cannot
be attributed. The recomputed values are the ones that reproduce from a fixed script and are what
the manuscript now reports; `SIZE_BIAS.md` is corrected at the point of claim.

The correction moves both values in the direction that flatters us — Vina's residual down
(0.5382 → 0.5187), Boltz-2's up (0.6418 → 0.6567), widening the gap from 0.104 to 0.138. That is
worth stating explicitly: a correction that happens to help the authors deserves more scrutiny
than one that hurts them, which is why the seed sweep and the worst-case separation above are
reported rather than a single pair of point estimates.

## Effect on the manuscript

The old band 0.507–0.559 is superseded by **0.494–0.573**, wider at both ends and computed one
way. The hedge that the band was "indicative rather than like-for-like" is no longer needed and
is removed. The single-target limitation stands: there is still no Factor Xa Boltz-2 residual.
