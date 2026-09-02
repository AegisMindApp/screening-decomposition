# Pre-registration — pose-ensemble scoring

**Written 2 September 2026, committed before any outcome number was computed.**

## Motivation

Search effort does not matter (ρ = 0.959 at 1× vs 8×, top-10 identical), which suggests the
best pose may be the wrong unit of analysis. Binding is an ensemble property: a compound with
one good pose and four terrible ones is physically different from one with five mediocre
poses. Current practice discards that. This is the only proposed direction from the external
review that our existing data could not already refute.

## Data

`analysis/retrospective_benchmark/kaggle/exh4_poses/` — 753 Mpro compounds, all Vina modes
retained with affinity and RMSD-from-top for each. **Exhaustiveness 4, original (donor-
defective) receptor.** All comparisons are therefore within-receptor and internally
consistent; this experiment cannot make claims about the repaired receptor.

## Features tested

Top-pose score (current practice) versus ensemble features derived from the same file:
Boltzmann-weighted free energy `−kT·ln Σexp(−ΔGᵢ/kT)` at T = 298 K; mean and standard
deviation of pose scores; score gap between best and second pose; count of poses within
1 kcal/mol of the top ("degeneracy"); and RMSD spread of those near-degenerate poses (well
width).

## Reading rules — fixed now

1. **Primary.** Ensemble features beat top-pose score by **> 0.04 AUROC** with a paired-
   bootstrap 95% CI excluding zero → supported. The 0.04 threshold is the measurement floor
   established in `../docking_value/MARGINAL_VALUE.md`; a gain below it is not resolvable on
   this benchmark and will be reported as **refuted**, not as "promising".
2. **Usefulness.** Even if rule 1 passes, ensemble + descriptors must beat descriptors alone
   (AUROC 0.7653) by > 0.04 with CI excluding zero. Failing this, the result is a curiosity
   about docking internals, not a screening method.
3. **Control.** A permuted-feature arm — ensemble features shuffled across compounds — must
   show no gain. If the permuted arm also gains, the pipeline is leaking and both results are
   void.

All arms use identical compounds and identical StratifiedKFold(5) folds, out-of-fold
predictions only.

## What would make me abandon the direction

Rule 1 failing. That is the whole hypothesis: that information exists in the pose
distribution which the top pose discards. If it does not clear the noise floor, there is
nothing here regardless of how the other arms behave.
