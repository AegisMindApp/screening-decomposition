# Eight-fold search effort on Factor Xa: +0.0096, six of six seeds positive, INCONCLUSIVE

25 September 2026. Arm 2 of `PREREGISTRATION.md`. Six seeds at exhaustiveness 32 on GCP spot,
paired within seed against the existing exhaustiveness-4 seeds on 861 common compounds.
Re-derive: `python analysis/replication_fxa/search_effort_fxa.py` → `SEARCH_EFFORT_FXA.json`.

## Result

| | Mpro | **Factor Xa** |
|---|---|---|
| ΔAUROC, exh 32 − exh 4 | +0.0082 [+0.0016, +0.0148] | **+0.0096 [+0.0049, +0.0142]** |
| floor charged | Vina exh-32, 0.01025 | **Factor Xa's own exh-32, 0.00517** |
| verdict | INCONCLUSIVE | **INCONCLUSIVE** |

Per-seed deltas: +0.0094, +0.0030, +0.0116, +0.0145, +0.0131, +0.0057 — **positive in six of
six**, and the point estimates agree across targets to within 0.0014.

This is the closest agreement of any arm, and it still does not clear the bar. The interval's
lower bound (+0.0049) sits just under Factor Xa's own exh-32 floor (0.00517) — short by
**0.0003**. On the pre-registered rule that is INCONCLUSIVE, and it is recorded as such rather
than rounded into a replication.

Two honest observations about that margin, which point in opposite directions:

- The effect is almost certainly real. Six of six seeds positive, on two independent targets,
  with point estimates 0.0082 and 0.0096. Nothing about this looks like noise.
- It is nonetheless **not resolvable** by the benchmark it was measured on, which is the
  paper's entire subject. An effect can be real and beneath the resolution of the instrument
  used to claim it, and this arm is a clean instance.

The floor here is stricter than Mpro's (0.00517 against 0.01025) because Factor Xa's exh-32
AUROC sd over six seeds is 0.00247 — the run is more reproducible on this target. So Factor Xa
gave this intervention its *best* chance of clearing, and it still did not.

## The floor is Factor Xa's own, measured from these seeds

Not borrowed from Mpro, and not Factor Xa's exhaustiveness-4 floor. The six exh-32 seeds supply
the sd directly (0.00247), and the χ² one-sided 95% upper bound at n = 6 carries a 2.09× penalty,
giving 0.00517. At n = 6 that bound is loose in the permissive direction — more seeds would
tighten it and make clearing *harder*, not easier, so the INCONCLUSIVE verdict is not an artefact
of thin replication.

## The replication is now complete: four arms, none replicating

| intervention | Mpro | Factor Xa | verdict |
|---|---|---|---|
| scoring function (gnina CNN) | +0.1397 [+0.0913, +0.1867] | **−0.0583 [−0.0930, −0.0235]** | **DOES NOT REPLICATE** |
| pose ensemble | +0.0552 [+0.0434, +0.0680] | −0.0135 [−0.0315, +0.0043] | INCONCLUSIVE |
| receptor preparation repair | +0.0451 [+0.0236, +0.0671] | +0.0116 [+0.0001, +0.0234] | INCONCLUSIVE |
| eight-fold search effort | +0.0082 [+0.0016, +0.0148] | **+0.0096 [+0.0049, +0.0142]** | INCONCLUSIVE |

Four of the five interventions now have a second target. The fifth (binding site to a co-folding
model) remains single-target and is the one whose Mpro verdict was already INCONCLUSIVE.

**No arm replicates.** The pattern is not uniform and should not be reported as though it were:

- the **largest** Mpro effect reverses sign;
- the two mid-sized effects shrink to roughly a quarter and cannot be resolved, one of them
  pointing the wrong way;
- the **smallest** effect reproduces almost exactly — and neither target can resolve it.

The last row is the one a reader should find most interesting. An intervention whose effect
agrees to three decimal places across two independent targets is still, on both, indistinguishable
from re-running the same protocol with a different seed. That is the paper's thesis arriving in
its strongest form, on the paper's own data, without any appeal to the literature.

## Scope

Six seeds per arm, paired within seed, 861 compounds scored in all twelve runs. The exh-32 runs
were executed on GCP spot instances rather than Kaggle; Vina with a fixed seed and `--cpu 1` is
deterministic and the exh-4 arm's GCP/Kaggle agreement was verified bit-identical on all 886
compounds, so the platform is not a variable. Nothing here addresses whether either protocol is
any good: Vina sits at 0.67 on this panel throughout.
