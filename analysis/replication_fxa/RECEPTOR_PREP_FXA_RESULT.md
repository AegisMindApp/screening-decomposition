# Receptor-preparation repair on Factor Xa: same sign, a quarter the size, INCONCLUSIVE

25 September 2026. Arm 1 of `PREREGISTRATION.md`.
Re-derive: `python analysis/replication_fxa/receptor_prep_fxa.py` → `RECEPTOR_PREP_FXA.json`.

## Result

| | Mpro | **Factor Xa** |
|---|---|---|
| ΔAUROC, protonated − unprotonated | **+0.0451 [+0.0236, +0.0671]** | **+0.0116 [+0.0001, +0.0234]** |
| floor charged | Mpro's | Factor Xa exh-4, **0.00753** at n=772 |
| per-seed deltas | — | +0.0028, +0.0193, +0.0128 |

**INCONCLUSIVE.** The effect has the **same sign** as Mpro's and the interval excludes zero by a
hair, but its lower bound (+0.0001) sits below the floor (0.00753), so the data do not establish
that the effect exceeds what reseeding alone produces. Per the pre-registration this is not
support.

This is the closest any arm has come to replicating. The direction agrees, the mechanism is the
same single documented change (`obabel -xr -p 7.4`, md5-verified as the only difference between
the two receptors), and the three per-seed deltas are all positive. What is missing is
resolution, not agreement: at a quarter the Mpro effect size, three seeds on 772 compounds cannot
separate it from the floor.

## The panel is truncated, and that was checked rather than waved through

The three unprotonated seeds hit **Kaggle's 12-hour session limit** and were cancelled at
875, 800 and 775 of 894 compounds. `run_seeds` writes its partial output as it goes, so the work
survived, but the surviving panel is *whatever finished first* — a deterministic docking order,
not a random sample.

| | value |
|---|---|
| common to both receptors and all three seeds | **772** of 886 |
| dropped | 114 |
| active fraction, kept vs dropped | 46.4% vs 41.2%, **−5.1pp** |

The pre-registered rule (borrowed from the Mpro Boltz-2 arm) voids a comparison whose dropout is
class-selected beyond 10pp. At −5.1pp it passes, and the script asserts it before computing
anything. Worth stating plainly: the dropout is **not** negligible, it is merely within the
declared tolerance, and a 114-compound truncation is a real limitation of this arm rather than a
formality.

## Where the replication now stands

| intervention | Mpro | Factor Xa | verdict |
|---|---|---|---|
| scoring function (gnina CNN) | +0.1397 [+0.0913, +0.1867] | **−0.0583 [−0.0930, −0.0235]** | **DOES NOT REPLICATE** |
| pose ensemble | +0.0552 [+0.0434, +0.0680] | −0.0135 [−0.0315, +0.0043] | INCONCLUSIVE |
| receptor preparation repair | +0.0451 [+0.0236, +0.0671] | +0.0116 [+0.0001, +0.0234] | INCONCLUSIVE |
| eight-fold search effort | +0.0082 [+0.0016, +0.0148] | *arm 2 running* | — |

Three arms report and **none replicates**: one reverses sign, two cannot be resolved on the
second target. The two inconclusive arms fail for opposite-looking reasons — pose ensemble points
the wrong way with a wide interval, receptor prep points the right way with too small an effect —
and it would be a mistake to average that into a single narrative.

The honest summary of what a second target has bought so far: **the one intervention large enough
to be decided on it went the other way**, and the smaller ones are below what a second panel of
this size can resolve. That is a statement about the reach of the replication as much as about
the interventions, and the paper's own framework is the reason it can be said precisely rather
than as "failed to replicate".

## Scope

Three seeds per receptor, exhaustiveness 4, paired within seed on identical compounds. The two
receptors differ in one documented step, with `receptor_md5_before_protonation` in the protonated
bundle matching the unprotonated bundle's `receptor_md5` exactly. Completing the truncated seeds
would need a resumable re-push; at this effect size it is unlikely to change the verdict, since
resolving +0.0116 against a 0.0075 floor needs a much larger panel rather than 114 more compounds.
