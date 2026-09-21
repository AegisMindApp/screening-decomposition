# Pre-registration: the run-to-run resolution limit of the Mpro docking protocol

Written before any seed replicate was run. Answers the defect that got `ci-2026-03151v`
declined by JCIM on 17 Sep 2026 without external review.

## What the editor said, and why it is right

The manuscript's decision rule rests on two "resolution limits". Neither measures what the rule
needs:

- **0.020** is the disagreement between two *implementations* of a scoring function. That is
  implementation variance.
- **0.039** comes from comparing two *search-effort conditions*. That is an intervention effect.
- And the fatal one: **the exhaustiveness comparison is both one of the six interventions and
  the threshold those interventions are judged against.** The intervention is evaluated against
  itself.

No rebuttal is appropriate. The floors are the decision rule for all six interventions, so this
is the spine of the paper.

## A third problem, found while building this

`kaggle_shard_worker.py:46` builds the Vina command with receptor, ligand, box, exhaustiveness,
num_modes and `--cpu` — and **no `--seed`**. AutoDock Vina without `--seed` draws one from
system entropy. So every docking result in the manuscript was produced at an *uncontrolled*
random seed, and seed variance is not merely unmeasured: it is already mixed into both reported
"limits" alongside the quantity each was supposed to isolate. That strengthens the editor's
objection rather than softening it, and it is reported rather than quietly fixed.

## The measurement

Hold **everything** fixed and vary **only** `--seed`:

- receptor 7VU6, `receptor_md5 = d1bfc9efe1cbf8f78592f72df9be99c4`
- vina binary, `vina_md5 = 0c5d02550bdb3e661a77ee0badcfbe78`
- box centre (9.05, 8.90, −1.51), size (22, 22, 22)
- `num_modes 5`, `--cpu 1`
- the same 753 compounds, the same labels (257 active / 496 inactive at pchembl < 5)

Both md5s were re-verified against the live bundle before this was written. The per-seed AUROC
is computed exactly as the manuscript computes it, and **the spread of AUROC across seeds is the
resolution limit** — a property of an unchanged protocol, contaminated by no intervention.

## Staging, and why

Measured on this machine: ~200 s of CPU time per ligand at exhaustiveness 32, so one full-panel
seed replicate costs ~42 core-hours and ten cost ~420. That is days of wall clock, so the work
is staged and the cheap stage is run first.

**Stage 1 — exhaustiveness 4, full panel, 10 seeds.** ~8× cheaper (~5 core-hours per seed).
Exhaustiveness is Vina's count of independent Monte-Carlo runs, so *less* search is *more*
stochastic: the exh=4 seed floor is an **upper bound** on the exh=32 floor. A conservative
number, measured, available in hours rather than days — and it applies directly to the
manuscript's exh=4 arm.

**Stage 2 — exhaustiveness 32, full panel, 10 seeds.** The floor for the protocol the headline
results use. Run across Kaggle and local CPU over several days.

Stage 2 also **tests the monotonicity assumption stage 1 relies on**. If the exh=32 floor is not
below the exh=4 floor, stage 1 was not the upper bound it is claimed to be, and that gets
reported.

## Statistic, fixed here

The point estimate is the **standard deviation of AUROC across seeds**. It is reported with its
**95% one-sided upper confidence bound**, from the χ² distribution on *n*−1 degrees of freedom,
and **the upper bound is what serves as the resolution limit**.

That choice is deliberate. An SD from 10 replicates carries roughly 24% relative standard error,
and quoting the point estimate as a threshold would overstate the precision of the floor itself —
the same species of error the editor identified. Using the upper bound means fewer seeds costs a
*wider, more conservative* floor rather than an unjustified claim, and any intervention that
clears it clears it honestly.

Also reported, because they are free once the runs exist: per-ligand score SD across seeds, the
number of ligands whose rank changes, and the full per-seed AUROC list.

## What this does and does not fix

It replaces the yardstick. **None of the six interventions needs re-running** — they are
re-read against a floor that is independent of them, which dissolves the circularity.

It does **not** address the editor's third point, chemical-series-controlled splits. That is a
separate piece of work and is not claimed here.

## What the floor will be compared against, fixed before it exists

Taken from the submitted manuscript so the comparison cannot be chosen after seeing the number.

| intervention | ΔAUROC | 95% CI | paper's verdict | floor it used |
|---|---|---|---|---|
| Scoring function (gnina CNN) | **+0.140** | [+0.091, +0.187] | directionally supported | 0.020 |
| Pose ensemble vs top pose | **+0.055** | [+0.043, +0.068] | supported | 0.020 |
| Receptor preparation repair | **+0.045** | [+0.024, +0.067] | supported | 0.039 (marginal) |
| Binding site supplied to Boltz-2 | −0.013 | [−0.026, +0.000] | **not demonstrated** | 0.039 |
| Eight-fold search effort | −0.019 | [−0.032, −0.006] | inside the floor | 0.039 |

### The risk is NOT one-directional, and that is the point

It is tempting to assume a smaller measured floor is good news. It is not uniformly:

- **A smaller floor rescues the positive findings.** +0.045 and +0.055 are currently marginal
  against a 0.039 floor; against a genuinely small floor they clear comfortably.
- **A smaller floor destroys the negative ones.** The paper dismisses −0.013 (Boltz-2 site) and
  −0.019 (eight-fold search) as *inside the noise*. That dismissal is only valid if the noise
  really is ~0.039. If the run-to-run floor is, say, 0.002, then **−0.019 is a real effect** —
  eight-fold more search measurably *hurts* — and "not demonstrated" becomes a false negative
  rather than a safe hedge.
- **A larger floor** (above ~0.055) removes the two supported findings and leaves only the
  +0.140 scoring result standing.

So there is no outcome that is simply favourable, and I am writing that down now precisely
because the temptation on seeing a small number would be to report only the half that helps.
Both directions get reported whichever way it lands.

A calibration, not a prediction: on a 400-compound synthetic panel, 0.05 kcal/mol of per-ligand
jitter produced an AUROC SD of 0.0009. If real seed jitter is of that order the floor is far
below every effect above — which is the branch that breaks the paper's two null results.

## OPEN QUESTION, 18 Sep 2026: is the seed floor a lower bound on reproducibility error?

**Stated as a question, not a finding — an earlier draft of this section overstated it.**

The gap below was first reported as "4.7σ", computed against the **point** SD. Against the
**95% upper bound** — the statistic this pre-registration commits to precisely because an SD
from few replicates is imprecise — the same gap is **2.1σ**, which is unremarkable. Using the
point estimate to call something anomalous is the identical overstatement of precision that got
the manuscript declined, committed here by me, two sections below where I wrote the rule
forbidding it.

At n = 10 the bound tightens from 2.4× to 1.65× the point estimate and this resolves either way.
Until then the material below is a hypothesis under test, not a result.

The manuscript's own exh=4 run is still on disk per-compound. Compared against my five seed
replicates on the **same 743 compounds with the same labels**:

| | value |
|---|---|
| pairwise per-compound score sd, mine vs mine | **0.2283 eV** |
| pairwise per-compound score sd, mine vs theirs | **0.2286 eV** — identical |
| AUROC, my five seeds | 0.4182 – 0.4237 |
| AUROC, the manuscript's run | **0.4059** |

At the score level their run is an ordinary seed replicate: the scatter is the same to three
decimals. But its AUROC sits **~4.7 seed-SDs below** my five, and the mechanism is a class
asymmetry — their run scores actives 0.036 eV worse *relative to inactives* than mine, which is
about 2σ on that differential.

**So seed variance is not the whole of run-to-run variability.** My floor varies only `--seed`
within one runner on one machine; two independent executions of the nominally identical protocol
differ by ~0.012–0.018 AUROC, several times that. The floor measured here is a **lower bound on
total reproducibility error**, and quoting it as *the* resolution limit would repeat the
manuscript's own error in the opposite direction — it named implementation variance as the
limit and ignored seed; I have measured seed and would be ignoring implementation.

**This changes the verdicts.** Against a seed-only floor of 0.0074 both negative findings sit
outside the noise. Against an implementation-inclusive figure nearer 0.018:

| intervention | vs seed floor 0.0074 | vs ~0.018 |
|---|---|---|
| +0.140 scoring | clears | clears |
| +0.055 pose ensemble | clears | clears |
| +0.045 receptor repair | clears | clears, marginally |
| −0.013 Boltz-2 site | outside | **inside — back to "not demonstrated"** |
| −0.019 search effort | outside | **marginal** |

The three positive findings survive either way. The two negative ones do not have a stable
answer yet, and **the paper's headline −0.019 search-effort effect is the same size as the
irreproducibility between two runs of its own exh=4 condition** — which is the more serious
observation of the two.

Resolving it needs replicate runs of the *whole pipeline* on independent machines, not just
`--seed` variation, and that is a larger measurement than stage 2 as scoped. It is not claimed
as done.

## Reading rule

- If the measured floor is **below** the effects the manuscript calls real, those findings
  survive with a defensible threshold.
- If it is **above** them, those findings do not survive, and the paper's conclusions change.
  That outcome is to be reported, not re-run with different settings until it is not.

Recorded now, before any seed replicate exists, precisely because the second branch is the one
that would be tempting to relitigate.
