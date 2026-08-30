# Is our docking failure in the pose search or in the scoring function?

**Pre-registered 31 August 2026, before any DiffDock or gnina run.** Committed ahead of execution
so the reading rule cannot move once the numbers exist.

## The hypothesis

> On the SARS-CoV-2 Mpro retrospective benchmark (753 ChEMBL compounds: 257 measured actives and
> 496 measured inactives at pChEMBL<5, receptor PDB 7VU6), AutoDock Vina's below-chance ranking
> (AUROC 0.427) is a **scoring** failure and not a pose-**sampling** failure. Therefore:
> substituting DiffDock-L pose generation while retaining Vina scoring will change AUROC by less
> than 0.05, whereas rescoring poses with the gnina CNN scoring function will raise AUROC by at
> least 0.10 over the 0.427 baseline. Neither arm will reach the 0.763 AUROC obtained from seven
> physicochemical descriptors alone on the same 753 compounds.

## Why this is worth running

Docking has been measured on this benchmark twice and is below chance both times — 0.427 at
exhaustiveness 32, 0.408 at exhaustiveness 4. **Seven physicochemical descriptors reach 0.763 on
the same compounds.** So the physics is 0.34 AUROC worse than counting properties.

The obvious explanation is a size confound, and it is **wrong** — checked, not assumed. Actives
are 2.3 heavy atoms *larger*, larger molecules receive better Vina scores (r = −0.588), and larger
correlates weakly with being active (+0.134). Size should push docking to rank actives *better* on
both legs. It still ranks them worse (corr(score, active) = +0.104). That is a genuine anti-signal
and it is not obviously a property of the scoring function rather than the search.

DiffDock/EquiBind were named as the required fix twice in `paper_draft.md` (§4.3, §4.7) and never
run. This is that experiment.

## Arms

All four cells share the same 753 compounds, labels, receptor (7VU6) and box.

| arm | poses from | scored by | status |
|---|---|---|---|
| 1 | Vina (exh 32) | Vina | **done: 0.427** |
| 2 | **DiffDock-L** | Vina | to run |
| 3 | Vina (exh 32) | **gnina CNN** | to run |
| 4 | **DiffDock-L** | **gnina CNN** | to run |

Arm 1 vs 2 isolates **sampling**. Arm 1 vs 3 isolates **scoring**. Arm 4 is the combination the
LIT-PCBA benchmarks report as strongest, included so a null in 2 and 3 cannot be attributed to
having tested them only in isolation.

## Pre-registered reading

Primary contrasts, each with a 95% CI from a paired bootstrap over compounds (10,000 resamples):

| result | reading |
|---|---|
| ΔAUROC(2−1) CI contains 0 **and** ΔAUROC(3−1) CI entirely above +0.10 | **hypothesis supported** — the failure is in scoring |
| ΔAUROC(2−1) CI entirely above +0.05 | **hypothesis refuted** — pose sampling was a real part of the failure |
| both CIs contain 0 | **neither component is the problem.** Report it: the target or the benchmark is what defeats docking, and no docking method fixes it |
| any arm ≥ 0.763 | the descriptor bar is beaten — report which arm and treat as the headline |

Secondary, reported regardless: Spearman ρ between each arm's scores and Vina's, and EF@1%.

## What no outcome licenses

One target, one receptor, one library. **A negative here is not evidence that DiffDock or gnina
are poor methods** — it is evidence about this benchmark. Published benchmarks report DiffDock's
own confidence score at AUROC 0.538–0.56 for screening, so we do **not** use it as a ranking
score; DiffDock supplies poses only. Using its confidence head would test a quantity its authors
do not claim is an affinity estimate.

Beating 0.5 is not the bar. **0.763 is the bar**, and an arm that lands at 0.55 has still lost to
a logistic regression on molecular weight.

## Validity conditions, declared now

- **≥90% of the 753 must produce a pose and a score in every arm.** Below that the surviving
  subset is not the same benchmark; report the shortfall rather than the AUROC.
- All contrasts computed **only on compounds present in all four arms.**
- Failures recorded with reasons, never silently dropped.
- The receptor and the compound set are asserted identical to the exh=32 run by md5 before any
  arm is analysed.
