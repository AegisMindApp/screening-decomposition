# Does fixing the receptor's hydrogen-bond donors change the Mpro result?

**Pre-registered 31 August 2026, before the re-dock was run.** Committed ahead of execution so
the reading rule cannot move once the number exists.

## The defect being repaired

Every receptor in this project was built by a converter that mapped **every nitrogen to `NA`
(hydrogen-bond acceptor)** and never added hydrogens. AutoDock Vina identifies a donor by its
attached polar hydrogen, so the donor set was empty by construction: the scorer could not form a
protein-donor → ligand-acceptor hydrogen bond, and counted 597 backbone amide NH groups as
acceptors.

## What changed, and what did not

Repaired with `obabel -xr -p 7.4` applied to the **exact coordinates that were docked**, not to a
fresh download of 7VU6 — a re-fetch could differ in waters, altlocs or occupancy and would
confound protonation with a change of structure.

| | before | after |
|---|---|---|
| `HD` polar hydrogens | **0** | **946** |
| `N` donor nitrogens | **0** | **727** |
| `NA` acceptor nitrogens | 754 | 27 |
| heavy atoms | 4520 | 4520 |

**Asserted before use: the heavy-atom coordinate multiset is identical.** obabel reordered atoms
and moved none, so protonation is the only difference. Everything else is held: the same 753
compounds, the same labels, the same box, the same vina binary, exhaustiveness 4.

## Baseline

Arm 1 at exhaustiveness 4 on this compound set: **AUROC 0.408** (0.3992 on the 743-compound
subset that also carries a gnina score). The comparison is paired within compound.

## Pre-registered reading

**ΔAUROC = AUROC(protonated) − AUROC(0.408)**, 95% CI from a paired bootstrap over compounds
(10,000 resamples).

| ΔAUROC 95% CI | reading |
|---|---|
| entirely **above 0** | receptor preparation was a real part of the failure. **Every docking number this project has published is affected** and must be re-run or re-labelled |
| entirely **below 0** | the unprotonated receptor was *better*, which would be evidence the H-bond term is not doing useful work on this target |
| **spans 0** | the defect is real but not the cause. Vina's H-bond term is not what is wrong here |

Secondary, reported regardless: Spearman ρ between the two score vectors, the top-10 shortlist
overlap, and whether the protonated run clears **0.5** and **0.763**.

## What no outcome licenses

**0.763 remains the bar.** Seven physicochemical descriptors reach it on these same compounds,
and a protonated run landing at 0.55 has still lost to a logistic regression.

One target, one receptor, one library. A positive result here would say preparation matters on
Mpro — it would *not* retroactively validate the Factor Xa or PD-L1 numbers, which were produced
with the same defect and would each need their own re-run.

And a null here does **not** clear the defect. Backbone amides typed as acceptors is wrong
regardless of whether Vina's scoring happens to be sensitive to it; the fix stays in either way.

## Validity conditions

- **≥95% of the 753 must dock**, or the surviving subset is not the same benchmark.
- ΔAUROC computed only on compounds scored in both runs.
- Failures recorded with reasons.
- The heavy-atom multiset assertion above must pass, or the run is void.
