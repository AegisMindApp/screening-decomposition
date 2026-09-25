# The Mpro floor re-measurement is sound — for the pocket-conditioned panel, which is not the one the paper uses

23 September 2026. `PREREGISTRATION_FLOOR_MPRO2.md` (amended to 35 pairs) ran to completion:
70/70 scored, zero failures, per-compound sd **0.02865**.

> **This file replaces an earlier version that declared the measurement void.** That was wrong.
> The run did replicate a different protocol from the *paper's* Mpro panel, which is the part
> worth knowing, but a panel under the matching protocol exists and the measurement is valid
> against it. Correcting rather than deleting, since the first reading is the one a reader
> would reach too.

## The gate fired, and chasing why found a protocol mismatch

The pre-registration rejects a sample unrepresentative of the panel on distance-from-bound, the
variable Boltz-2's noise scales with. The 35 compounds were stratified against the paper's Mpro
panel and match it almost exactly on that basis — 0.2022 against 0.2024 — yet the same compounds
*as measured in the replicate run* sit at 0.2557, a ratio of 1.263 and a clear fail. Simulating
the panel perturbed by the measured within-run noise gives 0.2039 [0.1952, 0.2124], so the gap is
not an artefact of the concave `min(p, 1−p)` transform.

The cause is that two Mpro panels exist and they were run under different conditioning:

| panel | pocket | mean distance-from-bound | is it the paper's? |
|---|---|---|---|
| arm 1, `analysis/boltz2/results/` | **blind** | 0.2024 | **yes** — matches `boltz2__Mpro.json` to 0.000000 |
| arm 2, `analysis/boltz2/results_arm2/` | **conditioned** | 0.2539 | no |

The floor run is pocket-conditioned: `push_boltz2_kaggle.py` sets `BOLTZ_POCKET=1` and the bundle
ships a `pocket.json`. So it replicates arm 2's protocol, not arm 1's.

## Against the right panel it passes every check

| check | vs arm 1 (blind) | vs arm 2 (conditioned) |
|---|---|---|
| representativeness ratio | 1.263 **fail** | **1.007 pass** |
| bias, replicate mean − panel | **+0.0952** | **+0.0004** |
| sd of that difference | 0.0824 | 0.0345 |
| within-run per-observation sd | 0.02865 | 0.02865 |
| correlation | 0.9746 | 0.9936 |

Against arm 2 the between-run difference is 1.2× the within-run noise — the same relationship
Factor Xa shows (0.0098 against 0.01059, 0.93×), which is what a matched protocol looks like.
Against arm 1 it is 2.9× with a systematic shift.

**So: 0.02865 is a sound per-compound floor for pocket-conditioned Mpro. It is not the floor for
the panel the manuscript scores.**

## What this settles, and what it leaves open

**Settled — the cross-target transfer verdict, and it strengthens.** All four floor runs were
pocket-conditioned, and each sample is representative of its own matching panel (Factor Xa 1.006,
Mpro 1.007), so the comparison is like for like:

| target | pairs | per-compound sd |
|---|---|---|
| Factor Xa | 100 | 0.01059 |
| Factor Xa | 24 | 0.01542 |
| Mpro | 23 | 0.01937 |
| **Mpro** | **35** | **0.02865** |

**2.71×** between the two best estimates, against 1.83× from the 23-pair figure and a
pre-registered ±25% band. `DOES NOT TRANSFER` was right and is now further from the band, not
nearer it. Boltz-2's noise constant is target-specific — unlike Vina's, whose 1/n panel-size
rescaling transfers between targets to within 9% (`analysis/seed_floor/RESULT_FXA_SCALING.md`).
Those two facts are about different things and both hold: panel-size rescaling transfers, the
per-target noise level does not.

**Open — the floor for the manuscript's blind Mpro panel.** Nothing measures it. Re-running the
same 35 compounds with `BOLTZ_POCKET=0` would, at roughly 9 GPU-hours; the stratification is
already correct on the blind panel's basis and only the conditioning changes.

## This has a consequence for the manuscript

§3 rescales the Factor Xa Boltz-2 floor by √(n_FXa/n) and charges it against Mpro and ALDH1. The
rescaling corrects for panel *size*. It does not correct for target, and the table above says the
target term is **2.71×** — larger than the size term being corrected for. Applying a Factor Xa
floor to Mpro is the same class of substitution the paper warns about for panel size, and the
paper should say so where it does it.

## A latent bug found on the way, now fixed

`boltz2_worker.py` fell back, when a bundle had no `pocket.json`, to
`[("A", r) for r in (164,145,41,189,165,143,142,49,144,166,141,26)]` — **PDB residue numbers**,
where the YAML field takes **1-based indices into the supplied sequence** (offset 2 here, as
`pocket.json` states: *"Boltz-2 receives a sequence and never sees PDB numbering"*). Those indices
select Ser, Ile, Pro and Ala rather than the Cys145/His41 dyad.

No published result used it — the Mpro panel ran blind and every current bundle ships a
`pocket.json` — but it would have silently mis-conditioned the next target added without one. The
fallback now raises instead of guessing, because a correct default cannot be written for a target
that does not exist yet.


## Postscript: what the propagation control says, including about itself

`floor_propagate.py` was written to make these floors re-derivable from committed code, and its
control run reports:

| noise model | quantity | sd | 95% upper |
|---|---|---|---|
| uniform | residual | 0.00407 | 0.00443 |
| uniform | **margin** | **0.00934** | **0.01018** |
| level-aware | residual | 0.00401 | 0.00437 |
| level-aware | margin | 0.00898 | 0.00979 |

Two things follow.

**The published 0.00440 / 0.00472 is the sd of the residual, not of the margin.** The margin —
residual minus its descriptor null, which is the quantity every verdict is actually charged
against — has an sd **2.2× larger**, because the residual and its null do not move together when
the scores are re-drawn. The manuscript charged margins against the residual's floor; corrected
in `docs/papers/resolution_v3_draft.md` §3.7, where it costs Mpro its verdict.

**The control does not reproduce the published figure exactly** — 0.00407 against 0.00440, 7%
low, on a noise model that is a reasonable reconstruction rather than the original. So the
published floor is still a number produced by code that is not in the repository. This file
narrows that gap and does not close it, and by the standard in `paper-audit` step 2d a headline
number that cannot be re-derived from git is unverified.
