# Changing only the scoring function moves Mpro docking from 0.399 to 0.539 — across chance, but nowhere near the descriptor bar

**Run 31 August 2026** on Kaggle GPU (`mpro-gnina-rescore`). Reading rule pre-registered at
`cd8f55639` **before** these numbers existed. `analyse_arm3.py`, results in `ARM3_RESULT.json`.

## What was held constant

Everything except the scoring function. Arm 3 rescored **the same pdbqt files Vina wrote** —
identical atom coordinates — against the same receptor pdbqt (md5 asserted equal to the Vina
run). No re-docking, no re-preparation, no format conversion.

**745 of 753 scored (98.9%)**; 743 carry a Vina score, a gnina score and a label. Coverage 98.7%
against a pre-registered floor of 90%. The 8 failures are `CG0` macrocycle glue atoms that Vina
emits and gnina rejects, recorded with reasons.

**Positive control:** gnina reports its own Vina-style affinity for the pose it is handed. Against
our cached scores that correlates at **r = +0.970** — so gnina is demonstrably scoring the pose we
think it is. The script asserts this before reporting any AUROC.

## Result

| arm | poses | scorer | AUROC |
|---|---|---|---|
| 1 | Vina exh=4 | Vina | **0.3992** |
| 3 | **the same poses** | **gnina CNN** | **0.5389** |

**ΔAUROC = +0.1397, 95% CI [+0.0913, +0.1867]** (paired bootstrap, 10,000 resamples).

## Reading it against the rule, strictly

The pre-registered threshold was a CI **entirely above +0.10**. The point estimate (+0.140)
exceeds it; the CI's lower bound (+0.091) does not. So:

> **Directionally supported, not at the pre-registered strictness.**

The effect is unambiguously real — the CI excludes zero by a wide margin — and it is large enough
to carry docking across chance for the first time on this benchmark. But the specific ≥0.10 claim
is not confirmed at 95%, and reporting it as confirmed would be reading the point estimate and
ignoring the interval the rule was written around.

**What is confirmed outright** is the hypothesis's third clause: *neither arm will reach 0.763*.
Arm 3 lands at 0.539. Seven physicochemical descriptors still beat the best docking configuration
we have by **0.22 AUROC**.

## The finding underneath

Two scoring functions, given **the identical atom positions**, agree at only **r = +0.489**.

That is the substantive result. Half the variance in a docking score is not the pose — it is the
function used to read the pose. A pipeline can hold geometry completely fixed and still move from
below chance to above it by changing nothing but the scorer. Whatever Vina's affinity term is
measuring on this target, it is measuring it in a way that anti-correlates with activity, and the
poses are not to blame.

## What this does not settle

Arm 2 (DiffDock-L poses + Vina score) has not run, so *sampling* has not been ruled out as a
contributor — only shown to be unnecessary for a substantial gain. The hypothesis predicted arm 2
would move AUROC by less than 0.05; that remains open.

One target, one receptor, one library. And 0.539 is not a usable screening tool — it is a
below-descriptor result that happens to clear chance.
