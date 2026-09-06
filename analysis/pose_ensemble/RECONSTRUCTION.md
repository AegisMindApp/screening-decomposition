# The published pose-ensemble refutation was fold-luck. On a repaired receptor the arm clears the floor.

**7 September 2026.** `reconstruct.py`, results in `RECONSTRUCTION.json`. Rules fixed in
`PREREGISTRATION_RECONSTRUCTION.md` (`824270441`) before any estimate was read. 40
cross-validation fold seeds, every one reported, none selected.

## All four pipeline controls pass

Registered after the previous run's control validated the docking but not the feature pipeline —
which was where the risk was.

| control | repaired | defective | |
|---|---|---|---|
| **Degenerate identity** — ensemble restricted to its top-score feature must equal the top-pose arm | max \|Δ\| **0.00e+00** | **0.00e+00** | **PASS** |
| **Permuted labels** must give ≈0.5 | 0.4982 | 0.5047 | **PASS** |
| **Permuted features** must not gain | −0.0457 | −0.0980 | **PASS** |
| **Descriptor baseline** range must span the published 0.7651 | [0.7595, 0.7692] | same | **PASS** |

The degenerate-identity control is the one that was missing before. A pipeline that cannot
reproduce a result it contains as a special case cannot be trusted for the general case; this
one reproduces it exactly.

## H1 is confirmed: the disagreement was fold assignment

| | mean | sd | 95% seed range | full range |
|---|---|---|---|---|
| **Defective receptor** | +0.0312 | 0.0059 | [+0.0197, +0.0418] | [+0.0168, +0.0433] |
| **Repaired receptor** | **+0.0552** | 0.0070 | **[+0.0434, +0.0680]** | [+0.0400, +0.0681] |

**Both prior numbers are draws from the same distribution.** The published **+0.0194** sits at
the **2nd percentile** of the defective-receptor distribution; the reimplementation's **+0.0393**
sits at the **92nd**. Neither was a coding error. Two people ran the same analysis with different
fold seeds and drew from opposite tails.

## Verdicts against the pre-registered rules

| receptor | verdict |
|---|---|
| Defective | **UNRESOLVED** — the seed range [+0.0197, +0.0418] straddles the 0.0201 scoring floor |
| Repaired | **RESOLVABLE, ABOVE FLOOR** — the entire range sits above 0.0201 |

**The published refutation does not stand.** It was a single draw from a distribution that
straddles the floor, reported as though it settled the question. The honest reading of the
defective-receptor data was always "unresolved", not "refuted".

**On a correctly prepared receptor the arm clears.** +0.0552, with all 40 seeds above the floor
and the lowest draw at +0.0400 — twice the floor.

## The receptor effect, paired per seed

| | |
|---|---|
| repaired − defective, paired on identical folds | **+0.0240** ± 0.0094 |
| 95% range | [+0.0112, +0.0392] |
| seeds where repaired > defective | **39/40** |
| is the receptor effect itself above the scoring floor? | **no** — its lower bound is +0.0112 |

So the repair helps consistently but by an amount that is itself at the edge of resolvability.
The arm clearing on the repaired receptor is driven by the arm's baseline effect (+0.031) plus a
repair contribution (+0.024), neither of which alone is decisive.

## Consequence for the manuscript — the paper gets weaker

§3.2 reports the pose-ensemble arm as **+0.019, not resolvable**, and the abstract counts it
among the interventions falling inside the floor. That is now wrong on both counts.

Corrected, the arm is **+0.055 [+0.043, +0.068], resolvable, above the scoring floor.** Two
interventions clear rather than one, and this changes the paper's story rather than a number:

> "The two interventions that clear the floor concern *what is computed*, not *where or how
> thoroughly the ligand is placed*."

Pose-ensemble scoring **is** a placement-side intervention — it reads the pose distribution the
search produced. It clears. The clean dichotomy in §3.2 does not survive, and the honest
statement is narrower: search *effort* and pose *source* do not matter, but how the pose
distribution is *read* does.

The paper is more interesting for this and less tidy. It should not be submitted until §3.2, the
abstract and Figure 1 are corrected.

## Why this happened, recorded so it does not recur

The original analysis script was never committed — only its pre-registration and its write-up.
A number quoted forward from code that no longer exists cannot be checked, and this one was
wrong. `reconstruct.py` is committed **with** its results and its controls, and every seed it
drew is in `RECONSTRUCTION.json`.

Single-seed cross-validation estimates have now produced two misleading results in two days: the
residual band (0.01–0.02 seed sensitivity) and this arm (0.027 between the 2nd and 92nd
percentile). Any future AUROC difference of this size should be reported as a seed distribution,
not a point estimate.
