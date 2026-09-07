# Pre-registration: the DiffDock pose-source arm on the repaired receptor

**Written 7 September 2026, before any score is computed on the repaired receptor.**

## Why

§3.2 reports the pose-source arm at **+0.017 [−0.027, +0.064]**, measured on the
donor-defective receptor, and it is now the **only** conditional null left in the manuscript.
The pose-ensemble arm flipped under exactly this treatment — a published refutation turned out
to be a 2nd-percentile fold draw on a broken receptor — so this null cannot be left untested.

## No GPU run is needed

DiffDock consumes a heavy-atom receptor PDB; protonation adds hydrogens and does not change its
pose generation. The 746 `rank1.sdf` poses from the original 31 August run were recovered from
the Kaggle kernel outputs, matching that arm's published coverage of 746 exactly. Only the
downstream Vina minimisation and scoring change, and both are CPU work.

## Protocol

Each DiffDock rank-1 pose is minimised in place against the **protonated** receptor
(`md5 28657daa…`) with `vina --local_only`, then scored — the same two-step treatment the
published arm used, changing only the receptor. The Vina-pose comparator is the protonated
exhaustiveness-4 run already completed (751/751, zero failures).

## Controls — the run is void unless both pass

Registered at pipeline level, not just at the docking level. The previous arm's control validated
the docking and missed a feature-pipeline problem.

1. **Receptor identity** — asserted equal to `28657daa2db1856b36d66bf70651209b` before scoring.
2. **Comparator reproduction** — the Vina-pose arm on the protonated receptor must reproduce
   its known AUROC of 0.4530 to within ±0.005 on the shared compound subset.

Coverage floor: ≥ 90% of the 746 recovered poses must minimise and score.

## Reading rules — fixed now

Primary: **ΔAUROC(DiffDock poses − Vina poses)**, both scored by Vina on the repaired receptor,
paired over the shared compounds. Because this arm changes *placement*, it is gated against the
**protocol floor, 0.0393**, not the scoring floor.

Reported across **40 cross-validation fold seeds** where any fitted model is involved, and as a
paired bootstrap over compounds (10,000 resamples) for the raw AUROC contrast. Single-seed point
estimates are not read: fold assignment has already produced two misleading results in this
project.

- **|Δ| > 0.0393 with a CI excluding zero** → the pose source matters on a correctly prepared
  receptor, and §3.2's null is overturned.
- **otherwise** → the null stands, and it becomes unconditional rather than conditional.

## Pre-committed disclosure

Reported either way. If this overturns the null, §3.2 loses a second row and the manuscript's
"pose source is inert" claim goes with it.

---

## Amendment 1 — 7 September 2026, after diagnosing the failures, before reading any effect

**The 90% coverage floor was unachievable by construction and is mine to correct.**

The failures are not conversion errors. Vina 1.2.7 rejects them with *"The ligand is outside the
grid box"* — DiffDock's pocket-conditioned top pose still lands outside the 22 Å box for roughly
a fifth of compounds. The run's failure classifier looked for the string "not in the search
space" and so mislabelled every one of them as `no_score`; that is corrected.

**The published arm has the same property.** `minv_diffdock.json` holds **604 of 751** minimised
poses — **80.4%** coverage. A 90% floor could never have been met by this arm on this data, and
setting it reflected my not having checked what the original achieved.

**Amendment:** the coverage floor becomes **75%**, below the published arm's own 80.4% and above
the ~79% this run is tracking. No reading rule, floor or comparison changes.

**A control this makes necessary, added now.** Excluding a fifth of the panel is not random: the
compounds DiffDock places outside the pocket may differ systematically from those it places
inside, and this project has previously been bitten by size-selected dropout. Before any effect
is read:

- the **active fraction among excluded compounds** must be within 10 percentage points of the
  panel's 34.2%, otherwise the exclusion is label-biased and the arm is reported as void;
- the primary comparison is **paired on the shared subset only** — DiffDock and Vina poses for
  the same compounds — and the subset size is reported alongside the effect;
- the excluded set's molecular-weight distribution is reported against the panel's, so a
  size-selected exclusion is visible rather than assumed absent.

This is a floor corrected because it was set without checking the achievable value, not because
the observed number was inconvenient. The distinction is the same one drawn in Amendment 1 of
`../docking_audit/PREREGISTRATION.md`, and it is the only kind of amendment that is legitimate
after data exists.
