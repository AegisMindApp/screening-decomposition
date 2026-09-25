# Pre-registration: does the floor hold for the metrics the field actually reports, and at AUROC levels other than ours?

Written 24 September 2026 **before either measurement runs**. Both use only data already
committed; neither costs compute. They exist because an examiner can refuse the paper's
conclusion on two grounds it currently does not address.

## Objection 1 — the paper measures one metric, and it is not the one most screening papers report

The manuscript measures a floor for **AUROC** and concludes about "virtual screening
improvements". Most screening papers report **enrichment** at an early cutoff — EF@1%, EF@5% —
or BEDROC, precisely because AUROC weights the whole ranking and screening only cares about the
top. A floor for AUROC says nothing about a floor for EF@1%, and EF@1% on a 755-compound panel
is decided by roughly **eight compounds**, so it should be far noisier.

**Measurement.** From the same ten exhaustiveness-4 seeds and ten exhaustiveness-32 seeds,
compute the seed-to-seed sd and χ² 95% upper bound for EF@1%, EF@5% and BEDROC(α=20) alongside
AUROC. Nothing is re-docked.

**Prediction, fixed now:** EF@1%'s floor **relative to its own scale** will be several times
AUROC's. Concretely, the coefficient of variation (sd ÷ mean) of EF@1% will exceed that of
AUROC by more than 3×. If it does not, the paper's AUROC-only scope is less of a limitation than
claimed and that should be said.

**What either outcome licenses.** A large EF floor does **not** weaken the paper's AUROC results;
it bounds their scope, and the scope claim in the title and discussion must then be narrowed to
AUROC or the EF floor reported alongside. A small EF floor would broaden the claim.

## Objection 2 — the floor was measured where the method has no signal

Vina's mean AUROC on this panel is **0.4175** at exhaustiveness 4 and 0.4531 at 32 — **below
chance**. A referee running a benchmark where docking reaches 0.75 can say the floor measured
here does not apply to them, and they would be raising a real question: AUROC's sensitivity to
score perturbation depends on how many compound pairs sit close together in score, which depends
on the level.

**Measurement.** Hold the panel, the class balance, the score *multiset* and the measured
per-compound seed noise all fixed, and vary only the AUROC level: reorder the observed scores
across compounds so the panel attains a target AUROC of 0.42, 0.50, 0.60, 0.70, 0.80 and 0.90,
then perturb with the measured noise (σ = 0.096 eV at exh 4) over many draws and record the
induced AUROC sd at each level. This isolates level from every other variable, which no
comparison between our real panels can do — they differ in *n* and level at once.

It is a simulation and is labelled as one. Its input noise is measured, not assumed.

**Prediction, fixed now:** the induced sd will be **largest near AUROC 0.5 and fall as the level
moves away**, because perturbation only flips pairs that are nearly tied and overlap is greatest
at chance. Therefore the floor reported in this paper, measured at 0.42, is **conservative** for
a better-performing benchmark.

**If the prediction is wrong** — if the floor is flat in level, or rises with it — the paper's
floor cannot be presented as transportable to benchmarks with real signal, and §3.2's scaling
result must be qualified as holding only near the level at which it was measured.

## What neither addresses

Neither touches the fact that Vina is below chance on this panel. That is a separate disclosure
the manuscript owes its reader regardless of these results, and it is being added independently.
