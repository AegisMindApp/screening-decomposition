# Pre-registration — Boltz-2 on the Mpro benchmark

**Written 3 September 2026, committed before any Boltz-2 output exists.**

## Why this test

Two 2026 papers disagree. Boltzina (arXiv 2508.17555) reports Boltz-2 discriminating true from
false positives on ultralarge-screening hits *"overcoming by a large margin all scoring
functions tested so far on raw docking poses."* A reliability evaluation (arXiv 2603.05532)
reports it failing to predict affinity or binding probability correctly for selected actives.

Unlike every candidate improvement we have tested this session, the claimed effect is large
enough to sit **above our ~0.04 measurement floor**, so this benchmark can actually resolve it.

Boltz-2 predicts structure and affinity from **protein sequence + ligand SMILES**. It therefore
never touches our receptor preparation, and is immune to the donor-typing defect that affected
every docking number we have produced. That makes it a genuinely independent read.

## The bar is the descriptor baseline, not Vina

On these 751 compounds: Vina (repaired receptor) **0.4530**, seven free physicochemical
descriptors **0.7653**. Beating Vina is not interesting — we already know Vina is below chance
here. The question is whether Boltz-2 beats descriptors that cost microseconds.

## Primary output — fixed now to prevent post-hoc selection

Boltz-2 emits both `affinity_probability_binary` and `affinity_pred_value`. Labels here are
binary active/inactive, so **`affinity_probability_binary` is the primary endpoint**.
`affinity_pred_value` is reported as a secondary, clearly labelled. Choosing whichever scores
better after the fact is forbidden.

## Reading rules

1. **Primary.** ΔAUROC(Boltz-2 − descriptors) **> +0.04** with a paired-bootstrap 95% CI
   excluding zero → Boltz-2 is the first method we have tested that adds real value over free
   descriptors. Anything less is reported as **not demonstrated**, not as "promising".
2. **Secondary.** ΔAUROC(Boltz-2 − Vina) > +0.04 with CI excluding zero → replicates the
   literature's weaker claim that it beats docking scores.
3. **Combination.** Descriptors + Boltz-2 versus descriptors alone, same rule. This is the
   marginal-value test applied to a new model class.

## Controls — each can void the result

- **Dropout.** If Boltz-2 fails on some compounds, the active fraction among failures must be
  within 10 percentage points of the active fraction among successes. Outside that, the run is
  **class-selected and void**. This has already caught us twice: arm 2's 149 size-selected
  dropouts and PD-L1's 38 missing compounds at −14.3pp.
- **Positive control.** The pipeline must reproduce the known descriptor baseline of 0.7653
  ± 0.005 on the same compounds *before* any Boltz-2 number is read. If it cannot reproduce a
  number we already know, no null from it is trustworthy.
- **Permutation.** Shuffled labels must give AUROC within [0.45, 0.55]. Outside that, the
  evaluation is leaking and both readings are void.
- **No silent failures.** Every per-compound error is recorded with its message. A compound
  that produces no score must appear in the failures file, never be silently absent. A bare
  `try/except` that swallows a parse error is what made an earlier analysis report zero scored
  compounds without raising.

## Compute — abort condition, fixed in advance

The paper reports **~20 s/ligand**, which would be 4.17 GPU-hours for 751 compounds. That
figure is on the authors' hardware; Kaggle provides T4 or P100, which may be several times
slower. We have **under 8 GPU-hours** of quota, resetting Saturday.

**Stage 1 is a timing probe on 12 compounds.** It measures actual seconds/ligand on the
allocated card and extrapolates. **If the projected full-panel cost exceeds 6.0 GPU-hours, the
run is abandoned and reported as not attempted** — we do not spend the week's quota discovering
that it did not fit. Stage 2 runs the full panel only if the probe clears that bar.

## What would make me abandon the direction

Rule 1 failing with a clean dropout control. That would extend the descriptor result from
classical scoring functions to co-folding affinity models, which is a stronger and more
general finding than anything the docking arms produced — and it would close this avenue too.

## Scope

One target, one benchmark, one model. A negative here is evidence about Boltz-2 **on this
panel**, not evidence that Boltz-2 is a poor method. Its own papers evaluate it on tasks this
benchmark does not pose.

---

## Addendum, same day, still before any output

**Construct: the dimer.** Our 7VU6 receptor contains chains A (298 aa, residues 3–300) and B
(299 aa, 3–301), identical in sequence. Mpro is catalytically active only as a homodimer and
the S1 subsite of each protomer is shaped by the partner's N-terminal residues, so the monomer
is the wrong construct biologically *and* would not match the receptor we docked against.
Both chains are therefore supplied to Boltz-2.

This roughly doubles the token count, and the pairwise stack does not scale linearly, so the
compute risk is higher than the 4.17 GPU-hour estimate above. That is precisely what the
Stage-1 probe measures, and the 6.0 GPU-hour abort threshold is unchanged. Recording the choice
here so it cannot be revisited after seeing a cost.

**MSA is computed once.** The protein is identical across all 751 predictions, so the MSA is
generated a single time and reused for every ligand. Per-ligand MSA generation would dominate
runtime and is not part of what we are measuring.

**Provenance note.** This file was first committed at 3 Sep 2026 08:25 as `43ae0af8` in the
`solver-press` repository, written there by mistake because the shell working directory had not
been reset after a website push. It was moved here intact and removed from that repo in
`2bf4bae3`. Both commits predate any Boltz-2 execution, so the pre-registration guarantee is
unaffected and independently checkable in either history.
