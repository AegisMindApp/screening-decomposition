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

---

## Amendment 2 — resource envelope changed to paid GCP, 3 Sep 2026

**Written before any GCP run. This changes a resource constraint, not a scientific rule.**

The 6.0 GPU-hour ceiling in the original pre-registration was explicitly a **Kaggle quota**
limit — "we do not spend the week's quota discovering that it did not fit". It fired correctly:
229.4 s/ligand projected 47.9 GPU-hours and the run abandoned itself.

The user has authorised paid GCP compute, so the binding constraint is now **money, not quota**.
Re-using a 6.0-hour ceiling in that setting would be arbitrary — it encodes a limit that no
longer exists.

**New ceiling: 20 GPU-hours**, fixed now, before the L4 probe runs. Justification: at L4
on-demand (~USD 0.71/hr ≈ AUD 1.10/hr) that is ≈ AUD 22, inside both the pre-existing AUD 30
monthly budget alert and the AUD 100 project budget created today. On spot it is ≈ AUD 9.

**What is explicitly NOT changed:** every scientific reading rule. The primary endpoint remains
`affinity_probability_binary`, the bar remains the descriptor baseline of 0.7653, the required
margin remains **> +0.04 AUROC with a CI excluding zero**, and the dropout, positive and
permutation controls are untouched. No threshold that bears on the *answer* has moved.

**Hardware substitution.** The user asked for A100 spot. This project has
`NVIDIA_A100_GPUS = 0`, `PREEMPTIBLE_NVIDIA_A100_GPUS = 0` and `A2_CPUS = 0` in us-central1,
us-east1 and europe-west4 — the A100 family cannot be created at all without a quota-increase
request under human review. `NVIDIA_L4_GPUS = 1`, so **L4 is the largest card actually
available** and is used instead: an estimated ~50 s/ligand against the A100's ~20, at comparable
total cost. The estimate is unverified — the same two-batch probe measures it before the full
panel is allowed to start.

---

## Amendment 3 — overage authorised by the budget owner, 3 Sep 2026

The L4 probe measured 101.0 s/ligand, projecting **21.07 GPU-hours** against the 20.0 ceiling
set in Amendment 2 — a 5% overshoot, about AUD 2.

**This line was not moved by me.** I reported the measurement and the three options (authorise
the overage, run a random stratified 700 inside the ceiling, or stop) and the budget owner chose
to authorise the full 751-compound panel. That distinction matters: a resource ceiling relaxed
by the person who owns the resource, after being shown the number, is a decision; the same
ceiling relaxed by me after seeing the number would be moving a goalpost.

Ceiling raised to **25.0 GPU-hours** — enough headroom for the measured 21.07 plus variance,
not an open budget. ~AUD 27 at L4 on-demand, inside the AUD 100 project budget.

**Unchanged, again:** every scientific reading rule. Primary endpoint
`affinity_probability_binary`; bar = descriptor baseline **0.7653**; required margin **> +0.04
AUROC with a CI excluding zero**; dropout, positive-control and permutation checks all intact.
The full 751-compound panel is retained precisely so Boltz-2 is scored on the identical
compounds and labels as Vina (0.4530) and the descriptors — option 2 would have broken that
comparability for AUD 2.
