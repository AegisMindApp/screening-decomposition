# What actually moves virtual screening performance: two measured resolution limits and six pre-registered interventions

**John Goodman** · OceanSparx Pty Ltd, Sydney, Australia · john.goodman@oceansparx.com
ORCID [0009-0008-4080-8043](https://orcid.org/0009-0008-4080-8043)

**PREPRINT — 7 September 2026**


## Abstract

Reported improvements to structure-based virtual screening are usually differences in AUROC or
enrichment between two protocols on one benchmark. We asked what a benchmark of typical size can
actually resolve, and then measured six interventions against that limit with the reading rule for
each committed to version control before the data existed.

Two limits, not one. Scoring the *same* Vina poses with a second implementation of the same
scoring function — agreeing at r = 0.970, no compound differing by more than 2 kcal/mol — moves
AUROC by 0.020 [95% CI 0.011, 0.029]. Re-running the whole protocol at a different search effort
and scoring it the same way moves AUROC by 0.039 [0.024, 0.056]. The first bounds a comparison of
scoring functions on fixed poses; the second bounds a comparison of protocols, which is what most
published improvements are. Both fall inside the range in which docking improvements are typically
reported.

On SARS-CoV-2 Mpro: changing the scoring function on identical poses moved AUROC +0.140,
repairing a receptor preparation defect +0.045, scoring the pose ensemble rather than the top
pose +0.055, supplying the correct binding site to a co-folding model −0.013, and eight-fold
search effort −0.019. Two clear the relevant floor, one is marginal and two fall inside. The
sixth — substituting DiffDock-L pose generation — is **not evaluable** on this panel: DiffDock
places a fifth of the compounds outside the binding site, and that fifth is heavier and more
often active than the panel, so any comparison on the remainder is biased. These are distributions over cross-validation fold
assignment, not point estimates: we previously published the pose-ensemble arm as +0.019 and
refuted, and it proved to be the 2nd percentile of a distribution centred at +0.031. The search-effort result is the sharpest either way: its measured
effect is smaller than the run-to-run reproducibility of the protocol it was measured in.

Docking adds nothing demonstrable over seven free physicochemical descriptors on either panel we
assembled: on Mpro it significantly subtracts (−0.002 [−0.004, −0.000]) and on Factor Xa its
+0.015 [−0.003, +0.034] crosses zero. Applying the same descriptors to LIT-PCBA gives a median
AUROC of 0.7260 across 15 targets against published AutoDock Vina at 0.581 and GNINA at
0.611–0.616. These descriptor baselines are supervised on each benchmark's own labels while
docking is not, so this is a bias control in the sense of Sieg et al., not a claim that
descriptors are a deployable substitute. On the AVE-debiased variant the descriptor median falls
to 0.6510.

One method exceeded the floor. Boltz-2, an open-weights co-folding model scoring compounds from
sequence and SMILES alone, does not beat the descriptor baseline on its own (Mpro +0.026, Factor
Xa +0.019, both intervals spanning zero) but adds +0.093 [+0.062, +0.125] and +0.075 [+0.051,
+0.099] respectively when combined with those descriptors. Its advantage is not that it is less
property-driven — descriptors explain 63% of its score and 63% of Vina's — but that on Mpro its
residual, the part orthogonal to those descriptors, ranks actives at 0.657 where six docking
scores residualised the same way sat between 0.494 and 0.573. Both descriptor comparisons are
supervised on each target's own labels; transferred cold to the other target the descriptor model
ranks actives *below chance*, and the standalone Boltz-2 score — which requires no labels at all —
is the only configuration usable on a target with no known actives. Across twelve cross-validation fold
seeds the two do not overlap: the highest docking residual observed is 0.581, the lowest Boltz-2
residual 0.644.

One point outside a band is not a generalisation, so we pre-registered a third scoring approach
to test it. FlashBind, an EGNN affinity head over learned docking poses, residualises to 0.538
[0.529, 0.546] — inside the docking band — with all three pre-registered controls passing and a
pocket-agreement check ruling out the confound that the band is a property of our docking box.
The claim that this ceiling belongs to the scoring approach rather than to these panels therefore
remains untested rather than supported.


## 1. Introduction

A virtual screening method is normally shown to work by reporting enrichment on a retrospective
benchmark. Two properties of that evidence are rarely reported alongside it. The first is what
difference the benchmark could distinguish from noise. The second is which component of a
multi-part protocol produced the difference — absent that, a pipeline change gets attributed to
whichever part its authors were working on.

Analogue and decoy bias in retrospective sets is well documented [1,2], and LIT-PCBA [3] was
constructed to control it — though LIT-PCBA has itself since been audited and found to carry
leakage and redundancy across its splits [8]. That literature establishes that reported enrichment
can measure property matching. It does not establish how much of a *given* reported improvement is
resolvable at all, which is the question here.

This paper reports a resolution limit, measured rather than assumed and in two variants; a set of
six single-factor comparisons against it, each with its reading rule fixed in advance; and a
descriptor baseline on LIT-PCBA, to test whether our own panels are unrepresentative.

## 2. Methods

**Benchmarks.** Two panels assembled from ChEMBL with measured actives and measured inactives:
SARS-CoV-2 Mpro (PDB 7VU6, 751 compounds, 34.2% active) and Factor Xa (886 compounds, 45.7%
active; an earlier 854-compound freeze of the same panel is used for the marginal-value
comparison in §3.4 and is labelled where quoted). LIT-PCBA [3] full and AVE-debiased variants, 15
targets, used as distributed.

**Descriptor baseline.** Molecular weight, clogP, hydrogen-bond donors, hydrogen-bond acceptors,
rotatable bonds, topological polar surface area, formal charge. Logistic regression, out-of-fold
predictions under 5-fold stratified cross-validation. On the AVE variant the authors' own
train/validation split is honoured rather than re-folded, so the AVE and full-set figures differ
in estimator as well as in compound selection. Inactives capped at 25,000 per LIT-PCBA target by
seeded random sampling.

**Docking.** AutoDock Vina 1.2.5, 22 Å box, receptor passing a redocking control at the screening
protocol. Exhaustiveness was 4 for the scoring-function, receptor-repair and pose-ensemble arms;
the search-effort arm compares exhaustiveness 4 against 32. gnina 1.3 CNN rescoring applied to the
identical Vina output poses. DiffDock-L for the pose-source arm, with poses minimised before
scoring. Supplementary Table S2 gives the receptor state, exhaustiveness and baseline for each arm.

**Boltz-2.** Version 2.2.1, protein sequence plus ligand SMILES, `affinity_probability_binary` as
the endpoint, fixed in advance. `diffusion_samples = 1` for cost. Pocket conditioning, where
applied, supplied as 1-based sequence indices derived from residues within 8 Å of the docking box
centre.

**Pre-registration.** Every reading rule — endpoint, threshold, control, abort condition — was
committed to a public git repository before the corresponding data existed. Supplementary Table S1
lists each arm's pre-registration file, commit hash and commit date. Amendments are appended rather than
edited, with the reason recorded. Two comparisons were not pre-registered and are identified as
such in S1.

**Statistics.** Paired bootstrap over compounds, 95% percentile intervals: 10,000 resamples for
the interventions in §3.2 and for the resolution floors, 4,000 for the Boltz-2 contrasts. A
difference is reported as demonstrated only when it exceeds the applicable resolution floor *and*
its interval excludes zero.

## 3. Results

### 3.1 Two resolution limits, 0.020 and 0.039 AUROC

gnina's `--score_only` mode re-evaluates the Vina scoring function on a pose file it is handed.
Applied to the exact poses Vina wrote at exhaustiveness 4, it reproduces the cached Vina score at
r = 0.9700, mean difference −0.108 kcal/mol, with no compound differing by more than 2 kcal/mol.
The two score vectors nonetheless give AUROC 0.3992 and 0.3791 — **0.020 apart [+0.011, +0.029]**.
Call this the **scoring floor**: what two implementations of one scoring function differ by when
placement is held fixed.

Comparing instead against the exhaustiveness-32 cache — same compounds, same scoring function, an
independent stochastic search — gives r = 0.9372 and AUROC 0.4186 against 0.3793, **0.039 apart
[+0.024, +0.056]**. Call this the **protocol floor**: run-to-run reproducibility of the pipeline
end to end.

Neither pair is a meaningfully different description of binding, yet the metric separates each by
an amount comparable to published improvements. We report a difference as unresolvable when it
falls below the floor that matches the comparison being made: the scoring floor where poses are
held fixed, the protocol floor where placement or search is allowed to vary.

### 3.2 One intervention clears its floor; four fall inside it

![**Five of the six interventions against the two measured resolution limits.** The
pose-source arm is omitted because it is not evaluable on this panel (§3.2); plotting a void
estimate would misrepresent it. Points are
ΔAUROC with 95% paired-bootstrap intervals; dashed lines are the two floors from §3.1 with their
own bootstrap intervals shaded. Each intervention is coloured by the floor that applies to it —
the scoring floor where poses were held fixed, the protocol floor where placement or search was
allowed to vary. The two floor intervals overlap between 0.024 and 0.029. Note that the
eight-fold search-effort comparison *is* the measurement that defines the protocol floor, so its
effect cannot exceed it by construction.](figures/figure1_floors.png)

| intervention | ΔAUROC | 95% CI | pre-registered bar | pre-registered verdict |
|---|---|---|---|---|
| Scoring function, identical poses (gnina CNN) | **+0.140** | [+0.091, +0.187] | CI entirely above +0.10 | directionally supported, **not at full strictness** (lower bound 0.091) |
| Receptor preparation repair | **+0.045** | [+0.024, +0.067] | CI entirely above 0 | **supported** |
| Pose ensemble rather than top pose | **+0.055** | [+0.043, +0.068] | > +0.04, CI excluding 0 | **supported** — corrects a previously published refutation, see below |
| Pose generator substituted (DiffDock-L) | — | — | CI above +0.05 would refute "scoring is the problem" | **not evaluable** — biased exclusion, see below |
| Binding site supplied to Boltz-2 | −0.013 | [−0.026, +0.000] | > +0.04, CI excluding 0 | **not demonstrated** |
| Eight-fold search effort | −0.019 | [−0.032, −0.006] | *not pre-registered* | — |

**On re-gating, and why we do not do it silently.** Four of these bars were fixed before the
resolution limits in §3.1 were separated, and three of them cite a single "≈0.04 measurement
floor" that we now know conflates the scoring and protocol cases. A pre-registration cannot be
moved after the fact, so the verdicts above are read against the bars as written. The corrected
floors are reported alongside as a secondary observation, and they change the picture in one place
only:

| intervention | poses fixed? | applicable floor | reading against it |
|---|---|---|---|
| Scoring function (gnina CNN) | yes | scoring, 0.020 | clears |
| Receptor preparation repair | receptor only; ligands re-docked | protocol, 0.039 | **marginal** — the point estimate exceeds the floor but the two intervals overlap across nearly their whole width |
| Pose ensemble | yes | scoring, 0.020 | **clears** — all 40 fold seeds above it, lowest +0.040 |
| Pose generator (DiffDock-L) | no | protocol, 0.039 | **not evaluable** |
| Binding site to Boltz-2 | n/a — different model | protocol, 0.039 | inside |
| Eight-fold search effort | no | protocol, 0.039 | inside |

The receptor-repair row is the one place the two readings diverge: it passes its pre-registered
bar cleanly and is only marginal against a floor fixed afterwards. We report both and claim
neither over the other.

These are six single-factor comparisons, not a partition of one pipeline. Four were run on the
donor-defective receptor described below, before it was repaired; the receptor-repair row measures
what that condition costs, and Supplementary Table S2 gives each arm's baseline in full. Their
nulls are therefore conditional on a receptor now known to be broken.

Two interventions clear their floor, and they do not divide as cleanly as *what is computed*
versus *where the ligand is placed*. Changing the scoring function clears by a wide margin; so
does scoring the pose ensemble rather than the top pose — a placement-side intervention, in that
it reads the distribution the search already produced. What does **not** matter is generating the
poses differently, or generating more of them. The narrower statement the data supports is that
search *effort* and pose *source* are inert, while how the pose distribution is *read* is not.

The search-effort row makes the point most economically: its CI
excludes zero, so the effect is real and negative, but it is smaller than the protocol floor — and
the exhaustiveness-4/32 pair *is* the protocol floor, so eight times the compute is by
construction indistinguishable from running the same protocol twice. **The pose-source arm is not evaluable on this panel, and we withdraw the number we published
for it.** DiffDock's top-ranked pose falls outside the 22 Å box for about a fifth of compounds
even under pocket conditioning, so the comparison can only be made on the remainder — and that
remainder is not a random sample. The 146 excluded compounds average 474.6 Da against 397.5 Da
for those retained, and are 44.5% active against the panel's 34.2%. Excluding them enriches the
surviving set in inactives and biases any ranking comparison built on it.

We previously reported this arm at +0.017 [−0.027, +0.064]. That figure was computed on the same
80% subset, produced by the same mechanism, and the exclusion was not tested for. It is withdrawn
rather than re-reported: the arm is uninterpretable on this panel on either receptor, which is a
property of DiffDock's behaviour on this target rather than of the analysis. What survives is the
descriptive observation that where both methods do place a ligand, their coordinates are nearly
uncorrelated (r = 0.139).

**The pose-ensemble arm carries a correction.** We previously published this arm at +0.019 with a
pre-registered verdict of *refuted*, measured on the donor-defective receptor. Re-run on the
repaired receptor across 40 cross-validation fold seeds it gives **+0.055 [+0.043, +0.068]**,
every seed above the floor. The earlier figure is not a contradictory result: on the defective
receptor the same pipeline gives +0.031 with a full seed range of [+0.017, +0.043], and the
published +0.019 sits at its 2nd percentile. That refutation was a single draw from a
distribution straddling the floor, and the honest verdict on that data was always *unresolved*.
Paired per seed, repairing the receptor contributes +0.024 [+0.011, +0.039] to this arm and helps
in 39 of 40 seeds.

**The receptor defect is worth reporting for its own sake.** Every receptor we had prepared
carried zero hydrogen-bond donors: the converter typed all 754 Mpro nitrogens as acceptors and
added no polar hydrogens, so the scoring function could not form a protein-donor hydrogen bond at
all. Repairing it moved AUROC +0.045 and changed six of the top ten shortlisted compounds. A
defect invisible in every summary statistic reordered most of the answer.

### 3.3 Seven descriptors reach 0.726 on LIT-PCBA, against published docking at 0.58–0.62

| | median AUROC | targets |
|---|---|---|
| descriptors, LIT-PCBA full | **0.7260** | 15 |
| descriptors, full, ≥50 actives | 0.6830 | 10 |
| descriptors, AVE-debiased | **0.6510** | 15 |
| AutoDock Vina, published [5] | 0.581 | 15 |
| GNINA, published [5] | 0.611–0.616 | 15 |

Per-target values are in Supplementary Table S4. Both published figures come from one source,
Sunseri & Koes [5], Table 2; we are not aware of an independent second measurement of Vina on this
benchmark. An independent 2026 benchmarking effort on LIT-PCBA [4] reports ROC-AUC "hovering near
0.6" for its strongest classical scorer, but its docking arm is AutoDock-GPU rather than Vina, so
it is not a substitute for the Vina figure. On the full set, seven descriptors costing microseconds
exceed published Vina by 0.145,
or 0.102 restricting to the ten targets with ≥50 actives. Our own panels' descriptor baselines
(Mpro 0.7654, Factor Xa 0.7041) fall inside LIT-PCBA's 0.577–0.864 range, so they are not
unrepresentative.

**What this comparison is and is not.** Our descriptor baselines are supervised on each target's
own actives and inactives; docking scores require no labels. The comparison is therefore a bias
control in the sense of Sieg et al. [2] — it bounds how much of a benchmark's apparent enrichment
is reachable from ligand properties alone — not a claim that descriptors are a deployable
substitute for docking on a novel target. Sunseri & Koes [5] reach the opposite conclusion, that
their CNN models significantly outperform simple-descriptor models; their descriptor models were
fit off-benchmark to PDBbind and CrossDocked affinity data and transferred, and the difference in
conclusion is a difference in what the baseline is allowed to see.

**What was not matched.** Our baselines use up to 25,000 inactives per target drawn by seeded
random sampling, against published figures computed on the full inactive sets and, in [5], taking
the maximum score over all reference receptors for the 13 of 15 targets with multiple templates.
AUROC is insensitive to negative-class size in expectation, so we expect no systematic direction
to the sampling difference, but this is not a matched comparison and we do not claim one.

AVE debiasing reduces the descriptor median to 0.6510, robust to excluding targets with few
validation actives (0.6510, 0.6497, 0.6510 at ≥0, ≥20, ≥50 actives). Two things change between
0.7260 and 0.6510 — the compound selection and, because the AVE variant ships its own
train/validation split, the estimator — so the 0.075 drop is attributable to the AVE variant as a
whole rather than to debiasing alone. What survives is that 0.651 from seven trivial properties on
a deliberately debiased set is well clear of chance.

**A comparison we cannot make.** No published Vina number exists for the AVE-debiased set.
Comparing our AVE descriptor result against a full-set docking result would not be like for like
and we do not do so. Vina would presumably also fall under debiasing; by how much is unknown to us.

### 3.4 A co-folding model exceeds the floor in combination — but only the standalone score is
deployable, and a third approach does not exceed it at all

| | Mpro | Factor Xa |
|---|---|---|
| descriptors | 0.7654 | 0.7041 |
| Boltz-2 alone | 0.7913 | 0.7227 |
| descriptors + Boltz-2 | **0.8588** | **0.7791** |
| Boltz-2 − descriptors | +0.026 [−0.023, +0.074] | +0.019 [−0.030, +0.066] |
| combination − descriptors | **+0.093 [+0.062, +0.125]** | **+0.075 [+0.051, +0.099]** |

Both halves replicate across a viral protease and a coagulation factor with independently
assembled panels: the standalone comparison fails on both, the combination clears on both with
overlapping intervals.

**Both rows are fitted on each target's own labels.** The combination in the table above is an
out-of-fold logistic regression over the target's actives and inactives, and so carries exactly
the supervision caveat §3.3 attaches to the descriptor baseline. On a target with no known
actives there is nothing to fit it with. We tested whether the gain survives when the weights
come from a different target — fitting on one panel in full and applying the model cold to the
other, in both directions.

Boltz-2 is not less property-driven than docking — descriptors explain 63.1% of its score and
62.8% of Vina's. The difference is the residual. Regressing the seven descriptors out of each
score out-of-fold and ranking actives on the remainder:

| score | raw AUROC | residual AUROC | 12-seed range |
|---|---|---|---|
| Vina, repaired receptor (Mpro) | 0.4530 | 0.5187 | [0.5131, 0.5259] |
| Vina, raw receptor (Mpro) | 0.4079 | 0.5092 | [0.4980, 0.5222] |
| gnina Vina term (Mpro) | 0.3791 | 0.4941 | [0.4887, 0.5017] |
| gnina CNNaffinity (Mpro) | 0.5389 | 0.4994 | [0.4908, 0.5091] |
| gnina CNNscore (Mpro) | 0.5041 | 0.5426 | [0.5320, 0.5592] |
| Vina, repaired receptor (Factor Xa) | 0.6775 | 0.5733 | [0.5634, 0.5810] |
| **Boltz-2 `prob_binary` (Mpro)** | 0.7913 | **0.6567** | **[0.6439, 0.6648]** |
| FlashBind `binary` (Mpro) | 0.4185 | 0.5384 | [0.5287, 0.5456] |

Six docking scores fall between 0.494 and 0.573; Boltz-2 sits at 0.657. The residual estimate
moves with the cross-validation fold assignment, so each value is a mean over twelve fold seeds
with its full range given, and the separation holds at the worst case — the highest docking
residual observed across all seeds is 0.581, the lowest Boltz-2 residual 0.644. It also survives
the choice of estimator: under linear rather than gradient-boosted residualisation the docking
band is 0.418–0.614 and Boltz-2 is 0.748. The residual analysis is single-target — we did not
compute a Factor Xa Boltz-2 residual.

**The increment transfers; the model does not.** Cold, the Boltz-2 increment clears the floor in
both directions (+0.088 [+0.056, +0.119] Mpro → Factor Xa; +0.248 [+0.210, +0.287] the other
way), but the transferred *models* are largely unusable. Sorted by what a screener actually has:

| configuration | Mpro | Factor Xa | labels required |
|---|---|---|---|
| **Boltz-2 alone, rank by score** | **0.7913** | **0.7227** | **none** |
| descriptors, in-target out-of-fold | 0.7673 | 0.6994 | the target's own actives |
| descriptors + Boltz-2, in-target out-of-fold | 0.8647 | 0.7803 | the target's own actives |
| descriptors, cold from the other target | 0.4553 | 0.4011 | another target's actives |
| descriptors + Boltz-2, cold from the other target | 0.7036 | 0.4890 | another target's actives |

Two things follow. **Seven descriptors fit on one target and applied to another rank actives
below chance in both directions** — inverted transfer, not weak transfer. The property signature
separating actives on one panel is close to the opposite of the other's, which is the strongest
evidence here that the descriptor baseline measures target-specific bias rather than chemistry
that generalises, and is why §3.3 treats it as a bias control.

**And on a novel target the only usable configuration is the standalone score**, which needs no
labels, no receptor preparation and no docking box. That is the configuration reported as *not
demonstrated* at the top of this section. Both readings are correct and they answer different
questions: the in-target comparison asks whether Boltz-2 beats a baseline fitted on the answers;
the prospective question asks what can be run when there are no answers to fit to. Under the
first it fails; under the second it is the only candidate that works.

Residualising against ligand descriptors is a weaker test than receptor ablation. Binding free
energy genuinely covaries with size and lipophilicity, so a physically correct scorer would also
lose signal under this operation — that is why ligand efficiency exists. The defensible claim is
narrow: the part of these docking scores orthogonal to seven descriptors does not rank actives in
these sets, and on these two targets the ceiling appears to be a property of the scoring approach
rather than of the task.

**We then tested that generalisation, and it did not hold.** A third scoring approach clearing
the ceiling would have been the evidence for it, so we pre-registered one against FlashBind
(EGNN affinity head over FABind+ poses), chosen because it sits architecturally between the two
families and claims early enrichment competitive with co-folding at roughly fiftyfold lower cost.
Reading rules, controls and abort conditions were committed before any score existed. All three
controls pass: the descriptor baseline reproduces the panel at 0.7664 against a required
0.7654 ± 0.005, coverage is 750 of 751, and shuffled labels give 0.4978.

Its residual is **0.5384 [0.5287, 0.5456]** — inside the docking band, with the whole twelve-seed
range inside it, and inside the narrower Mpro-only band of 0.494–0.543 as well. Its increment
over the descriptors is +0.0247 [+0.0217, +0.0271] across the same twelve seeds, positive but
clear below the lower end of Boltz-2's +0.093 [+0.062, +0.125]. Because FABind+ selects its own
pocket rather than receiving our box, we pre-registered a pocket-agreement check to separate
"property of the scoring approach" from "property of being box-constrained": its pocket centroids
sit a median 8.00 Å from our box centre with 749 of 750 inside the box, so that alternative does
not apply.

Boltz-2 therefore remains the only point outside the band, and one point is not a generalisation.
The claim that the ceiling is a property of the scoring approach rather than of this panel stands
untested rather than supported.

FlashBind's raw AUROC on this panel is 0.4185, below chance — but so are three of the five
docking scores in the table above, and we draw no conclusion about the method from that. Its
binary head is trained on MF-PCBA high-throughput data, where activity is defined differently
from the ChEMBL potency threshold used here, and no protease assay appears anywhere in that
training set. This panel is out of domain for it in both target and label definition.

## 4. Limitations

**Two self-assembled targets** for the six comparisons, and one of them (Mpro) is a target where
our docking protocol scores below chance. LIT-PCBA broadens the descriptor result to 15 targets
but not the comparisons, which remain a single-benchmark measurement.

**Two of the six arms sit on the donor-defective receptor**, as set out in §3.2 and
Supplementary Table S2 — the search-effort comparison and the scoring-function arm. The
pose-ensemble arm was re-run on the repaired receptor and its verdict changed. The pose-source
arm was also re-run and proved **not evaluable** on this panel in either condition.

**One of the six interventions cannot be measured here at all.** That is a limitation of the
benchmark-plus-method combination rather than a null result, and we report it as such rather than
as evidence that pose source does not matter.

**Four of the six pre-registered bars predate the floor correction.** Three of them cite a single
"≈0.04 measurement floor" that §3.1 shows was conflating two quantities. The bars stand as written
and the verdicts in §3.2 are read against them; the corrected floors are reported beside those
verdicts, never substituted for them.

**Cross-validation fold assignment moves AUROC differences on these panels by 0.01–0.03**, which
is the same size as the effects being measured. Two of our own results were misled by single-seed
estimates: the residual band in §3.4, and the pose-ensemble arm, where a published refutation
proved to be a 2nd-percentile draw. Every difference of this magnitude is reported here as a
distribution over fold seeds with its range, and we would treat any single-seed AUROC difference
below about 0.05 on a benchmark of this size as uninterpretable.

**One analysis script was never committed.** The original pose-ensemble implementation exists
only as its pre-registration and its write-up, so its published number could not be re-derived
and had to be reconstructed from a prose description. Every script behind the results reported
here is in the repository with its outputs.

**The resolution floors are measured on one benchmark of one size**, on the unprotonated receptor.
They should be re-measured rather than assumed elsewhere; the method costs one extra scoring pass.
An earlier version of this analysis reported a single 0.039 floor and described it as a comparison
of identical poses; it is not, and the two variants are separated here.

**The descriptor baselines are supervised and docking is not**, as stated in §3.3. This is the
limitation most likely to be mistaken for a stronger claim than we make.

**No Vina number on AVE-debiased LIT-PCBA**, and only one source for the full-set figure.

**The third-approach test is one method on one target.** FlashBind was pre-registered as the
test of whether the residual ceiling generalises, and it does not clear the ceiling — but a panel
whose actives skew large is unfavourable to a model whose score is strongly anti-correlated with
size, and its binary head was trained on a different activity definition with no protease assay
in its training set. A negative from one out-of-domain target leaves the generalisation open; it
does not close it in the other direction either.

**The Boltz-2 evaluation is contested territory.** A 2025 and a 2026 paper disagree about its
virtual screening performance [6,7]. Our contribution is the pre-registration and the residual
decomposition, not the discovery that it works.

**A caveat we withdrew.** Boltz-2's first run was blind while its docking comparator was given the
binding site, and we recorded that as a handicap making +0.093 a lower bound. Supplying the pocket
produced 0.7783 against the blind run's 0.7913, with the two runs correlating at r = 0.952. The
caveat was unnecessary and is withdrawn. It is reported here because a pre-registered caveat that
turns out to be wrong is part of the record.

## 5. Data and code availability

All pre-registration files, amendments, analysis code and raw outputs are public at
**https://github.com/AegisMindApp/screening-decomposition**, with the commit hash and date for
each pre-registration given in Supplementary Table S1. The FlashBind arm of §3.4 is the
exception at the time of writing: its pre-registration, amendments, worker, controls and raw
scores exist in the origin repository at the hashes given in S1 and are marked there as pending
export. The history is intact rather than a
snapshot, so any cited commit can be resolved and dated directly
(`git log -1 --format=%ad --date=iso <commit>`).

The dates are git metadata carried through from the repository this was exported from; they are
not a third-party timestamp, and a reader can confirm ordering and internal consistency but
cannot from this alone exclude authored metadata. The resolution-floor analysis is reproducible
from that repository with no dependencies beyond the Python standard library, and asserts a
positive control against previously published values before reporting; that control failed once,
on 6 September 2026, and the defect it caught is recorded in the repository alongside the result.

Two comparisons — the search-effort arm and the LIT-PCBA descriptor baseline — were not
pre-registered and are identified as such in S1.

## References

1. Chen L, Cruz A, Ramsey S, Dickson CJ, et al. (2019) Hidden bias in the DUD-E dataset leads to misleading performance of deep learning in structure-based virtual screening. *PLOS ONE* 14:e0220113. doi:10.1371/journal.pone.0220113

2. Sieg J, Flachsenberg F, Rarey M (2019) In Need of Bias Control: Evaluating Chemical Data for Machine Learning in Structure-Based Virtual Screening. *J Chem Inf Model* 59:947–961. doi:10.1021/acs.jcim.8b00712

3. Tran-Nguyen VK, Jacquemard C, Rognan D (2020) LIT-PCBA: An Unbiased Data Set for Machine Learning and Virtual Screening. *J Chem Inf Model* 60:4263–4273. doi:10.1021/acs.jcim.0c00155

4. Abo-Dahab Y, Xiang X, Chun J, Zhao L (2026) Benchmarking Single-Pose Docking, Consensus Rescoring, and Supervised ML on the LIT-PCBA Library: A Critical Evaluation of DiffDock, AutoDock-GPU, GNINA, and DiffDock-NMDN. arXiv:2605.01681

5. Sunseri J, Koes DR (2021) Virtual Screening with Gnina 1.0. *Molecules* 26(23):7369. doi:10.3390/molecules26237369

6. Furui K, Ohue M (2025) Boltzina: Efficient and Accurate Virtual Screening via Docking-Guided Binding Prediction with Boltz-2. arXiv:2508.17555

7. Wan S, Zhang X, Xue X, Coveney PV (2026) On the Reliability of AI Methods in Drug Discovery: Evaluation of Boltz-2 for Structure and Binding Affinity Prediction. arXiv:2603.05532

8. Huang A, Knight IS, Naprienko S (2025) Data Leakage and Redundancy in the LIT-PCBA Benchmark. arXiv:2507.21404

\newpage

# Supplementary material

**S1. Pre-registration index.** Every reading rule — endpoint, threshold, control, abort
condition — was committed before the corresponding data existed. Each row gives the file, the
commit that introduced it, and that commit's date. Amendments are appended to the file rather
than editing the original text, with the reason recorded in place.

Paths are relative to `analysis/`; the file is `PREREGISTRATION.md` except where named.

| arm | directory | public | origin | committed |
|----------------------|------------------------------------|------------|------------|------------------|
| Scoring function (gnina rescore) | `three_arm_docking/` | `f5e42306a` | `cd8f55639` | 2026-08-31 08:46 |
| Receptor preparation repair | `receptor_prep/` | `05d9f78ec` | `1f68b20c4` | 2026-08-31 14:26 |
| Pose ensemble vs top pose | `pose_ensemble/` | `d5d08c1ac` | `448d954e5` | 2026-09-02 22:52 |
| Pose generator (DiffDock-L) | `three_arm_docking/` | `f5e42306a` | `cd8f55639` | 2026-08-31 08:46 |
| Boltz-2, Mpro, arm 1 (blind) | `boltz2/` | `8df29ba21` | `9e3083d02` | 2026-09-03 08:27 |
| Boltz-2 analysis script, written before the data | `boltz2/analyse_boltz2.py` | `42f94969b` | `5afb5bec8` | 2026-09-03 13:51 |
| Boltz-2, Mpro, arm 2 (pocket-conditioned) | `boltz2/PREREGISTRATION_ARM2.md` | `064969caf` | `e96d4fea1` | 2026-09-04 11:28 |
| Boltz-2, Factor Xa | `boltz2/PREREGISTRATION_TARGETS.md` | `00baa6db9` | `68bd9b94a` | 2026-09-04 11:40 |
| FlashBind, Mpro (third scoring approach) | `flashbind/` | pending export | `c6b4504f9` | 2026-09-08 10:06 |
| — amendment: pocket-agreement control | `flashbind/` | pending export | `24bc26573` | 2026-09-08 11:36 |
| — amendment: abort decision on torchdrug | `flashbind/` | pending export | `453e2031f` | 2026-09-08 12:04 |

The public repository was produced from the origin repository with `git filter-repo`, which
rewrites commit hashes but preserves author and committer dates. Both hashes are given so a
reader can check either. Verify a date with
`git log -1 --format=%ad --date=iso <commit>`.

Commit dates above are git metadata carried through from the project's own history; they are not
a third-party timestamp, and a reader should treat them as the authors' record rather than as
independent verification of ordering.

The three FlashBind rows are one pre-registration and two amendments appended to it before any
score existed. The first amendment adds the pocket-agreement control, because FABind+ selects its
own binding site and a residual above the band would otherwise have been confounded with not
being box-constrained. The second records, while the environment work was still failing, that the
arm would be reported blocked rather than assembled from workarounds if a required dependency
could not run — a condition that was then tested and not met, which is also recorded in place.

Bars citing a "≈0.04 measurement floor" — the pose-ensemble arm and both Boltz-2 arms — were fixed
before the floor correction described in §3.1 and S3. They are read as written; see §3.2.

The search-effort comparison (exhaustiveness 4 vs 32) and the LIT-PCBA descriptor baseline were
not separately pre-registered; both are reported here as such. The search-effort arm reuses two
runs that pre-dated the pre-registration programme, and the LIT-PCBA baseline applies a method
already fixed for the Mpro and Factor Xa panels to a third-party dataset.

**S2. Per-arm experimental conditions.** The six interventions are single-factor comparisons
against differing baselines, not a partition of one pipeline. Four were run before the receptor
defect described in §3.2 was found and repaired.

| intervention | baseline AUROC | receptor | exhaustiveness | n | poses held fixed? |
|---|---|---|---|---|---|
| Scoring function (gnina CNN) | 0.3992 | unprotonated | 4 | 743 | **yes** |
| Receptor repair | 0.4079 | unprotonated → protonated | 4 | 751 | receptor coordinates fixed; ligands re-docked |
| Pose ensemble | — | unprotonated | 4 | 751 | **yes** (same run, 5 modes) |
| Pose generator (DiffDock-L) | — | unprotonated | 4 (Vina arm) | 751 | no |
| Pocket conditioning (Boltz-2) | 0.7913 | none — co-folding from sequence | n/a | 750 | n/a |
| Eight-fold search effort | 0.4079 (exh 4) | unprotonated | 4 → 32 | 753 | no |

**S3. Resolution floor, both variants.** `analysis/docking_value/floor_interval.py`,
`FLOOR_INTERVAL.json`. Paired bootstrap over compounds, 10,000 resamples, unprotonated receptor.

| | value | 95% CI | n | r between the two score vectors |
|---|---|---|---|---|
| Protocol floor — exhaustiveness-32 Vina cache vs gnina `--score_only` on the exhaustiveness-4 poses | 0.0393 | [+0.0235, +0.0557] | 745 | 0.9372 |
| Scoring floor — exhaustiveness-4 Vina cache vs gnina `--score_only` on the same exhaustiveness-4 poses | 0.0201 | [+0.0112, +0.0286] | 743 | 0.9700 |

Neither score vector contains a compound differing from its counterpart by more than 2 kcal/mol.
An earlier description of the protocol floor as a comparison of "identical poses" was incorrect
and is corrected in `analysis/docking_value/FLOOR_CORRECTION.md`.

**S4. LIT-PCBA per-target results.** `analysis/litpcba/descriptor_baseline_full.json` and
`descriptor_baseline_ave.json`.

| target | full: actives | inactives | AUROC | AVE: actives | inactives | AUROC |
|---|---|---|---|---|---|---|
| ADRB2 | 17 | 25,000 | 0.7008 | 4 | 25,000 | 0.6016 |
| ALDH1 | 7,168 | 25,000 | 0.6169 | 1,343 | 25,000 | 0.6123 |
| ESR1_ago | 13 | 5,583 | 0.7260 | 3 | 908 | 0.8737 |
| ESR1_ant | 102 | 4,948 | 0.7525 | 25 | 794 | 0.6497 |
| FEN1 | 369 | 25,000 | 0.7596 | 91 | 25,000 | 0.7471 |
| GBA | 166 | 25,000 | 0.7456 | 41 | 25,000 | 0.7390 |
| IDH1 | 39 | 25,000 | 0.7921 | 9 | 25,000 | 0.7838 |
| KAT2A | 194 | 25,000 | 0.6656 | 48 | 25,000 | 0.6032 |
| MAPK1 | 308 | 25,000 | 0.7285 | 77 | 15,250 | 0.7129 |
| MTORC1 | 97 | 25,000 | 0.5945 | 24 | 8,243 | 0.5686 |
| OPRK1 | 24 | 25,000 | 0.8640 | 6 | 25,000 | 0.7672 |
| PKM2 | 546 | 25,000 | 0.6743 | 136 | 25,000 | 0.6068 |
| PPARG | 27 | 5,211 | 0.8459 | 6 | 861 | 0.7195 |
| TP53 | 79 | 4,168 | 0.5769 | 16 | 795 | 0.5050 |
| VDR | 884 | 25,000 | 0.6918 | 165 | 25,000 | 0.6510 |
| **median** | | | **0.7260** | | | **0.6510** |

**S5. Panels.** Mpro: 751 compounds, 257 active (34.2%), PDB 7VU6. Factor Xa: 886 compounds, 405
active (45.7%). The Factor Xa marginal-value comparison in §3.4 uses an earlier 854-compound
freeze of the same panel; where a number from that freeze is quoted it is labelled.
