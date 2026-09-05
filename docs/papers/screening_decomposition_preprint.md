# What actually moves virtual screening performance: two measured resolution limits and a pre-registered decomposition of six interventions

**John Goodman** · OceanSparx Pty Ltd, Sydney, Australia · john.goodman@oceansparx.com

**PREPRINT — 6 September 2026**

---

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
repairing a receptor preparation defect +0.045, scoring the pose ensemble rather than the top pose
+0.019, substituting DiffDock-L pose generation +0.017, supplying the correct binding site to a
co-folding model −0.013, and eight-fold search effort −0.019. Read against each arm's own
pre-registered bar, one is supported at full strictness, one directionally, and the rest are not.
Read against the resolution limits measured here — which were fixed after four of the six bars,
and so are reported alongside those readings rather than in place of them — one clears, one is
marginal, and four fall inside. The search-effort result is the sharpest either way: its measured
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
residual, the part orthogonal to those descriptors, ranks actives at 0.642 where five docking
scores we measured sat between 0.507 and 0.559.

---

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
committed to a git repository before the corresponding data existed. Supplementary Table S1 lists
each arm's pre-registration file, commit hash and commit date. Amendments are appended rather than
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

| intervention | ΔAUROC | 95% CI | pre-registered bar | pre-registered verdict |
|---|---|---|---|---|
| Scoring function, identical poses (gnina CNN) | **+0.140** | [+0.091, +0.187] | CI entirely above +0.10 | directionally supported, **not at full strictness** (lower bound 0.091) |
| Receptor preparation repair | **+0.045** | [+0.024, +0.067] | CI entirely above 0 | **supported** |
| Pose ensemble rather than top pose | +0.019 | [−0.017, +0.055] | > +0.04, CI excluding 0 | **refuted** |
| Pose generator substituted (DiffDock-L) | +0.017 | [−0.027, +0.064] | CI above +0.05 would refute "scoring is the problem" | scoring hypothesis **retained** |
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
| Pose ensemble | yes | scoring, 0.020 | inside |
| Pose generator (DiffDock-L) | no | protocol, 0.039 | inside |
| Binding site to Boltz-2 | n/a — different model | protocol, 0.039 | inside |
| Eight-fold search effort | no | protocol, 0.039 | inside |

The receptor-repair row is the one place the two readings diverge: it passes its pre-registered
bar cleanly and is only marginal against a floor fixed afterwards. We report both and claim
neither over the other.

These are six single-factor comparisons, not a partition of one pipeline. Four were run on the
donor-defective receptor described below, before it was repaired; the receptor-repair row measures
what that condition costs, and Supplementary Table S2 gives each arm's baseline in full. Their
nulls are therefore conditional on a receptor now known to be broken.

The one intervention that clears its floor outright concerns *what is computed*, not *where or how
thoroughly the ligand is placed*. The search-effort row makes the point most economically: its CI
excludes zero, so the effect is real and negative, but it is smaller than the protocol floor — and
the exhaustiveness-4/32 pair *is* the protocol floor, so eight times the compute is by
construction indistinguishable from running the same protocol twice. The pose-source result is
next sharpest: Vina and DiffDock poses correlate at r = 0.139, nearly uncorrelated coordinates,
and produce rankings differing by less than the floor.

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

### 3.4 A co-folding model exceeds the floor, but only in combination

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

Boltz-2 is not less property-driven than docking — descriptors explain 63.1% of its score and
62.8% of Vina's. The difference is the residual. After regressing out the seven descriptors
out-of-fold, Boltz-2's remainder ranks actives at **0.6418** on Mpro; Vina's at 0.5382. Across
five docking scores we residualised — Vina raw and repaired, two gnina heads, and Vina on Factor
Xa — the remainder fell between 0.507 and 0.559. Two residualisation procedures were used across
those five (out-of-fold logistic and out-of-fold gradient-boosted), so the band is indicative
rather than a like-for-like interval; Boltz-2's 0.6418 sits outside it under either. The residual
analysis is single-target: we did not compute a Factor Xa Boltz-2 residual.

Residualising against ligand descriptors is a weaker test than receptor ablation. Binding free
energy genuinely covaries with size and lipophilicity, so a physically correct scorer would also
lose signal under this operation — that is why ligand efficiency exists. The defensible claim is
narrow: the part of these docking scores orthogonal to seven descriptors does not rank actives in
these sets, and on these two targets the ceiling appears to be a property of the scoring approach
rather than of the task. A third scoring approach clearing it would be the test that generalises
this.

## 4. Limitations

**Two self-assembled targets** for the six comparisons, and one of them (Mpro) is a target where
our docking protocol scores below chance. LIT-PCBA broadens the descriptor result to 15 targets
but not the comparisons, which remain a single-benchmark measurement.

**Four of the six arms sit on the donor-defective receptor**, as set out in §3.2 and Supplementary
Table S2. Whether the pose-ensemble and DiffDock nulls survive on a repaired receptor is untested.

**Four of the six pre-registered bars predate the floor correction.** Three of them cite a single
"≈0.04 measurement floor" that §3.1 shows was conflating two quantities. The bars stand as written
and the verdicts in §3.2 are read against them; the corrected floors are reported beside those
verdicts, never substituted for them.

**The resolution floors are measured on one benchmark of one size**, on the unprotonated receptor.
They should be re-measured rather than assumed elsewhere; the method costs one extra scoring pass.
An earlier version of this analysis reported a single 0.039 floor and described it as a comparison
of identical poses; it is not, and the two variants are separated here.

**The descriptor baselines are supervised and docking is not**, as stated in §3.3. This is the
limitation most likely to be mistaken for a stronger claim than we make.

**No Vina number on AVE-debiased LIT-PCBA**, and only one source for the full-set figure.

**The Boltz-2 evaluation is contested territory.** A 2025 and a 2026 paper disagree about its
virtual screening performance [6,7]. Our contribution is the pre-registration and the residual
decomposition, not the discovery that it works.

**A caveat we withdrew.** Boltz-2's first run was blind while its docking comparator was given the
binding site, and we recorded that as a handicap making +0.093 a lower bound. Supplying the pocket
produced 0.7783 against the blind run's 0.7913, with the two runs correlating at r = 0.952. The
caveat was unnecessary and is withdrawn. It is reported here because a pre-registered caveat that
turns out to be wrong is part of the record.

## 5. Data and code availability

All pre-registration files, amendments, analysis code and raw outputs are in the project
repository, with the commit hash and date for each pre-registration given in Supplementary Table
S1. At the time of writing that repository is private, so a reader cannot independently confirm
that a given commit predates the corresponding data; the hashes and dates in S1 are our assertion,
and the repository will be made public with its history intact. Readers should weight the
pre-registration claim accordingly until then.

## References

[1] Chen L, Cruz A, Ramsey S, Dickson CJ, et al. (2019) Hidden bias in the DUD-E dataset leads to
misleading performance of deep learning in structure-based virtual screening. *PLOS ONE*
14:e0220113. doi:10.1371/journal.pone.0220113
[2] Sieg J, Flachsenberg F, Rarey M (2019) In Need of Bias Control: Evaluating Chemical Data for
Machine Learning in Structure-Based Virtual Screening. *J Chem Inf Model* 59:947–961.
doi:10.1021/acs.jcim.8b00712
[3] Tran-Nguyen VK, Jacquemard C, Rognan D (2020) LIT-PCBA: An Unbiased Data Set for Machine
Learning and Virtual Screening. *J Chem Inf Model* 60:4263–4273. doi:10.1021/acs.jcim.0c00155
[4] Abo-Dahab Y, Xiang X, Chun J, Zhao L (2026) Benchmarking Single-Pose Docking, Consensus
Rescoring, and Supervised ML on the LIT-PCBA Library: A Critical Evaluation of DiffDock,
AutoDock-GPU, GNINA, and DiffDock-NMDN. arXiv:2605.01681
[5] Sunseri J, Koes DR (2021) Virtual Screening with Gnina 1.0. *Molecules* 26(23):7369.
doi:10.3390/molecules26237369
[6] Furui K, Ohue M (2025) Boltzina: Efficient and Accurate Virtual Screening via Docking-Guided
Binding Prediction with Boltz-2. arXiv:2508.17555
[7] Wan S, Zhang X, Xue X, Coveney PV (2026) On the Reliability of AI Methods in Drug Discovery:
Evaluation of Boltz-2 for Structure and Binding Affinity Prediction. arXiv:2603.05532
[8] Huang A, Knight IS, Naprienko S (2025) Data Leakage and Redundancy in the LIT-PCBA Benchmark.
arXiv:2507.21404
