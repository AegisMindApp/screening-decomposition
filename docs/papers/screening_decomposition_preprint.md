# What actually moves virtual screening performance: a pre-registered decomposition, a measured resolution limit, and seven descriptors that beat docking on LIT-PCBA

**John Goodman** · OceanSparx Pty Ltd, Sydney, Australia · john.goodman@oceansparx.com

**PREPRINT — 6 September 2026**

---

## Abstract

Reported improvements to structure-based virtual screening are usually differences in AUROC or
enrichment between two protocols on one benchmark. We asked what a benchmark of typical size
can actually resolve, and then measured six interventions against that limit with the reading
rule for each committed to version control before the data existed.

Two scoring functions applied to identical poses, agreeing at r = 0.937 and never differing by
more than 2 kcal/mol, produce AUROC values 0.039 apart. That sets a resolution floor of roughly
0.04 AUROC for a ~750-compound benchmark — the range in which most reported docking
improvements fall. Measured against it on SARS-CoV-2 Mpro: changing the scoring function moved
AUROC +0.140, repairing a receptor preparation defect +0.045, scoring the pose ensemble rather
than the top pose +0.019, substituting DiffDock-L pose generation +0.017, supplying the correct
binding site to a co-folding model −0.013, and eight-fold search effort −0.019 with an
unchanged top-ten shortlist. Five of six fall inside the floor.

Against seven free physicochemical descriptors, docking loses on both targets we assembled
(Mpro 0.4530 vs 0.7654; Factor Xa 0.6775 vs 0.7041). To test whether that reflects our panels,
we applied the same descriptors to LIT-PCBA: median AUROC 0.7260 across 15 targets, against
published AutoDock Vina at 0.61 and GNINA at 0.61–0.62 on the same benchmark. On the
AVE-debiased variant descriptors fall to 0.6510, so part of the full-set figure is property
bias that debiasing removes — and part is not.

One method exceeded the floor. Boltz-2, an open-weights co-folding model scoring compounds from
sequence and SMILES alone, does not beat the descriptor baseline on its own (Mpro +0.026,
Factor Xa +0.019, both intervals spanning zero) but adds +0.093 [95% CI +0.062, +0.125] and
+0.075 [+0.051, +0.099] respectively when combined with those descriptors. Its advantage is not
that it is less property-driven — descriptors explain 63% of its score and 63% of Vina's — but
that its residual, the part orthogonal to those descriptors, ranks actives at 0.642 where every
docking score we measured sat between 0.507 and 0.555.

---

## 1. Introduction

A virtual screening method is normally shown to work by reporting enrichment on a retrospective
benchmark. Two properties of that evidence are rarely reported alongside it: what difference the
benchmark could distinguish from noise, and which component of a multi-part protocol produced
the difference. Without the first, an improvement cannot be separated from resampling. Without
the second, a pipeline change is attributed to whichever part its authors were working on.

Analogue and decoy bias in retrospective sets is well documented [1,2], and LIT-PCBA [3] was
constructed to control it. That literature establishes that reported enrichment can measure
property matching. It does not establish how much of a *given* reported improvement is
resolvable at all, which is the question here.

This paper reports three things. A resolution limit, measured rather than assumed. A
decomposition of six interventions against it, each with its reading rule fixed in advance. And
a descriptor baseline on LIT-PCBA, to test whether our own panels are unrepresentative.

## 2. Methods

**Benchmarks.** Two panels assembled from ChEMBL with measured actives and measured inactives:
SARS-CoV-2 Mpro (PDB 7VU6, 751 compounds, 34.3% active) and Factor Xa (854–886 compounds,
45.8% active). LIT-PCBA [3] full and AVE-debiased variants, 15 targets, used as distributed.

**Descriptor baseline.** Molecular weight, clogP, hydrogen-bond donors, hydrogen-bond
acceptors, rotatable bonds, topological polar surface area, formal charge. Logistic regression,
out-of-fold predictions under 5-fold stratified cross-validation. On the AVE variant the
authors' own train/validation split is honoured rather than re-folded. Inactives capped at
25,000 per LIT-PCBA target by seeded random sampling.

**Docking.** AutoDock Vina 1.2.5, exhaustiveness 32, 22 Å box, receptor passing a redocking
control at the screening protocol. gnina 1.3 CNN rescoring applied to the identical Vina output
poses. DiffDock-L for the pose-source arm, with poses minimised before scoring.

**Boltz-2.** Version 2.2.1, protein sequence plus ligand SMILES, `affinity_probability_binary`
as the endpoint, fixed in advance. `diffusion_samples = 1` for cost. Pocket conditioning, where
applied, supplied as 1-based sequence indices derived from residues within 8 Å of the docking
box centre.

**Pre-registration.** Every reading rule — endpoint, threshold, control, abort condition — was
committed to a public git repository before the corresponding data existed, and commit hashes
are cited in the supplementary material. Amendments are appended rather than edited, with the
reason recorded.

**Statistics.** Paired bootstrap over compounds, 4,000 resamples, 95% percentile intervals. A
difference is reported as demonstrated only when it exceeds the measured resolution floor *and*
its interval excludes zero.

## 3. Results

### 3.1 A ~750-compound benchmark cannot resolve below about 0.04 AUROC

gnina's `--score_only` re-evaluation of the identical Vina poses reproduces the cached Vina
score at r = 0.9369, with a mean difference of −0.043 kcal/mol and no compound differing by
more than 2 kcal/mol. The two score vectors nonetheless give AUROC 0.4182 and 0.3793 — **0.039
apart**.

Two scores this close are not meaningfully different descriptions of the same poses, yet the
metric separates them by an amount comparable to most published improvements. We therefore
treat 0.04 AUROC as the resolution floor for a benchmark of this size and report any smaller
difference as not demonstrated, irrespective of its confidence interval.

### 3.2 Five of six interventions fall inside the floor

| intervention | ΔAUROC | 95% CI | resolvable |
|---|---|---|---|
| Scoring function, identical poses (gnina CNN) | **+0.140** | [+0.091, +0.187] | yes |
| Receptor preparation repair | **+0.045** | [+0.024, +0.067] | marginal |
| Pose ensemble rather than top pose | +0.019 | [−0.017, +0.055] | no |
| Pose generator substituted (DiffDock-L) | +0.017 | [−0.027, +0.064] | no |
| Binding site supplied to Boltz-2 | −0.013 | [−0.026, +0.000] | no |
| Eight-fold search effort | −0.019 | top-ten unchanged | no |

The two interventions that clear the floor concern *what is computed*, not *where or how
thoroughly the ligand is placed*. The pose-source result is the sharpest: Vina and DiffDock
poses correlate at r = 0.139 — essentially unrelated coordinates — and produce the same
chance-level ranking.

**The receptor defect is worth reporting for its own sake.** Every receptor we had prepared
carried zero hydrogen-bond donors: the converter typed all 754 Mpro nitrogens as acceptors and
added no polar hydrogens, so the scoring function could not form a protein-donor hydrogen bond
at all. Repairing it moved AUROC +0.045 and changed six of the top ten shortlisted compounds.
A defect invisible in every summary statistic reordered most of the answer.

### 3.3 Seven descriptors beat published docking on LIT-PCBA

| | median AUROC | targets |
|---|---|---|
| descriptors, LIT-PCBA full | **0.7260** | 15 |
| descriptors, full, ≥50 actives | 0.6830 | 10 |
| descriptors, AVE-debiased | **0.6510** | 15 |
| AutoDock Vina, published [4] | 0.61 | 15 |
| GNINA, published [5] | 0.61–0.62 | 15 |

On the full set, like for like, seven descriptors costing microseconds exceed published Vina by
0.116. Our own panels' descriptor baselines (Mpro 0.7654, Factor Xa 0.7041) fall inside
LIT-PCBA's 0.577–0.864 range, so they are not unrepresentative.

AVE debiasing reduces the descriptor median to 0.6510, and this is robust to excluding targets
with few validation actives (0.6510, 0.6497, 0.6510 at ≥0, ≥20, ≥50 actives). **Roughly 0.075
of the full-set figure is therefore property bias that the AVE procedure removes.** The
remainder is not: 0.651 from seven trivial properties on a deliberately debiased set is well
clear of chance.

**A comparison we cannot make.** No published Vina number exists for the AVE-debiased set.
Comparing our AVE descriptor result against a full-set docking result would not be like for
like and we do not do so. Vina would presumably also fall under debiasing; by how much is
unknown to us. The descriptors-beat-docking comparison stands on the full set only.

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
out-of-fold, Boltz-2's remainder ranks actives at **0.6418**; Vina's at 0.5382. Every docking
score we measured — Vina raw and repaired, two gnina heads, across both targets — produced a
residual between 0.507 and 0.555. The ceiling we observed for docking is a property of the
scoring approach, not of the task.

## 4. Limitations

**Two self-assembled targets** for the decomposition, and one of them (Mpro) is a target where
our docking protocol scores below chance. LIT-PCBA broadens the descriptor result to 15 targets
but not the decomposition, which remains a single-benchmark measurement.

**The resolution floor is measured on one benchmark of one size.** It should be re-measured
rather than assumed elsewhere; the method costs one extra scoring pass.

**No Vina number on AVE-debiased LIT-PCBA**, as above.

**The Boltz-2 evaluation is contested territory.** Two 2026 papers disagree about its virtual
screening performance [6,7]. Our contribution is the pre-registration and the residual
decomposition, not the discovery that it works.

**A caveat we withdrew.** Boltz-2's first run was blind while its docking comparator was given
the binding site, and we recorded that as a handicap making +0.093 a lower bound. Supplying the
pocket produced 0.7783 against the blind run's 0.7913, with the two runs correlating at
r = 0.952. The caveat was unnecessary and is withdrawn. It is reported here because a
pre-registered caveat that turns out to be wrong is part of the record.

## 5. Data and code availability

All reading rules, amendments, analysis code and raw outputs are in the project repository,
with the commit hash for each pre-registration cited at the point of use. The analysis script
for the principal Boltz-2 result was committed before the data existed.

## References

[1] Chen et al. (2019) *PLOS ONE* 14:e0220113 — hidden bias in DUD-E.
[2] Sieg, Flachsenberg & Rarey (2019) *J Chem Inf Model* 59:947 — bias control.
[3] Tran-Nguyen, Jacquemard & Rognan (2020) *J Chem Inf Model* 60:4263 — LIT-PCBA.
[4] arXiv:2605.01681 — benchmarking single-pose docking on LIT-PCBA.
[5] Sunseri & Koes (2021) — GNINA scoring on DUD-E and LIT-PCBA.
[6] arXiv:2508.17555 — Boltzina.
[7] arXiv:2603.05532 — reliability of AI methods in drug discovery.
[8] arXiv:2507.21404 — data leakage and redundancy in LIT-PCBA.
