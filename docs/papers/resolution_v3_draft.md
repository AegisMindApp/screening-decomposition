# What a virtual screening benchmark can resolve, and which interventions survive it

**John Goodman** · OceanSparx Pty Ltd, Sydney, Australia · john.goodman@oceansparx.com

**DRAFT v3 — 21 September 2026.** Written fresh. Every number below is traceable to a named
artifact in the project repository; numbers carried from earlier analysis without
re-verification in this cycle are marked **[carried]** and must be re-checked before submission.

---

## Abstract

Improvements to structure-based virtual screening are normally reported as a difference in AUROC
between two protocols on one benchmark. Whether a benchmark of that size can resolve the
difference is rarely asked. We measure the resolution directly, by replicating an *unchanged*
protocol across independent search seeds, and then judge five interventions against it under a
formal equivalence framework.

**The floor, measured.** Ten replicates of an identical AutoDock Vina protocol on a
755-compound SARS-CoV-2 Mpro panel, varying only `--seed`, give an AUROC standard deviation of
**0.00497** (one-sided 95% upper bound **0.00817**). Six replicates at eight-fold higher search
effort give **0.00641** (bound **0.01338**) on a 343-compound panel. Matched on identical
compounds and seeds the two are indistinguishable (Pitman–Morgan p = 0.62), so search effort
does not measurably change reproducibility.

**Floors are panel-size specific.** AUROC variance scales as 1/n. The stage-1 figure predicts
0.00742 on a 339-compound panel and we measure 0.00812 — a ratio of 1.09 from a prediction using
data the check does not otherwise touch. A floor measured on 755 compounds is therefore ~1.5×
too small to charge against a claim computed on 339, and must be rescaled.

**Three interventions survive, one is negligible, one is unresolved.** Against the floor
rescaled to the panel in use: changing the scoring function on identical poses **+0.140**
[+0.089, +0.191] survives; scoring the pose ensemble rather than the top pose **+0.055**
[+0.035, +0.075] survives; repairing a receptor preparation defect **+0.045** [+0.014, +0.076]
survives. Eight-fold search effort, re-measured paired within seed across six seeds, is
**+0.004** [−0.006, +0.014] and is **statistically equivalent to zero** within the resolution
(TOST p = 0.003 / 0.033). Supplying a co-folding model with the correct binding site
(**−0.013** [−0.039, +0.013]) is **inconclusive** — its interval is wider than the margin, so
the data do not decide.

**Chemical series matter for absolutes, not for these differences.** The panel holds 388 Butina
series over 749 compounds with a Kish effective n of 136.5 — about one fifth of the compound
count. Resampling series rather than compounds widens an *absolute* AUROC interval by **1.63×**,
but widens a *paired difference* by only **1.04×**, because both arms carry the same series
structure and it cancels. Every intervention here is a paired difference.

**The decision rule detects effects of known size.** Five effects planted at +0.000, +0.004,
+0.050, +0.140 and −0.050 are recovered with the verdict fixed in advance, including under
realistic reseed noise.

**Method floors are not interchangeable.** gnina CNN rescoring is bitwise deterministic
(745/745 identical) and contributes no run-to-run noise. Boltz-2 is not: zero of 24
duplicated-ligand pairs agreed to five decimals, giving a residual-AUROC floor of **0.00592** on
an 885-compound panel. Borrowing one method's floor for another is an error of a factor of
four in this data set.

---

## 1. Introduction

A published virtual screening improvement is typically a difference of 0.02–0.10 AUROC between
two protocols evaluated once each on one benchmark. The implicit claim is that the difference is
a property of the protocols. It is only that if the benchmark can resolve a difference of that
size — and the resolution is a measurable quantity that is almost never measured.

The right measurement is narrow. Take one protocol, change nothing about it except the random
seed, run it repeatedly, and record how much the evaluation metric moves. Anything smaller than
that spread is not a difference between protocols; it is the same protocol disagreeing with
itself.

This is easy to get wrong in ways that look reasonable. A floor taken from the disagreement
between two *implementations* of a scoring function measures implementation variance, not
reproducibility. A floor taken from comparing two *search-effort conditions* measures an
intervention effect, and if search effort is also one of the interventions being judged, the
intervention is being evaluated against itself. Neither quantity answers the question, and both
are easily mistaken for it.

We therefore measure the floor by seed replication of an unchanged protocol, establish how it
scales with panel size, and then judge interventions against it with a rule that distinguishes
three outcomes rather than two: an effect larger than the resolution, an effect measured to be
negligible, and an effect the data cannot resolve. The third category is usually reported as if
it were the second.

## 2. Methods

### 2.1 Panel and protocol

SARS-CoV-2 main protease (Mpro), receptor 7VU6, ChEMBL target CHEMBL4523582. 755 compounds
(257 active / 498 inactive) by pChEMBL threshold with a fixed cutoff year. AutoDock Vina, box
centred on the co-crystallised ligand, exhaustiveness 4 (stage 1) and 32 (stage 2). Everything
fixed but `--seed`.

Compounds that Vina scores as a whole number are a parsing hazard: a regular expression
requiring a decimal point silently records them as docking failures, at roughly 1 in 250. The
runner asserts receptor and binary checksums against a manifest and records per-compound
failures explicitly.

### 2.2 The resolution floor

Stage 1: ten seeds at exhaustiveness 4 on all 755 compounds. Stage 2: six seeds at
exhaustiveness 32 on a fixed 350-compound subset (343 scored in every seed), sharded four ways
per seed. Each seed's AUROC is computed on the intersection of compounds scored in *every* seed,
because a partial seed truncates the panel for all of them.

The floor is reported as the χ² one-sided 95% upper bound on the standard deviation, not the
point estimate. A point sd from ten replicates understates the spread about half the time.

### 2.3 Panel-size scaling

AUROC variance scales as 1/n. This is verified rather than assumed: the stage-1 sd measured on
755 compounds predicts a value on a 339-compound panel, and the matched re-measurement is
compared against that prediction (§3.2). Floors are rescaled by √(n_measured/n_claim) before
being charged against any claim computed on a different panel.

### 2.4 Equivalence framework

With the floor as the equivalence margin *d*, each effect is classified:

| verdict | criterion | meaning |
|---|---|---|
| SURVIVES | 95% CI lower bound > +*d* | exceeds the resolution |
| NEGLIGIBLE | 90% CI ⊂ (−*d*, +*d*) | TOST: equivalent to zero within resolution |
| INCONCLUSIVE | otherwise | the data do not decide |

The 90% interval is the standard two-one-sided-tests construction at α = 0.05. Using the 95%
interval there is a common error that makes equivalence harder to declare than it should be.

The distinction between NEGLIGIBLE and INCONCLUSIVE is the point of the framework. Reporting
both as "no significant effect" claims knowledge in the second case that the data do not supply.

### 2.5 Chemical-series structure

Butina clustering on Morgan fingerprints (radius 2, 2048 bits) at Tanimoto 0.65. Intervals are
computed two ways — resampling compounds, and resampling whole series — and compared, for both
an absolute AUROC and a paired difference.

### 2.6 Positive control

Effects of known size are planted by shifting active scores to hit a target AUROC delta, then
passed through the identical analysis. Expected verdicts are fixed before running, and include a
sign-flipped case so the rule cannot be one-sided. A second tier carries each planted effect on
a *different* seed's run, so the comparison includes genuine reseed noise rather than only
planting precision.

## 3. Results

### 3.1 The floor

| protocol | seeds | panel | AUROC sd | 95% upper bound |
|---|---|---|---|---|
| Vina, exhaustiveness 4 | 10 | 755 | 0.00497 | **0.00817** |
| Vina, exhaustiveness 32 | 6 | 343 | 0.00641 | **0.01338** |

Per-compound behaviour at exhaustiveness 4: score sd **0.096 eV**, median rank span **70 places
of 755**, maximum 549. At exhaustiveness 32 the per-compound reproducibility is far better —
score sd **0.035 eV**, median rank span **10** — while the AUROC sd is similar, because AUROC sd
is set by panel size rather than by per-compound precision.

Matched on the same 339 compounds and the same six seeds, exhaustiveness 4 gives sd 0.00812 and
exhaustiveness 32 gives 0.00626 (ratio 0.77). A Pitman–Morgan paired test — paired rather than
F, because the arms share seeds and their AUROCs are correlated — gives **p = 0.62**. The two
floors are not distinguishable at this replication, and the common assumption that less search
is more stochastic is not supported.

### 3.2 Floors are panel-size specific

The stage-1 exhaustiveness-4 sd of 0.00497 on 755 compounds predicts **0.00742** on a
339-compound panel under 1/n scaling. The matched re-measurement gives **0.00812**, a ratio of
**1.09**, from a prediction using data the check does not otherwise use.

A floor measured on 755 compounds is therefore about 1.5× too small to charge against a claim
computed on 339, and using it there would admit effects that are inside the noise. This is a
general point about reporting: a resolution limit is a property of *(method, protocol, panel
size)*, not of the method.

### 3.3 Five interventions

Against the exhaustiveness-32 floor rescaled to the panel in use:

| intervention | Δ AUROC | 95% CI | verdict |
|---|---|---|---|
| scoring function (gnina CNN) on identical poses | **+0.140** | [+0.089, +0.191] | SURVIVES |
| pose ensemble rather than top pose | **+0.055** | [+0.035, +0.075] | SURVIVES |
| receptor preparation repair | **+0.045** | [+0.014, +0.076] | SURVIVES |
| eight-fold search effort | **+0.004** | [−0.006, +0.014] | **NEGLIGIBLE** |
| correct binding site to a co-folding model | −0.013 | [−0.039, +0.013] | **INCONCLUSIVE** |

The search-effort result is re-measured **paired within seed** across six seeds that were run at
both exhaustiveness settings, which is the comparison the design supports: pairing cancels what
the two arms share. TOST against the margin gives p = 0.0033 and p = 0.0328, so the effect is
equivalent to zero within the resolution. The verdict is close to its boundary — the 90% upper
bound is +0.0121 against a margin of 0.01346 — and should be read as such rather than as a
comfortable null.

Pairing was worth doing and that too is measured, not assumed: the paired sd of the difference
is 0.00968 against an unpaired expectation of √2 × sd = 0.01148.

### 3.4 Chemical series affect absolutes, not paired differences

388 series over 749 compounds; 293 singletons; largest 24; **Kish effective n 136.5**.

| quantity | compound bootstrap | cluster bootstrap | widening |
|---|---|---|---|
| paired difference (exh32 − exh4) | 0.04075 | 0.04230 | **1.04×** |
| absolute AUROC (exh=4 arm) | 0.1338 | 0.2180 | **1.63×** |

A paired difference sees the same series in both arms and the clustering cancels. An absolute
AUROC has no such protection: its compound-resampled interval is 63% too narrow. All five
interventions are paired differences, so their intervals stand; any absolute AUROC reported from
a clustered panel requires the cluster interval.

This is offered in place of a scaffold-disjoint split, which would answer a different question.
The subject here is the size of an intervention on a fixed panel, not generalisation to unseen
chemotypes — and for that quantity pairing, not splitting, is what controls series confounding.
No claim is made about performance on a new chemical series.

### 3.5 The decision rule detects planted effects

| planted | recovered | verdict | expected |
|---|---|---|---|
| +0.000 | +0.00001 | NEGLIGIBLE | NEGLIGIBLE |
| +0.004 | +0.00402 | NEGLIGIBLE | NEGLIGIBLE |
| +0.050 | +0.05001 | SURVIVES | SURVIVES |
| +0.140 | +0.14003 | SURVIVES | SURVIVES |
| −0.050 | −0.05000 | SURVIVES | SURVIVES |

Carried on a *different* seed's run, so that the comparison includes real reseed noise, a
planted +0.004 becomes INCONCLUSIVE while +0.050 and +0.140 remain SURVIVES. That is the
expected behaviour and it quantifies what the pairing buys: the same effect size is resolvable
paired and unresolvable unpaired.

The control failed on first run and the failure was informative: Vina affinities are
lower-is-better, the planting code assumed the opposite, and it recovered −0.446 for a planted
+0.004 and inverted the sign of a planted −0.050. A control that cannot fail is not a control.

### 3.6 Method floors are not interchangeable

| method | floor | basis |
|---|---|---|
| Vina, exh 4, 755 compounds | 0.00817 | 10 seed replicates |
| Vina, exh 32, 343 compounds | 0.01338 | 6 seed replicates |
| gnina CNN rescoring | **0.0** | 745/745 bitwise identical, measured |
| Boltz-2, 885-compound panel | **0.00592** | 24 duplicated-ligand pairs, propagated |

gnina's zero is measured, not assumed: repeat runs agree bitwise, and so do runs at different
seeds, because `--score_only` performs no stochastic search. Boltz-2 is stochastic — zero of 24
duplicated pairs agreed to five decimals, maximum difference 0.0905 in binding probability, with
noise heteroscedastic in the predicted value (r = +0.80 with distance from the 0/1 bounds).

Boltz-2's floor was previously taken to be 0.0201, which is the Vina scoring floor. It is about
four times too large. A borrowed threshold is not safe merely because it is conservative.

### 3.7 A co-folding model across three admissible targets **[carried, partially re-verified]**

Targets are admissible only if seven physicochemical descriptors do not already separate the
panel (5-fold CV AUROC ≤ 0.80); on an inadmissible panel nothing can demonstrate that a method
exceeds properties. Mpro (0.7654), Factor Xa (0.7041) and ALDH1 (0.5698) are admissible. PD-L1
(0.874) is not, and is excluded.

With residual AUROC measured against a per-target null — the residual of the method's own
out-of-fold descriptor prediction, which is **not** 0.5 — Boltz-2 gives margins of +0.0571
(Mpro), +0.0600 (ALDH1) and **+0.0042** (Factor Xa) against its measured floor of 0.00592. Two
targets clear comfortably; Factor Xa does not, and by a margin finer than the floor's own
precision.

Cold transfer is usable (absolute AUROC ≥ 0.60) in three of six directions. The target where the
model is weakest in-target, ALDH1 at 0.5629, is its best donor — transferring at 0.727 into Mpro
and 0.665 into Factor Xa — and its hardest recipient, with nothing transferring in above 0.58.
That is consistent with ALDH1 having the least property signal for a model to lean on.

## 4. Discussion

The measured resolution of this benchmark is roughly 0.008–0.013 AUROC. Published virtual
screening improvements are routinely reported in the 0.02–0.10 range, so the majority are
nominally resolvable — but the two smallest interventions here, at +0.004 and −0.013, are not,
and both would have been reported as findings under the usual convention of quoting a point
estimate and a bootstrap interval.

Three conclusions generalise beyond this panel.

**A resolution limit must be measured on the axis it polices.** Implementation variance and
search-effort effects are both real quantities and neither is run-to-run reproducibility. Using
either as a floor mismeasures the thing the decision rule depends on, and if the borrowed
quantity is itself one of the interventions, the intervention is judged against itself.

**A floor is specific to a panel size.** Since AUROC variance scales as 1/n, a floor is not a
property of a method but of a method on a panel. Transporting one between panels of different
size is a silent error of up to 1.5× in this work.

**Equivalence and non-detection are different results.** An interval inside the margin is a
measurement that the intervention does nothing worth having. An interval wider than the margin
is an absence of information. Collapsing them into "no significant difference" reports the
second as though it were the first.

## 5. Limitations

- One target for the intervention work. The floor is Mpro's, and there is no reason to expect a
  single number to transport to other targets or panels; §3.2 shows it does not even transport
  between panel sizes without rescaling.
- The stage-2 floor rests on six seeds and the exhaustiveness verdict sits near its margin
  boundary. More seeds would tighten both.
- Boltz-2's floor is derived by propagating per-compound noise measured on 24 duplicated pairs
  through the full panel, not by full-panel replication. The Factor Xa verdict turns on the third
  decimal place; 100 pairs would settle it and is inexpensive.
- Propagation assumes independent per-compound noise. Correlated run-level drift is largely
  irrelevant to a rank-based metric, but structured differences would not be captured.
- No claim is made about generalisation to unseen chemical series; §3.4 addresses series
  confounding of the reported differences, not transfer to new chemotypes.
- Descriptor-baseline comparisons and the combination results **[carried]** have not been
  re-verified in this cycle and must be before submission.

## 6. Data and code availability

All analysis code, pre-registrations and per-seed results are in the project repository:
`analysis/seed_floor/` (floor, equivalence, positive control, series structure),
`analysis/boltz2/` (co-folding runs and floor), `analysis/method_bench/` (admissibility and
transferability harness). Pre-registrations are committed before their corresponding data in
every case; the reader's ability to verify that depends on the repository history being
publicly inspectable, which must be arranged before submission.
