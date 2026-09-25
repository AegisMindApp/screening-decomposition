# Clearing a benchmark's resolution limit is necessary and not sufficient:
# measured floors for structure-based virtual screening, and five interventions on a second target

**John Goodman** · OceanSparx Pty Ltd, Sydney, Australia · john.goodman@oceansparx.com

**Version 3 — 25 September 2026.** Supersedes version 2 of this record (9 September 2026); what
changed, and why, is set out under *Relation to the earlier version*. Every number is recomputed
from data in the public repository, and every verdict follows a rule fixed in writing before the
measurement it decides. Where a claim from the earlier version did not survive, it is named as a
correction beside the result that replaces it — §3.3 for the resolution limits, §3.6 for the
scoring-function result — and both are collected under *Relation to the earlier version*.

---

## Abstract

A virtual screening improvement is normally reported as a difference in AUROC between two
protocols, each run once, on one benchmark. Two questions are almost never asked: whether the
benchmark can resolve a difference of that size, and whether the difference survives a second
target. We ask both, and the answers are not the same.

**The resolution is measurable and predictable.** Fifteen replicates of an identical AutoDock
Vina protocol on a 755-compound SARS-CoV-2 Mpro panel, varying only `--seed`, give an AUROC
standard deviation of 0.00441 and a χ² one-sided 95% bound of **0.00644**. Twelve replicates on
an 886-compound Factor Xa panel give **0.00682**. AUROC variance scales as 1/*n*, verified within
a target (755 → 339, ratio 1.14) and, on a prediction fixed in advance, across targets
(0.00459 predicted, 0.00427 measured, ratio 0.93). A floor can therefore be quoted for a panel
that has not been re-measured.

**It is specific to the method, the protocol and the panel, and borrowing is expensive.** Three
substitutions we made and corrected cost factors of 4.0 (one method's floor used for another),
2.2 (the sd of a residual charged against a margin), and 1.5 (a floor moved between panel sizes).
The floor is only weakly specific to the AUROC *level*: across simulated panels from 0.28 to 0.90
it varies by 12%, so a floor measured where a method performs poorly is conservative elsewhere.

**The metric matters more than any of it.** On the same replicates, EF@1% has a coefficient of
variation **8.6×** AUROC's — a floor of 0.38 against a mean of 2.28, and a third of the value on
a smaller panel. Enrichment at an early cutoff is what screening papers usually report, and it is
the quantity a benchmark of this size resolves worst. BEDROC(α=20) sits between the two.

**Clearing the floor does not make a result real, and that is the paper's central claim.** Of
five interventions measured on Mpro, three clear the floor by wide margins. Re-measured on a
second target, **none of the four testable ones replicates**. The largest — replacing the scoring
function on identical poses, **+0.140 on Mpro — twenty-two times that panel's floor, and its
whole 95% interval above it** — comes back at **−0.058** on Factor Xa: a reversal of sign, with
the interval clear of that target's own floor. Of the two mid-sized effects, one falls to about a
quarter of its Mpro value and one changes sign; neither is resolvable. The smallest effect
reproduces almost exactly (+0.0082 against +0.0096) and is unresolvable on both.

A resolution limit is a necessary filter and a poor certificate. We recommend measuring it,
reporting it with its panel size, and treating a result that clears it as *eligible* for belief
rather than established.

## 1. Introduction

A published virtual screening improvement is commonly a difference of order 0.02–0.10 AUROC
between two protocols evaluated once each on one benchmark. That range is an impression from
reading the field, not a systematic survey, and no conclusion here depends on it: every claim
below concerns quantities measured on our own panels.

The implicit claim in such a report is that the difference is a property of the protocols. It is
only that if two conditions hold. The benchmark must be able to resolve a difference of that
size — otherwise the number is the protocol disagreeing with itself. And the difference must be a
property of the method rather than of the one protein it was measured on.

The first condition has a cheap and exact test that is almost never run: take one protocol, change
nothing but the random seed, run it repeatedly, and record how far the metric moves. The second
has an obvious test that is run even less often on methodological claims: do it again somewhere
else.

This is easy to get wrong in ways that look careful. A floor taken from the disagreement between
two *implementations* of a scoring function measures implementation variance. A floor taken from
comparing two *search-effort settings* measures an intervention effect — and if search effort is
also one of the interventions being judged, the intervention is being scored against itself.
Neither is run-to-run reproducibility, and both are easily mistaken for it.

We measure the floor by seed replication of an unchanged protocol; establish how it scales with
panel size, method, protocol, metric and performance level; judge five interventions against it
with a rule distinguishing three outcomes rather than two; and then take four of those
interventions to a second target.

The two halves give different answers, and the paper exists because of the gap between them.

## 2. Methods

### 2.1 Panels and protocol

**Mpro:** SARS-CoV-2 main protease, receptor 7VU6, ChEMBL target CHEMBL4523582, 755 compounds
(257 active) by pChEMBL threshold with a fixed cutoff year. **Factor Xa:** ChEMBL244, 886
compounds (405 active), receptor protonated with `obabel -xr -p 7.4`. Both panels are built from
ChEMBL activity data by pChEMBL threshold — **not** from a curated unbiased benchmark such as
LIT-PCBA [6], and we do not claim the debiasing that set applies. Threshold-defined actives carry
a known risk that measured enrichment reflects property matching rather than molecular
recognition [7]. §3.5 addresses that directly: seven physicochemical descriptors reach 0.714 on
the Mpro panel against Vina's 0.416, so on this benchmark the properties win outright, and every
intervention result is a statement about paired differences on a panel where that is true. AutoDock Vina [1,2], box
centred on the co-crystallised ligand, exhaustiveness 4 and 32. Everything fixed but `--seed`.

Compounds Vina scores as a whole number are a parsing hazard: a regular expression requiring a
decimal point silently records them as failures, at roughly 1 in 250. The runner asserts receptor
and binary checksums against a manifest on every host and records per-compound failures
explicitly.

**Vina is below chance on Mpro** — mean AUROC 0.4163 — and above it on Factor Xa, at 0.6703. We
report this here rather than leaving it to be discovered, because it governs how the intervention
results should be read (§3.6) and because it is the mechanism behind the paper's main finding.

### 2.2 The floor

Seed replication of an unchanged protocol: fifteen seeds at exhaustiveness 4 on Mpro, ten at
exhaustiveness 32 on a fixed 350-compound subset (343 scored in every seed), twelve at
exhaustiveness 4 on Factor Xa. Each seed's metric is computed on the intersection of compounds
scored in *every* seed, because a partial seed truncates the panel for all of them.

The floor is reported as the **χ² one-sided 95% upper bound** on the standard deviation, not the
point estimate. A point sd from a handful of replicates understates the spread about half the
time, and the penalty is large: 2.09× at *n* = 6, 1.65× at 10, 1.46× at 15.

### 2.3 Scaling, and what a floor belongs to

A floor is a property of **(method, protocol, panel size, metric)**. Panel size is the one
dimension with a law: AUROC variance scales as 1/*n*, so floors are rescaled by
√(*n*_measured / *n*_claim) before being charged against a claim computed elsewhere. The other
dimensions have no law and must be measured; §3.3 quantifies what ignoring them costs.

### 2.4 Equivalence framework

With the floor as the equivalence margin *d*:

| verdict | criterion | meaning |
|---|---|---|
| SURVIVES | 95% CI lower bound > +*d* | exceeds the resolution |
| NEGLIGIBLE | 90% CI ⊂ (−*d*, +*d*) | TOST: equivalent to zero within resolution |
| INCONCLUSIVE | otherwise | the data do not decide |

The 90% interval is the standard two-one-sided-tests construction at α = 0.05 [5]. Reporting
NEGLIGIBLE and INCONCLUSIVE together as "no significant difference" claims knowledge in the
second case that the data do not supply.

### 2.5 Replication on a second target

Four of the five interventions were re-measured on Factor Xa under rules fixed before any of that
data existed, including which arms could be replicated and which could not. Each is charged
**Factor Xa's own floor**, never Mpro's. Verdicts: REPLICATES (same sign, interval clears that
target's floor), DOES NOT REPLICATE (opposite sign, clearing the floor), INCONCLUSIVE (straddles).
The pre-registration states that an INCONCLUSIVE replication is reported as prominently as a
success and **does not count as support**.

### 2.6 Controls

Planted effects of known size with verdicts fixed in advance, including a sign-flipped case; a
second tier carries each planted effect on a *different* seed's run so the comparison includes
real reseed noise. Chemical-series structure is assessed by Butina clustering [4] on Morgan
fingerprints (radius 2, 2048 bits, Tanimoto 0.65), with intervals computed by resampling
compounds and by resampling whole series. Each replication arm carries its own control: a
degenerate-identity check, a permuted-feature check, a dropout check, or a scorer-identity check,
as the arm requires.

## 3. Results

### 3.1 The floor

| protocol | target | seeds | panel | mean AUROC | sd | **95% bound** |
|---|---|---|---|---|---|---|
| Vina, exh 4 | Mpro | 15 | 755 | 0.4163 | 0.00441 | **0.00644** |
| Vina, exh 32 | Mpro | 10 | 343 | 0.4531 | 0.00619 | **0.01019** |
| Vina, exh 4 | Factor Xa | 12 | 886 | 0.6703 | 0.00440 | **0.00682** |

Per-compound behaviour at exhaustiveness 4 on Mpro: score sd **0.101 eV**, median rank span **84
places of 755**. At exhaustiveness 32 the per-compound reproducibility is far better — score sd
0.041 eV, median rank span 15 — while the AUROC sd is *similar*, because AUROC sd is set by panel
size rather than by per-compound precision. Matched on 339 compounds and shared seeds, the two
floors are not distinguishable (Pitman–Morgan p = 0.38): the common assumption that less search
is more stochastic has no support here.

**Replication buys the bound, not the estimate.** Adding seeds barely moves the sd and steadily
tightens the bound, because the χ² penalty shrinks:

| Mpro exh 4, *n* | 8 | 10 | 12 | 15 |
|---|---|---|---|---|
| sd | 0.00526 | 0.00497 | 0.00485 | **0.00441** |
| bound | 0.00946 | 0.00817 | 0.00751 | **0.00644** |

Since the bound is what gets charged against a claim, **a floor from few replicates is loose in
the direction that admits effects**. This is not a minor point of practice: one verdict in this
paper changed when its floor was re-measured with more seeds (§3.5).

### 3.2 Scaling: 1/n holds within and across targets

| check | *n* | predicted | measured | ratio |
|---|---|---|---|---|
| Mpro, within target (755 → 339) | 339 | 0.00742 | 0.00844 | 1.14 |
| **Factor Xa, across targets (755 → 886)** | 886 | **0.00459** | **0.00427** | **0.93** |

The cross-target prediction and its ±25% acceptance band were **fixed in writing before any
Factor Xa docking ran**, which is what makes it a test rather than a fit. It passed.

Two caveats, both stated rather than discovered. The sd estimator is imprecise at ten replicates —
the two-sided 95% interval on the Factor Xa value is [0.00294, 0.00780], **wider than the
acceptance band** — so the test is asymmetric: clearing a band narrower than the measurement's own
precision is strong evidence, failing it would have been weak. And the verdict turned on the tenth
seed: at *n* = 9 the sd sat 4×10⁻⁶ outside the lower edge. Two further seeds settle it — at
*n* = 12 the sd is 0.00440, ratio 0.96, sitting 21% of the prediction clear of that edge. The
pre-registered reading at *n* = 10 stands as published and is not re-read; the later seeds are
reported as robustness, because retro-fitting a prediction to a larger sample destroys what made
it worth running.

### 3.3 What a floor belongs to: three substitutions, and what each cost

Every one of these is an error we made and corrected, which is why the sizes are known.

| substitution | correct value | value used | error |
|---|---|---|---|
| one **method's** floor for another (Vina's for Boltz-2) | 0.00472 | 0.0201 | **4.0×** |
| the sd of a **residual** charged against a **margin** | 0.01018 | 0.00472 | **2.2×**, permissive |
| a floor moved between **panel sizes** without rescaling | 0.00844 | 0.00497 | **1.5×** |

The second is the subtlest and the most dangerous, because both quantities are legitimate and the
names are nearly the same. Our verdicts are charged against a *margin* — residual AUROC minus its
descriptor null — but the published floor was the sd of the *residual alone*. Propagating the same
measured per-compound noise to the margin gives a figure 2.2× larger, because the residual and its
null do not move together under re-running. Correcting it cost one intervention its verdict.

**One substitution turned out to be safe, and we could not have known without measuring.** Boltz-2's
per-compound noise is target-specific — 0.01059 on Factor Xa against 0.02341 on Mpro, a factor of
2.2 — and we predicted that carrying its floor between targets would be correspondingly wrong.
It is not: the measured Mpro margin floor is **0.01098** against **0.01107** for Factor Xa's
rescaled by panel size, **0.8% apart**. AUROC is rank-based, so what governs it is noise relative
to the spread of the scores it must reorder; Mpro's noise is 2.2× larger and its score
distribution 1.77× wider, leaving 1.25× in the units the metric responds to. A ratio measured on
one quantity is not a ratio on another, and the propagation connecting them is cheap to run.

### 3.4 The metric dominates: EF@1% is the least resolvable thing the field reports

Same seeds, same panels, four metrics:

| exh | metric | mean | 95% floor | CV = sd/mean |
|---|---|---|---|---|
| 4 | AUROC | 0.4174 | 0.00817 | **0.0119** |
| 4 | **EF@1%** | 2.2767 | **0.3821** | **0.1020** |
| 4 | EF@5% | 1.3142 | 0.1696 | 0.0784 |
| 4 | BEDROC(α=20) | 0.4489 | 0.0275 | 0.0372 |
| 32 | **EF@1%** | 2.5751 | **0.8931** | **0.2108** |

The pre-registered prediction was CV(EF@1%) > 3 × CV(AUROC). It is **8.6×** at exhaustiveness 4
and 15.4× at 32. In absolute terms, reseeding an unchanged protocol moves EF@1% by up to **0.38 of
an enrichment factor of 2.28** — a sixth of the quantity — and by a third of it on the
343-compound panel. The mechanism is not subtle: EF@1% on 755 compounds is decided by the **top
eight compounds**, and which eight they are is not stable under reseeding.

BEDROC(α=20), which weights early recognition smoothly rather than at a hard cutoff, sits at CV
0.037 — three times AUROC's and a third of EF@1%'s. **Where early recognition must be reported,
BEDROC is the more resolvable choice**, and that is a recommendation this measurement supports
rather than a preference.

### 3.5 Five interventions on Mpro: three clear the floor

Each row is charged the floor for **its own** method, protocol and panel.

| intervention | Δ AUROC | 95% CI | *n* | floor | *d* | verdict |
|---|---|---|---|---|---|---|
| scoring function (gnina CNN [3]) on identical poses | **+0.140** | [+0.089, +0.191] | 745 | Vina exh 4 | 0.00648 | **SURVIVES** |
| pose ensemble rather than top pose | **+0.055** | [+0.035, +0.075] | 751 | Vina exh 4 | 0.00646 | **SURVIVES** |
| receptor preparation repair | **+0.045** | [+0.014, +0.076] | 753 | Vina exh 4 | 0.00645 | **SURVIVES** |
| eight-fold search effort | +0.008 | [+0.002, +0.015] | 339 | Vina exh 32 | 0.01025 | INCONCLUSIVE |
| correct binding site to a co-folding model [8] | −0.013 | [−0.039, +0.013] | 750 | Boltz-2 margin, measured | 0.01098 | INCONCLUSIVE |

Two of these verdicts moved during this work, and both movements are instructive rather than
embarrassing. **Search effort** was NEGLIGIBLE at six replicates and is INCONCLUSIVE at ten: the
point estimate doubled and the margin tightened 24% as the χ² penalty fell, and an equivalence
verdict depends on the margin as much as on the effect. **The co-folding row** was charged a Vina
floor until we noticed that §3.3's own rule forbids it; charging Boltz-2's own measured floor
leaves the verdict unchanged but the number 1.6× larger.

Note what the NEGLIGIBLE column now contains: **nothing**. The category the framework exists to
distinguish is exercised only by the planted controls of §3.7. That is a weaker demonstration
than reporting a real measured-negligible effect, and it is the honest state of the evidence.

The intervention intervals charge run-to-run noise **twice** — once inside the interval, once at
the threshold, both built from the same standard deviation. The consequence runs one way only: it
makes SURVIVES harder and can never inflate a positive.

### 3.6 The same interventions on a second target: none replicates

| intervention | Mpro | **Factor Xa** | *n* | floor | verdict |
|---|---|---|---|---|---|
| scoring function (gnina CNN) | +0.1397 [+0.0913, +0.1867] | **−0.0583 [−0.0930, −0.0235]** | 849 | 0.00697 | **DOES NOT REPLICATE** |
| pose ensemble | +0.0552 [+0.0434, +0.0680] | −0.0135 [−0.0315, +0.0043] | 886 | 0.00682 | INCONCLUSIVE |
| receptor preparation repair | +0.0451 [+0.0236, +0.0671] | +0.0116 [+0.0001, +0.0234] | 772 | 0.00731 | INCONCLUSIVE |
| eight-fold search effort | +0.0082 [+0.0016, +0.0148] | **+0.0096 [+0.0049, +0.0142]** | 861 | 0.00517 | INCONCLUSIVE |

Each arm passed its own control before any AUROC was read: gnina's scorer-identity control
(its own Vina-style affinity against our cached scores) at **r = 0.9934**; the pose-ensemble
degenerate-identity control at **exactly 0.00e+00**, with a permuted-feature control costing
0.24 AUROC; receptor prep's dropout control at **−5.1pp** against a 10pp void threshold.

**The headline result reverses.** gnina rescoring gained +0.140 on Mpro — **twenty-two times that
panel's floor of 0.00648**, with even its 95% lower bound fourteen times the floor, so it is as
clearly resolved as anything in this paper — and loses 0.058 on Factor Xa, with the interval clear
of that target's floor. Both measurements are correct. The mechanism is visible in
the absolute numbers: on Mpro, Vina ranks **worse than random** (0.399) and gnina lifts it across
chance to 0.539; on Factor Xa, where Vina already ranks usefully (0.687), the same rescoring
degrades it to 0.629.

**We tested the obvious explanation and could not establish it.** The natural reading of those
numbers is that the Mpro gain was recovery from a broken baseline rather than a property of the
scoring function — a method that replaces a worse-than-random ordering with a near-random one
gains a great deal of AUROC while producing nothing a screen could use. That is testable without
new data, by asking whether gnina's *absolute* AUROC tracks Vina's across sub-panels: a
constant-quality ranker gives slope 0, a method that tracks target difficulty gives slope 1.
(Regressing the *delta* on the baseline, the tempting version, is worthless — its slope is −1 by
construction for any uncorrelated score.)

Across 19 scaffold-binned sub-panels spanning Vina AUROC 0.241–0.854, the slope is
**+0.577 [+0.173, +0.956]** — straddling the 0.5 midpoint fixed in advance. The positive control
(Vina at a second seed, true slope 1) returns **+0.958 [+0.918, +1.006]**, so the design can
detect tracking and the null is not an artefact of regression dilution; a permuted-gnina control
ranks at chance with a near-zero slope, as it should.

The verdict is INCONCLUSIVE and the explanation is therefore **withdrawn rather than softened**,
under a rule fixed before the analysis ran. **The reversal stands as a measured fact with no
established mechanism.** A reader may weigh the absolute numbers above for themselves; we decline
to supply a story the data do not carry.

**The pattern is not uniform and should not be reported as one.** The largest effect reverses; the
two mid-sized effects shrink to roughly a quarter and cannot be resolved, one pointing the wrong
way; the **smallest effect reproduces almost exactly** — +0.0082 against +0.0096, six of six seeds
positive — and is unresolvable on both targets. That last row is the paper's thesis in its
cleanest form: an effect agreeing to three decimal places across two independent targets, still
indistinguishable on either from re-running the same protocol with a different seed.

Factor Xa gave search effort its best chance and it still did not clear. That target's exh-32
floor is **stricter** than Mpro's (0.00517 against 0.01025) because the run is more reproducible
there, and the interval's lower bound missed it by 0.0003.

### 3.7 Controls, and one that failed usefully

Planted effects, verdicts fixed in advance, re-run against the current margin:

| planted | recovered | verdict | expected |
|---|---|---|---|
| +0.000 | +0.00001 | NEGLIGIBLE | NEGLIGIBLE |
| +0.004 | +0.00399 | NEGLIGIBLE | NEGLIGIBLE |
| +0.050 | +0.05000 | SURVIVES | SURVIVES |
| +0.140 | +0.14002 | SURVIVES | SURVIVES |
| −0.050 | −0.05002 | SURVIVES | SURVIVES |

Tier-1 intervals are ±0.00005, so tier 1 tests the planting arithmetic rather than the rule under
realistic noise. **Tier 2 is the control that can fail**: carrying each planted effect on a
*different* seed's run, a planted +0.004 returns INCONCLUSIVE while +0.050 and +0.140 remain
SURVIVES. That tier predicted the search-effort verdict change (§3.5) before the data did.

The control failed on first run and the failure was real: Vina affinities are lower-is-better, the
planting code assumed the opposite, and it recovered −0.446 for a planted +0.004 and inverted the
sign of a planted −0.050. A control that cannot fail is not a control.

**Chemical series** inflate absolute AUROC intervals by 1.63× but paired differences by only
1.04×, because both arms carry the same structure (388 series over 749 compounds, Kish effective
*n* 136.5). All interventions here are paired differences, so their intervals stand; any absolute
AUROC from a clustered panel requires the cluster interval. This is offered in place of a
scaffold-disjoint split, which answers a different question — generalisation to unseen chemotypes,
about which we claim nothing.

**The floor is only weakly level-specific.** Holding the panel, labels, class balance and score
multiset fixed and varying only which compound holds which score, then perturbing by resampling
each compound's own measured seed deviations, the induced sd peaks at AUROC 0.516 and falls away
— by only ~12% across 0.28 to 0.90. The simulation reproduces the real measurement at the real
panel (0.00487 against 0.00497, ratio 0.98), which is what makes the trend readable. A floor
measured where a method performs poorly is therefore conservative elsewhere, by a tenth rather
than a factor, and panel size dominates the level term by an order of magnitude.

## 4. Discussion

**A resolution limit is a filter, not a certificate.** Measuring it is cheap, the number is
predictable across panel sizes, and it correctly excludes effects that cannot be believed. It does
not, and cannot, tell you that an effect which clears it is real. Our largest intervention cleared
its floor by a factor of five and reversed sign on the next target we tried.

This is not an argument against measuring floors — it is an argument for the order of operations.
A result below the floor needs no further discussion. A result above it has earned the right to be
tested, and nothing more.

**The metric you report decides how much resolution you have.** EF@1% is the field's habitual
choice and the one this benchmark resolves worst, by a factor of 8.6 relative to its own scale. A
single-run EF@1% on a panel of this size carries irreproducibility of roughly a sixth of its
value. BEDROC is three times more resolvable and available at no cost.

**A floor must name its method, protocol, panel size and metric.** We made three substitutions and
each cost a factor: 4.0, 2.2 and 1.5. One further substitution we expected to be costly turned out
to be safe, which we could only establish by running the propagation rather than reasoning about
the ratio.

**Replication is the cheap experiment nobody runs on methodological claims.** Ours cost roughly
$8 of preemptible compute and one day, and it changed the status of every intervention in the
paper. It is difficult to argue that a methodological improvement is a property of the method
while declining to check it on a second protein.

### What this work does not show

Vina is below chance on Mpro and modest on Factor Xa; none of these interventions produces a
screen worth running, and SURVIVES throughout means "larger than the noise", never "useful". We
make no claim about generalisation to unseen chemical series. Two targets are two, and the
mechanism proposed for the gnina reversal — recovery from a broken baseline — is a hypothesis
consistent with the data, not a demonstrated cause.

## 5. Limitations, stated exhaustively

- **Two targets.** The floor's 1/*n* scaling is verified on both and within one; the intervention
  results are Mpro-measured and Factor-Xa-tested. Neither is a survey.
- **The fifth intervention is single-target.** Supplying a co-folding model with the correct
  binding site was not replicated: it needs GPU capacity we did not have, and its Mpro verdict was
  already INCONCLUSIVE, making it the least informative arm to add. It is excluded from the
  replication count, which is four, not five.
- **One replication arm ran on a truncated panel.** The unprotonated Factor Xa seeds hit a
  12-hour session limit at 875, 800 and 775 of 894 compounds. The surviving 772-compound panel is
  what finished first, not a random sample; its dropout control is −5.1pp against a 10pp
  threshold, which passes but is not negligible.
- **The Boltz-2 floor is not re-derivable from committed code.** It was computed inline;
  a later reconstruction reproduces it to 7% but not exactly. By our own standard that number is
  unverified, and it is flagged rather than presented as checked.
- **NEGLIGIBLE has no instance among the five interventions.** The framework's distinctive
  category is demonstrated only by planted controls.
- **The level-dependence result is a simulation**, validated against one real point. It has not
  been checked by re-docking a benchmark that genuinely reaches AUROC 0.8.
- **EF@1% and BEDROC floors are measured on two panels only**, and EF@1% in particular depends on
  the top-*k* cut interacting with panel size: 1% of 755 is eight compounds, and a larger panel
  would not be eight.
- **No systematic survey of reported effect sizes was conducted.** The 0.02–0.10 range in §1 is an
  impression; nothing here depends on it.

## Relation to the earlier version

This work is **version 3 of ChemRxiv record 10.26434/chemrxiv.15008496**
(`10.26434/chemrxiv.15008496/v3`) and supersedes **/v2**, posted 9 September 2026. It is a new
version of that record rather than a separate deposit.

**ChemRxiv registers no concept DOI for this record** — the base DOI does not resolve, only
`/v1` and `/v2` do — so a citation necessarily pins a version and cannot inherit a later
correction. That makes an explicit statement here, and a new version on the same record, the only
available means of supersession.

Two things in the earlier version do not survive, and both are corrected rather than quietly
dropped:

1. **Its two resolution limits measured the wrong quantity.** 0.020 was the disagreement between
   two *implementations* of a scoring function — implementation variance, not run-to-run
   reproducibility. 0.039 came from comparing two *search-effort conditions*, which is an
   intervention effect used as a noise floor; and because search effort was simultaneously one of
   the interventions being judged, that intervention was scored against itself. Both floors were
   the decision rule for every intervention in that version. This paper replaces them with seed
   replication of an unchanged protocol (§2.2), which is the measurement those numbers should
   have been.
2. **Its scoring-function result is now known not to replicate.** The +0.140 reported there stands
   as a correct Mpro measurement, but on a second target the same intervention returns −0.058
   (§3.6). The earlier version could not have known this; a reader of it should.

## Data and code availability

Every number in this paper is re-derivable from committed code and data at
**https://github.com/AegisMindApp/screening-decomposition**, which holds the pre-registrations,
the analysis scripts and the result files. Each result document names the script that produces it
and the JSON it writes. The one exception is stated in §5: the Boltz-2 floor of §3.3 was computed
inline and is reproduced to 7%, not exactly, by the committed reconstruction.

Pre-registrations are dated and were committed before the measurements they govern, including the
cross-target prediction of §3.2, the metric and level tests of §3.4 and §3.7, the replication
rules of §3.6, and the mechanism test whose INCONCLUSIVE result removed a claim from this version.

## References

1. Trott O, Olson AJ. AutoDock Vina: improving the speed and accuracy of docking with a new scoring function, efficient optimization, and multithreading. *J Comput Chem* 2010;31:455–461. doi:10.1002/jcc.21334
2. Eberhardt J, Santos-Martins D, Tillack AF, Forli S. AutoDock Vina 1.2.0: new docking methods, expanded force field, and Python bindings. *J Chem Inf Model* 2021;61:3891–3898. doi:10.1021/acs.jcim.1c00203
3. McNutt AT, Francoeur P, Aggarwal R, et al. GNINA 1.0: molecular docking with deep learning. *J Cheminform* 2021;13:43. doi:10.1186/s13321-021-00522-2
4. Butina D. Unsupervised data base clustering based on Daylight's fingerprint and Tanimoto similarity: a fast and automated way to cluster small and large data sets. *J Chem Inf Comput Sci* 1999;39:747–750. doi:10.1021/ci9803381
5. Lakens D. Equivalence tests: a practical primer for t-tests, correlations, and meta-analyses. *Soc Psychol Personal Sci* 2017;8:355–362. doi:10.1177/1948550617697177
6. Tran-Nguyen V-K, Jacquemard C, Rognan D. LIT-PCBA: an unbiased data set for machine learning and virtual screening. *J Chem Inf Model* 2020;60:4263–4273. doi:10.1021/acs.jcim.0c00155
7. Sieg J, Flachsenberg F, Rarey M. In need of bias control: evaluating chemical data for machine learning in structure-based virtual screening. *J Chem Inf Model* 2019;59:947–961. doi:10.1021/acs.jcim.8b00712
8. Passaro S, Corso G, Wohlwend J, et al. Boltz-2: towards accurate and efficient binding affinity prediction. *bioRxiv* 2025. doi:10.1101/2025.06.14.659707
