# Can we build a better docking algorithm? The best idea in our own data fails a control

**2 September 2026.** All numbers below are on the **repaired (protonated) receptor**, with
paired-bootstrap 95% CIs over the same compounds. An earlier version of this file reported
these on the donor-defective receptor and without CIs; both are corrected here.

## 1. What docking is worth over seven free descriptors

Descriptors: MW, clogP, HBD, HBA, RotB, TPSA, formal charge. Out-of-fold logistic
regression, StratifiedKFold(5), identical folds per arm.

| target | Vina alone | descriptors | desc + Vina | **docking adds** | 95% CI |
|---|---|---|---|---|---|
| Mpro (751, 257 act) | 0.4530 | 0.7653 | 0.7630 | **−0.0023** | [−0.0044, −0.0002] **sig.** |
| Factor Xa (854, 388 act) | 0.6775 | 0.7129 | 0.7279 | **+0.0150** | [−0.0033, +0.0339] **n.s.** |

**On neither target does docking add anything demonstrable over descriptors that cost
microseconds.** On Mpro it significantly *subtracts*. Factor Xa's +0.015 — the number the
previous version of this file published as +0.011 without an interval — **crosses zero**.

## 2. The one algorithmic idea our data actually motivated, and why it fails

Docking scores decompose as `score = f(compound) + g(compound, target)`. Descriptors explain
52–68% of the score, and the AcrB/CTSS control put cross-target correlation at r = 0.932 —
so `f` dominates. The natural derivative algorithm: **rank on `g`, the residual, discarding
the compound-intrinsic part.**

| target | raw AUROC | residual AUROC | **residual − raw** | 95% CI |
|---|---|---|---|---|
| Mpro | 0.4530 | 0.5300 | **+0.0770** | [+0.0353, +0.1189] **sig.** |
| Factor Xa | 0.6775 | 0.5589 | **−0.1186** | [−0.1486, −0.0903] **sig.** |

**Opposite directions, both significant.** Residual-ranking rescues the target where docking
is broken and destroys the target where docking works. Since you cannot know in advance which
regime a new target is in, this is a **symptom detector, not a method**. Idea killed by its
own control.

## 3. The measurement floor — the sharpest thing here, and it is ours

The gnina rescore worker embedded a control nobody had run: gnina's `--score_only` Vina term
should reproduce the cached docking score if the poses and receptor are what we think.

    Vina docking cache vs gnina --score_only:  r = 0.9369
    mean difference −0.043 kcal/mol; disagreements >2 kcal/mol: 0/745
    AUROC from cache 0.4182   AUROC from gnina 0.3793

The control **passes** on pose identity — and yet two scores agreeing at r = 0.937, never
differing by 2 kcal/mol, produce **AUROCs 0.039 apart**.

That is the finding with real teeth. A benchmark of this size cannot resolve a scoring
improvement smaller than about 0.04 AUROC, because sub-kcal numerical noise already moves the
metric that far. Most published docking improvements are reported in that range. **A genuinely
better scorer and a rounding difference are not distinguishable here** — which is the real
reason we cannot build and validate a better algorithm, and it has nothing to do with our
budget.

## 4. gnina's advantage, stated only as far as the evidence reaches

Same tool, same poses, same receptor, so the comparison is internally clean:

    gnina Vina term  raw 0.3793   residual 0.5070
    gnina CNNaffinity raw 0.5389   residual 0.5084   (raw +0.16, residual +0.0014)

**This is not a replication of Chen et al. 2019** and must not be written as one. Chen ablated
the *receptor* and got R² = 0.998 — a direct test of whether receptor information was used. We
residualised against seven ligand descriptors, which is a different and weaker test: a
physically *correct* scorer would also have a large descriptor-explainable component, because
binding free energy genuinely covaries with size and lipophilicity. That is why ligand
efficiency exists. Residualising strips a perfect scorer too.

The defensible claim is only: *the part of these docking scores orthogonal to seven
descriptors does not rank actives in these sets.* Not "docking ignores the target."

## 5. Prior art — the critique space is occupied

- Chen et al. 2019, PLOS One — DUD-E analogue/decoy bias; CNN ignores the receptor (R² = 0.998).
- Sieg, Flachsenberg & Rarey 2019 — descriptor baselines as required bias control.
- 2026, arXiv 2507.21404 — data leakage and redundancy in **LIT-PCBA**, DUD-E's bias-corrected
  successor.
- 2026, arXiv 2605.01681 — no docking method dominates on LIT-PCBA; only modest early enrichment.

The field has known the evaluation substrate is broken for seven years, and the replacement
benchmark is broken too.

## Conclusion

**No, we cannot construct a groundbreaking docking algorithm — and the binding constraint is
measurement, not method or budget.** The best idea our own data supports (residual ranking)
reverses sign between two targets. The benchmark cannot resolve improvements below ~0.04
AUROC. And the critique that would otherwise be our contribution was published in 2019 and
extended in 2026.

What survives from this week is a set of controls that each *could have failed and did*: the
receptor donor defect, the size-selected dropout that would have biased arm 2, the
exhaustiveness null, the cross-target r = 0.932, the pose-identity r = 0.937. That is a
capability. It is not a docking algorithm, and it should not be dressed up as one.
