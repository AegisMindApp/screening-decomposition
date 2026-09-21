# Chemical-series structure: it matters for absolutes, not for the paper's differences

21 Sep 2026. Addresses the third remedy the JCIM editor named — "splits controlling
chemical-series effects". Expectation stated in `series_structure.py` before the run.

## The panel is strongly clustered

Butina clustering on Morgan fingerprints (radius 2, 2048 bits) at Tanimoto 0.65, over the 749
Mpro compounds with both a score and a parseable SMILES:

| | |
|---|---|
| series | 388 |
| singletons | 293 |
| largest series | 24 |
| top-5 sizes | 24, 17, 16, 16, 14 |
| **Kish effective n** | **136.5** against 749 compounds |

So the panel carries about **one fifth** of the independent information its compound count
suggests. Taken alone that looks fatal for intervals built by resampling compounds — which is
what `RESULT.md` records the manuscript doing.

## But it depends entirely on which quantity you are bootstrapping

| quantity | compound bootstrap | cluster bootstrap | widening |
|---|---|---|---|
| **paired difference** (exh32 − exh4) | width 0.04075 | width 0.04230 | **1.04x** |
| **absolute AUROC** (exh=4 arm) | width 0.1338 | width 0.2180 | **1.63x** |

A paired difference sees the *same* series in both arms, so the clustering cancels almost
entirely. An absolute AUROC has no such protection, and there the effective-n penalty shows up
in full: intervals are **63% too narrow** when compounds are resampled independently.

## What this means for the manuscript

**The five interventions are paired differences.** Their reported intervals are therefore very
close to correct — multiply widths by ~1.04, which changes no verdict:

| intervention | delta | CI | clears the exh=32 floor (0.00902 at n=755)? |
|---|---|---|---|
| scoring function (gnina CNN) | +0.140 | [+0.089, +0.191] | survives |
| pose ensemble | +0.055 | [+0.035, +0.075] | survives |
| receptor preparation repair | +0.045 | [+0.014, +0.076] | survives |
| binding site to Boltz-2 | −0.013 | [−0.039, +0.013] | does not |
| eight-fold search effort | +0.004 | [−0.006, +0.014] | does not |

**Any absolute AUROC the paper quotes needs the cluster interval**, not the compound one. That
is a real correction and it should be made rather than argued around.

## Why this is a better answer than a series-disjoint split

The obvious response to the editor would have been to re-run everything on scaffold-disjoint
splits. That would answer a question the paper does not ask. The paper's subject is the *size of
each intervention's effect on a fixed panel*, not generalisation to unseen chemotypes — and for
that quantity the pairing, not the split, is what controls series confounding. Measuring the two
bootstraps separately demonstrates the mechanism rather than asserting it.

A series-disjoint split remains necessary for any claim about a *new* chemical series. The paper
should state explicitly that it makes no such claim.

## Caveat

The comparison uses the 340 compounds common to the exh=4 and exh=32 arms, not all 749, because
the exh=32 stage ran a 350-compound panel. The series statistics are computed on the full 749.
Repeating the bootstrap comparison on a full-panel pair would be stronger; the mechanism — that
pairing cancels shared structure — does not depend on panel size.
