# Factor Xa: repair helps, barely — and the effect shrinks again on a third target

**2 September 2026.** 854 of 886 paired (96.4%, clears the 95% floor). 32 failures are the
known high-torsion tail. Same compounds, labels, box, binary and exhaustiveness 32 — protonation
the only change.

## Result

| | AUROC |
|---|---|
| unprotonated | **0.6641** (published on the full set: 0.6657) |
| protonated (0 → 403 polar H, 0 → 300 donor N) | **0.6775** |
| **ΔAUROC** | **+0.0134**, 95% CI **[+0.0009, +0.0266]** |

The CI clears zero — but only just: the lower bound is **+0.0009**. By the pre-registered rule
this reads as *preparation helped*, and it is the weakest such verdict the rule can return.

## The published Factor Xa conclusion is unchanged

0.6641 → 0.6775 is above chance before and after. **Factor Xa remains the one panel where docking
works and passed every gate.** Nothing on the site needs correcting; the number becomes slightly
better, not different in kind.

## Three targets now, and the effect shrinks each time

| target | ΔAUROC | 95% CI | mean score shift | top-10 overlap |
|---|---|---|---|---|
| **Mpro** | **+0.045** | [+0.024, +0.067] | 0.379 kcal/mol | 4/10 |
| **Factor Xa** | **+0.013** | [+0.001, +0.027] | 0.222 | 7/10 |
| **AcrB** | (no threshold crossed) | — | 0.128 | 9/10 |

A clean ordering, and it argues directly against extrapolating from Mpro. I earlier estimated
repair "would move Factor Xa to about 0.71". It moved it to **0.678** — my estimate was **2.4×
the measured effect**. The AcrB result had already warned that the Mpro magnitude does not
transfer; Factor Xa confirms it on a third target with a different site.

**Where hydrogen bonding matters most, the defect costs most.** Mpro's catalytic cleft depends
heavily on it; Factor Xa's S1 pocket less so; AcrB's large hydrophobic efflux channel least. That
is a plausible reading of the ordering, not a tested claim.

## What still holds

**The defect was real and worth repairing on every target.** But its practical consequence is
target-dependent and, on two of three targets, small. The honest summary of the whole
receptor-repair line is: *a genuine protocol fault, correctly fixed, that changed no published
conclusion anywhere.*

The one place it bites hardest is the shortlist rather than the summary statistic — 6 of 10
compounds changed on Mpro, 3 of 10 here. For a screening campaign that synthesises the top ten,
that is still the number to care about.
