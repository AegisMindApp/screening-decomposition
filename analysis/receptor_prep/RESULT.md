# Restoring the receptor's hydrogen-bond donors improves Mpro docking — and does not rescue it

**Run 31 August 2026.** Reading rule pre-registered at `1f68b20c4` before these numbers existed.
`analyse.py`, results in `RESULT.json`. **751/751 scored, zero failures.**

## What changed

Only protonation. Same 753 compounds, same labels, same box, same vina binary, exhaustiveness 4,
and the same heavy-atom coordinates — asserted as an identical multiset before the run.

| | unprotonated | protonated |
|---|---|---|
| polar hydrogens | 0 | 946 |
| donor nitrogens | 0 | **727** |
| acceptor nitrogens | 754 | 27 |
| receptor md5 | `d1bfc9ef` | `28657daa` |

## Result

| | AUROC |
|---|---|
| unprotonated (0 donors) | **0.4079** |
| protonated (727 donors) | **0.4530** |
| **ΔAUROC** | **+0.0451**, 95% CI **[+0.0236, +0.0671]** |

The CI lies **entirely above zero**. By the pre-registered rule that reads as:

> **Receptor preparation was a real part of the failure. Every docking number this project has
> published was computed with the defect and is affected.**

## What "affected" does and does not mean

It does **not** mean any published conclusion reverses. The effect is **+0.045 AUROC**. Applied to
our three targets it would move Mpro 0.427 → roughly 0.47 (**still below chance**), PD-L1 0.5948
→ roughly 0.64, Factor Xa 0.6657 → roughly 0.71. Every qualitative reading survives: Mpro still
fails, PD-L1 is still confounded, Factor Xa is still the clean panel.

What it does mean is that those numbers were produced by a protocol with a known fault, so each
is a lower bound on what a correctly prepared receptor would give, and none should be quoted
again without that caveat or without re-running.

**It also does not rescue docking.** 0.4530 is still below chance, and still 0.31 AUROC short of
the **0.763** seven descriptors achieve on the same compounds.

## What it changes in practice, which is more than the AUROC suggests

| | |
|---|---|
| Spearman ρ between the two score vectors | **+0.893** |
| mean absolute score change | **0.379 kcal/mol** |
| **top-10 shortlist overlap** | **4 / 10** |
| top-50 overlap | 35 / 50 |

**Six of the top ten compounds change.** For a screening campaign, where the deliverable is a
shortlist someone will synthesise, that is a large practical difference produced by a one-line
preparation bug — much larger than the 0.045 AUROC shift implies.

## Where this sits against the other two interventions

Three things have now been varied on this identical benchmark, each with a pre-registered rule:

| intervention | ΔAUROC | top-10 overlap |
|---|---|---|
| **scoring function** (Vina → gnina, identical poses) | **+0.140** [+0.091, +0.187] | — |
| **receptor protonation** (0 → 727 donors) | **+0.045** [+0.024, +0.067] | 4/10 |
| **8× search effort** (exh 4 → 32) | −0.019 [−0.032, −0.006] | 10/10 |

A coherent ordering falls out: **the scoring function dominates, receptor preparation matters
second, and search effort does not matter at all.** The two things practitioners most often spend
compute on — exhaustiveness and sampling — are the two that moved the answer least here.

None of it clears 0.763.

## What no outcome licenses

One target, one receptor, one library. Factor Xa and PD-L1 were prepared with the same defect and
each needs its own re-run; the +0.045 measured here does not transfer to them by assumption.

And the defect stands independently of this number: backbone amide nitrogens typed as
hydrogen-bond acceptors is wrong whether or not Vina's scoring is sensitive to it. The fix stays
in regardless.
