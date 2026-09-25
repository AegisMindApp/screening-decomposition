# The blind Mpro floor, measured at last — and it refutes a claim I put in the manuscript this morning

24 September 2026. Answers the amendment to `PREREGISTRATION_FLOOR_MPRO2.md`. The run completed
on Kaggle GPU with **70/70 scored, zero failures**, and its log records `conditioning: BLIND` —
the line added hours earlier precisely so a future reader could tell, working on the first run
that needed it.

Re-derive: `python analysis/boltz2/floor_propagate.py --pairs analysis/boltz2/mproblind/boltz2_scores.json --target Mpro --panel-ref 0.2024`

## The measurement passes the gate that rejected two attempts

| check | pocket-conditioned (voided) | **blind** |
|---|---|---|
| representativeness ratio vs the manuscript's panel | 1.263 **fail** | **0.995 pass** |
| bias against that panel | +0.0952 | **−0.00237** |
| correlation | 0.9746 | **0.9964** |
| per-compound sd | 0.02865 | **0.02341** |

Protocol and panel now match, so this is the floor the manuscript's Mpro row actually needs.

| propagated onto the 748-compound blind panel | sd | 95% bound |
|---|---|---|
| residual, level-aware | 0.00435 | 0.00474 |
| **margin, level-aware** (the quantity verdicts use) | 0.01007 | **0.01098** |

## The part that corrects me

`MPRO_FLOOR3_RESULT.md` measured Boltz-2's per-compound noise as target-specific — 2.71×
pocket-conditioned, and the blind run independently confirms it at **2.21×** (0.02341 against
Factor Xa's 0.01059). From that I wrote into §3.6 of the manuscript that *"a Boltz-2 floor
carried between targets is therefore wrong by more than the panel-size correction being applied
to it,"* and told the operator the Mpro row's borrowed floor was "probably an underestimate."

**That inference is false, and this measurement is what shows it.**

| | value |
|---|---|
| Factor Xa margin floor, rescaled 885 → 748 (what the manuscript charges) | 0.01107 |
| **Mpro margin floor, measured blind** | **0.01098** |
| difference | **0.8%** |

Borrowing the floor across targets was accurate to within one percent, despite the per-compound
noise differing by a factor of 2.2.

The mechanism is not mysterious once measured. AUROC is rank-based, so what matters is noise
**relative to the spread of the scores it has to reorder**, not noise in absolute units:

| target | score sd | per-compound noise | noise / spread |
|---|---|---|---|
| Factor Xa | 0.1872 | 0.01059 | 0.0566 |
| Mpro | 0.3305 | 0.02341 | 0.0708 |

Mpro's noise is 2.21× larger and its score distribution is 1.77× wider, so in the units the
metric responds to the two targets differ by 1.25×, not 2.2× — and the propagated floors land on
top of each other.

**The error I made was quoting a quantity outside the scope of its definition**: a per-compound
noise ratio is a statement about scores, and I used it as a statement about an AUROC floor
without checking the propagation that connects them. The propagation is the check, it was
available, and I asserted before running it.

## What now stands

- **Per-compound Boltz-2 noise is target-specific.** Measured twice, 2.21× and 2.71×.
- **The propagated AUROC-margin floor is not**, at least between these two targets, and the
  panel-size rescaling alone gets within 1%. Two targets is two, and a pair of panels whose
  relative noise happens to be close is weak evidence that this holds generally — but it is
  the evidence there is, and it points the opposite way to what the manuscript said.
- **§3.3's Boltz-2 row is unchanged in verdict.** Its interval [−0.039, +0.013] straddles a floor
  of 0.01098 exactly as it straddled 0.01106. INCONCLUSIVE, now against a floor that was measured
  on the right target under the right conditioning rather than borrowed.

The manuscript is corrected at both places: §3.6's inference is withdrawn and replaced with what
was measured, and §3.3 now charges the measured floor and drops the "probably too small" note.
