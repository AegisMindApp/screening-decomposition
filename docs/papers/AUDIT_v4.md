# Audit of `resolution_v4_draft.md`, 25 September 2026

Run with the `paper-audit` skill. Four defects found, all fixed; verdicts unchanged throughout.

## Step 0 — cited artifacts

No supplementary tables, figures or appendices are cited, so nothing can be missing. All five
internal cross-references (§1, §3.3, §3.5, §3.6, §3.7) resolve.

**BLOCK (fixed):** the draft asserted "data in the public repository" and never named it. A
**Data and code availability** section now points at
`https://github.com/AegisMindApp/screening-decomposition` and states the one exception — the
Boltz-2 floor, reproduced to 7% rather than exactly.

## Step 1 — numbers

Twenty-two headline numbers checked against their committed source files.

| group | result |
|---|---|
| three floor rows (mean, sd, bound) | **VERIFIED**, exact |
| per-compound sd and rank span, both protocols | **VERIFIED** |
| eight metric-scope figures (§3.4) incl. CV ratio 8.6× and 15.4× | **VERIFIED** |
| four replication deltas and CIs | **VERIFIED** |
| gnina absolute AUROCs and control r = 0.9934 | **VERIFIED** |
| mechanism-test slopes and sub-panel count | **VERIFIED** |
| **four replication floor values** | **MISMATCH — fixed, see below** |

**BLOCK (fixed): the paper and its own code disagreed on every replication floor.** The draft
rescaled Factor Xa's **n = 12** bound (0.00682); the committed scripts hardcoded the superseded
**n = 10** value (0.00703). Neither was wrong arithmetic — they were different vintages of the
same measurement, and a reader re-running the scripts would not have reproduced the table.

Fixed at the source: the three replication scripts now share one documented constant
`FXA_EXH4_FLOOR = 0.00682` with its provenance, and the paper carries what they emit —
0.00697 (gnina), 0.00682 (pose ensemble), 0.00731 (receptor prep), 0.00517 (search effort,
measured from its own seeds and unaffected). **All four verdicts are unchanged.**

*In fixing this I first aligned the paper to the old script values, then updated the scripts —
briefly making the mismatch worse in the other direction. Recorded because it is the same
drift-between-two-copies failure the audit exists to catch.*

## Step 1b — plausibility

All AUROCs in [0, 1]; every χ² penalty recomputes correctly (1.460× at n = 15, 1.645× at 10,
1.551× at 12); EF@1% CV confirmed as sd/mean, not bound/mean; "a sixth" and "a third" both check
out; 1% of 755 is 8 compounds as stated. Vina is **below chance** on Mpro (0.4163) and this is
disclosed in §2.1 rather than buried.

**BLOCK (fixed): a repurposed magnitude.** The abstract and §3.6 both said gnina's +0.140 was
"five times the floor". Against the current Mpro exh-4 floor of 0.00648 it is **twenty-two
times**, with even its 95% lower bound fourteen times the floor. The figure five matched no
current quantity — it traces to the superseded borrowed floor of 0.0201. The error *understated*
the paper's own case, which is why it survived several readings.

## Step 3b — front matter against body

All eleven substantive abstract claims are supported in the body.

**BLOCK (fixed): one dropped qualifier.** The abstract read "Two mid-sized effects shrink to
roughly a quarter and cannot be resolved". Only receptor prep shrinks (+0.0451 → +0.0116). Pose
ensemble **changes sign** (+0.0552 → −0.0135); it is a quarter in magnitude only. The body carries
"one pointing the wrong way" and the abstract dropped it. Rewritten to distinguish the two.

Also removed in the same sentence: "resolved beyond any reasonable doubt", an unhedged absolute
now replaced by the interval that supports it.

Qualifier carry-through is otherwise intact — the abstract says "simulated panels" for the level
result and "the four testable ones" for the replication, flagging the fifth intervention.

## Step 3 — spine

> **This paper claims: clearing a measured resolution floor is necessary and not sufficient.**

Every substantive section serves it. §2.1 (panels and protocol) and §4's "What this work does not
show" are descriptive and scope-limiting respectively, which is appropriate rather than a break.

**Three sharpest objections, and where each is answered:**

1. *"Your docking is below chance — why should I care?"* — §2.1 states it, §3.5 says SURVIVES
   means larger than the noise and never useful, §5 repeats it. Answered.
2. *"Two targets is not a generalisation."* — §5 first bullet. Answered, and the paper does not
   claim otherwise.
3. *"Why did the largest effect reverse?"* — §3.6 tested the obvious explanation and withdrew it
   as INCONCLUSIVE under a pre-registered rule. Answered by declining to answer, which is the
   defensible position.

## Step 4 — prose

No throat-clearing openings, no runs of uniform sentence length, no unhedged absolutes remaining
after the §3b fix.

## Gate

```
Artifacts:    0 cited | 0 missing            PASS (availability statement added)
Numbers:      22 verified | 0 untraced | 0 contradicted | 4 mismatches FIXED
Plausibility: 1 repurposed magnitude FIXED   PASS
Front matter: 1 dropped qualifier FIXED      PASS
Spine:        clear, all sections serve it   PASS
Prose:        0 flags                        PASS
Reproduction: paper and committed scripts now agree on all four floors  PASS
References:   8 cited, 8 listed, 0 orphans   /verify-refs NOT re-run since v3
```

**READY TO SUBMIT: not yet — one item outstanding.**

`/verify-refs` has not been run against v4. No reference was added and no claim about a cited
work changed from v3, so the risk is low, but the claim-provenance step has not been repeated
since §3.2, §3.3, §3.6 and §3.7 were rewritten. That is the remaining gate.
