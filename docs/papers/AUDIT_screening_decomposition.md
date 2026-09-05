# Pre-submission audit — `screening_decomposition_preprint.md`

**6 September 2026.** `/paper-audit` against commit `ea9f00a2e`. Verdict: **NOT READY.**
Four blockers, all fixable without new experiments except one.

---

## Step 0 — Cited artifacts

| Cited artifact | Line | Exists? | Verdict |
|---|---|---|---|
| "commit hashes are cited in the **supplementary material**" | 82 | **NO** — no supplementary file exists | **BLOCK** |
| "committed to a **public** git repository" | 81 | **NO** — `tradingjohn/aegismind` is `PRIVATE` (`gh repo view`) | **BLOCK** |
| "All reading rules... are in the project repository" (§5) | 194 | no URL, no DOI, unreachable | **BLOCK** |
| `PREREGISTRATION.md` `9e3083d02`, `PREREGISTRATION_ARM2.md` `e96d4fea1`, `PREREGISTRATION_TARGETS.md` `68bd9b94a`, `analyse_boltz2.py` `5afb5bec8` | various | exist locally, unreachable publicly | see below |

**The private-repo finding is a validity problem, not a formatting one.** The paper's stated
contribution is that reading rules were fixed before the data existed. The discriminating
question is: *can a reader confirm `9e3083d02` predates the data?* Today, no. A Zenodo deposit
made now does not fix this — it timestamps the deposit, not the commit.

Two acceptable resolutions:
1. Make the pre-registration artifacts publicly inspectable **with history intact**, then the
   claim stands as written.
2. Keep them private, drop the word "public", and state precisely what a reader can and cannot
   verify. The pre-registration then becomes an author assertion, and the paper must say so.

Before making anything public, run the standing patent-disclosure check against the two
surviving filings (AU2026905785, AU2026906776).

---

## Step 1 — Numbers audit

### VERIFIED (recomputed or read from source)

| Claim | Value | Source | Status |
|---|---|---|---|
| Vina/gnina score agreement | r = 0.9369 | `docking_value/MARGINAL_VALUE.md:§3` | VERIFIED |
| mean difference | −0.043 kcal/mol | same | VERIFIED |
| disagreements > 2 kcal/mol | 0 / **745** | same | VERIFIED (paper says "no compound"; n=745 not stated) |
| AUROC cache / gnina | 0.4182 / 0.3793 | same | VERIFIED |
| gap | 0.039 | derived | VERIFIED |
| Scoring function ΔAUROC | +0.140 [+0.091, +0.187] | `three_arm_docking/ARM3_RESULT.md:27` (+0.1397) | VERIFIED |
| Receptor repair ΔAUROC | +0.045 [+0.024, +0.067] | `receptor_prep/RESULT.md:24` (+0.0451) | VERIFIED |
| "six of the top ten changed" | 4/10 overlap | `receptor_prep/RESULT.md:47` | VERIFIED |
| 754 nitrogens typed acceptor, 0 polar H | 754 / 0 | `receptor_prep/RESULT.md:12-16` | VERIFIED |
| Pose ensemble ΔAUROC | +0.019 [−0.017, +0.055] | `pose_ensemble/RESULT.md:22` (+0.0194) | VERIFIED |
| DiffDock-L ΔAUROC | +0.017 [−0.027, +0.064] | `three_arm_docking/ARM2_RESULT.md:24` (+0.0173) | VERIFIED |
| Pose-source correlation | r = 0.139 | `three_arm_docking/ARM2_RESULT.md:38` | VERIFIED |
| Pocket conditioning | −0.013 [−0.026, +0.000] | `boltz2/ARM2_RESULT.md` (−0.0130) | VERIFIED |
| Search effort | −0.019 | `receptor_prep/RESULT.md:66` | VERIFIED — **but CI omitted, see below** |
| Mpro Vina / desc / Boltz-2 / combined | 0.4530 / 0.7654 / 0.7913 / 0.8588 | `boltz2/RESULT.md` | VERIFIED |
| FXa Vina / desc / Boltz-2 / combined | 0.6775 / 0.7041 / 0.7227 / 0.7791 | `boltz2/FXA_RESULT.md` | VERIFIED |
| Mpro combination gain | +0.0934 [+0.0616, +0.1253] | `boltz2/RESULT.md` | VERIFIED |
| FXa combination gain | +0.0751 [+0.051, +0.099] | `boltz2/FXA_RESULT.md` | VERIFIED |
| Standalone gains | +0.0259 / +0.0186, both spanning zero | same | VERIFIED |
| Arm1↔arm2 correlation | r = 0.952 | `boltz2/ARM2_RESULT.md` | VERIFIED |
| descriptors explain 63.1% / 62.8% | same | `boltz2/SIZE_BIAS.md:9-10` | VERIFIED |
| Boltz-2 residual 0.6418, Vina 0.5382 | same | `boltz2/SIZE_BIAS.md:9-10` | VERIFIED |
| LIT-PCBA full median | **0.7260** | `litpcba/descriptor_baseline_full.json` — recomputed, 8th of 15 | VERIFIED |
| LIT-PCBA full ≥50 actives | **0.6830** | recomputed, 10 targets | VERIFIED |
| LIT-PCBA AVE median | **0.6510** | `descriptor_baseline_ave.json` — recomputed | VERIFIED |
| AVE robustness 0.6510 / 0.6497 / 0.6510 | ≥0 / ≥20 / ≥50 | recomputed | VERIFIED |
| Panel AUROC range | 0.577–0.864 | recomputed | VERIFIED |
| Beats Vina by 0.116; 0.075 property bias | derived | arithmetic | VERIFIED |

### CONTRADICTED / UNTRACED — these block

**1. The residual band "0.507 to 0.555" — CONTRADICTED and method-mixed.**

`0.555` appears nowhere in the repository except `SIZE_BIAS.md`'s own prose. The measured
docking residuals actually on disk are:

| score | residual AUROC | source | residualisation |
|---|---|---|---|
| gnina Vina term (Mpro) | 0.5070 | `MARGINAL_VALUE.md:61` | OOF logistic |
| gnina CNNaffinity (Mpro) | 0.5084 | `MARGINAL_VALUE.md:62` | OOF logistic |
| Vina (Mpro) | 0.5300 | `MARGINAL_VALUE.md:30` | OOF logistic |
| Vina repaired (Mpro) | 0.5382 | `SIZE_BIAS.md:9` | OOF **gradient-boosted** |
| **Vina (Factor Xa)** | **0.5589** | `MARGINAL_VALUE.md:31` | OOF logistic |

Two problems. **(a)** The Factor Xa value, 0.5589, lies **outside** the stated band. **(b)** The
band aggregates two different residualisation procedures, and §3.4's punchline compares a
*gradient-boosted* Boltz-2 residual (0.6418) against that mixed set.

**The conclusion survives; the interval does not.** The honest band from artifacts in hand is
0.507–0.559, and 0.6418 sits clearly outside it. Fix: recompute all five docking residuals under
one procedure (zero compute — the scores exist), or narrow the claim to the gradient-boosted set.

**2. Boltz-2's residual is Mpro arm 1 only, compared against a band claimed "across both targets."**
Abstract line 38 and §3.4 read as a matched comparison. There is no Factor Xa Boltz-2 residual.
Either compute it or say the residual analysis is single-target.

### Presentation errors

| Claim | Paper | Source | Fix |
|---|---|---|---|
| Mpro active fraction | 34.3% | 257/751 = **34.22%** | 34.2% |
| FXa active fraction | 45.8% | 405/886 = **45.71%** | 45.7% |
| "4,000 resamples" (§2, all statistics) | 4,000 | **10,000** for §3.2 arms (`receptor_prep/PREREGISTRATION.md:39`, `three_arm_docking/PREREGISTRATION.md:48`); 4,000 only for Boltz-2 (`analyse_boltz2.py:71`) | state both |
| "exhaustiveness 32" (§2, all docking) | 32 | **4** for receptor-prep and pose-ensemble arms; 4→32 *is* the search-effort arm | state per-arm |
| Search-effort CI | "top-ten unchanged" in CI column | **[−0.032, −0.006]** — excludes zero | print the CI; it is a significant small *negative*, and hiding it in the CI column obscures that |
| FXa panel size | "854–886 compounds" | two different panels: 854 (`MARGINAL_VALUE.md`, desc 0.7129) and 886 (`FXA_RESULT.md`, desc 0.7041) | see Step 3b |

---

## Step 1b — Plausibility

| Check | Result |
|---|---|
| Physically possible? | Yes — all values are AUROCs in [0,1] |
| Beats trivial baseline? | Yes, and the paper's most useful move is that several results *don't* — permutation controls landed 0.4952 and 0.5165, both inside [0.45, 0.55] |
| Units repurposed? | **One risk.** "0.04" serves as both the *measured* gnina/Vina gap and a *pre-registered threshold* in five separate arms. Same number, two roles |
| Inside training range? | Not applicable to AUROC, but see Boltz-2 scope note — Mpro scaffolds are likely in Boltz-2's training data; the paper says this only in §3.4's parent artifact, not the paper |

**The floor is n = 1.** This is the paper's lead contribution and rests on one observation of one
score pair. It is a valid *existence* argument — sub-kcal noise can move AUROC by 0.039 — but the
paper then uses 0.04 to gate six other results as if it were calibrated. Same failure shape as
the debias gate in this project's own record, which was used as an admissibility criterion for a
month before anyone measured its null.

**Cheap fix, no new compute:** paired-bootstrap the difference (AUROC_cache − AUROC_gnina) and
report the floor with an interval. That converts the headline from an anecdote into a
measurement. Also state which receptor and how many compounds it was measured on — 0.4182 matches
neither the unprotonated 0.4079 nor the repaired 0.4530, and the source records 0/745.

---

## Step 2 — References

Run `/verify-refs`. Two references carry the title and must be verified against primary sources,
not search summaries:

- **[4] arXiv:2605.01681, "Vina 0.61"** — confirm whether 0.61 is median or mean AUROC, over
  which 15 targets, against which inactive sets.
- **[5] Sunseri & Koes 2021, "GNINA 0.61–0.62"** — no journal or volume given. Likely *Molecules*
  2021, "Virtual Screening with Gnina 1.0"; confirm the exact table and metric.

**Comparability gap, unstated:** inactives were capped at **25,000 per target by seeded random
sampling**; the published figures are on LIT-PCBA's full inactive sets (hundreds of thousands).
AUROC is robust to that in expectation, but it needs one explicit sentence, and this project's own
record contains a case where gate admissibility flipped both ways under active subsampling.

---

## Step 3 — Spine

> **This paper claims:** most reported virtual-screening improvements are smaller than the
> benchmark can resolve, and when six interventions are measured against a measured resolution
> limit, only changing *what does the scoring* survives.

| Section | What it does | Serves spine? | Issue |
|---|---|---|---|
| Abstract | Floor, then decomposition, then LIT-PCBA, then Boltz-2 | Yes | overstates on two counts — Step 3b |
| §1 Introduction | Positions against known-bias literature | Yes | — |
| §2 Methods | Panels, descriptors, docking, Boltz-2, prereg, stats | **Partial** | resamples and exhaustiveness wrong; no repo URL |
| §3.1 Floor | Establishes the limit | Yes — this is the spine | n=1, no interval |
| §3.2 Decomposition | Six interventions vs the floor | Yes | — |
| §3.3 LIT-PCBA | External validation of descriptor baseline | Yes | comparability caveat missing |
| §3.4 Boltz-2 | The one thing that clears | Yes | residual band; single-target residual |
| §4 Limitations | Honest, including a withdrawn caveat | Yes | strongest section in the paper |
| §5 Availability | Points at an unreachable repository | **No** | BLOCK |

### The three sharpest objections

**Objection 1 — "Your floor is one observation, and you use it as a threshold six times."**
*Answered?* No. §3.1 asserts it; §4 says it should be re-measured elsewhere but not that it is
itself uncertain here.
*Insert at end of §3.1:* "The 0.039 gap is a point estimate from one score pair. A paired
bootstrap over compounds places it at [X, Y], which is the interval we use when calling a
difference unresolvable. The floor is an existence result — sub-kcal numerical noise *can* move
AUROC this far — not a calibrated constant, and we apply it as a screening criterion rather than
a test."

**Objection 2 — "Your descriptors used 25,000 sampled inactives; the published docking numbers
used the full sets. That is not like for like."**
*Answered?* No.
*Insert at end of §3.3 paragraph 1:* "Our descriptor baselines use up to 25,000 inactives per
target, drawn by seeded random sampling, against published docking figures computed on the full
inactive sets. AUROC is insensitive to the negative-class size in expectation, so we expect no
systematic direction to this difference, but it is not a matched comparison and we do not claim
one."

**Objection 3 — "You residualise against seven ligand descriptors and conclude docking carries no
target-specific signal. A physically correct scorer would also have a large descriptor-explainable
component — that is why ligand efficiency exists."**
*Answered?* No — and `MARGINAL_VALUE.md:§4` makes exactly this point about the project's own work,
so the paper is currently weaker than its source artifact.
*Insert after §3.4 paragraph 2:* "Residualising against ligand descriptors is a weaker test than
receptor ablation: binding free energy genuinely covaries with size and lipophilicity, so a
correct scorer would also lose signal under this operation. The defensible claim is narrow — the
part of these docking scores orthogonal to seven descriptors does not rank actives in these sets —
not that docking ignores the target."

---

## Step 3b — Front matter vs body

| Front matter says | Body / source says | Verdict |
|---|---|---|
| "docking **loses** on both targets we assembled (… Factor Xa 0.6775 vs 0.7041)" (l. 26) | Paired marginal value on FXa is **+0.0150 [−0.0033, +0.0339], n.s.** (`MARGINAL_VALUE.md:16`) | **OVERSTATED** — by the paper's own §2 rule ("demonstrated only when the interval excludes zero"), this is supported on one of two targets |
| Same sentence, "0.6775 vs 0.7041" | 0.6775 is from the **854**-compound panel where descriptors score **0.7129**; 0.7041 is from the **886**-compound panel | **MISMATCHED** — two panels in one comparison. Like-for-like is 0.6775 vs 0.7129 |
| "every docking score we measured sat between 0.507 and 0.555" (l. 38) | max measured is **0.5589** | **CONTRADICTED** |
| Title: "seven descriptors that **beat docking** on LIT-PCBA" | True on the full set only; §3.3 says so explicitly and §4 repeats it | **acceptable but load-bearing** — depends entirely on refs [4] and [5] verifying |
| "commit hashes are cited in the supplementary material" | no supplementary material | **UNSUPPORTED** |
| "public git repository" | private | **FALSE** |

Suggested abstract repair for the first row:

> "Against seven free physicochemical descriptors, docking adds nothing demonstrable on either
> target we assembled: on Mpro it significantly subtracts (−0.002 [−0.004, −0.000]), and on Factor
> Xa its +0.015 [−0.003, +0.034] crosses zero."

That is stronger than the current wording, not weaker — a significant *subtraction* is a sharper
result than "loses", and it survives the paper's own reading rule.

---

## Step 3c — Corrections

The paper carries one withdrawn caveat (§4, "A caveat we withdrew" — the arm-1 lower-bound
claim). Checked: the retracted "lower bound" framing appears **only** in §4 where it is struck.
The abstract and §3.4 state +0.093 without the lower-bound qualifier. **Correctly landed.**

Note that `analysis/boltz2/RESULT.md` still asserts the lower-bound framing as current in its
"What this does not show" section. That is an internal artifact, not the submission, but it should
carry a pointer to `ARM2_RESULT.md` so a reader of the repository does not read a withdrawn
caveat as live.

---

## Step 4 — Prose register

Low flag count. The paper is written in a specific register with short declaratives breaking up
longer sentences, and it does not read as generated.

| Flag | Location | Issue | Fix |
|---|---|---|---|
| **B** unhedged confidence | l. 168–169: "The ceiling we observed for docking is a property of the scoring approach, not of the task." | A general claim from two targets and one co-folding model | "On these two targets, the ceiling is a property of the scoring approach rather than of the task — a third scoring approach clearing it would be the test that generalises this." |
| **B** unhedged confidence | l. 116–117: "essentially unrelated coordinates — and produce the same chance-level ranking" | r = 0.139 is low but "essentially unrelated" overstates it, and "the same" hides ΔAUROC = +0.017 | "nearly uncorrelated coordinates — and produce rankings that differ by less than the floor" |
| **A** uniform density | l. 44–48, four consecutive parallel-structure sentences ("Two properties…", "Without the first…", "Without the second…") | rhythm is close to templated | break the parallelism on the third |
| **D** throat-clearing | l. 55: "This paper reports three things." | acceptable in a short preprint; borderline | keep |

No **C** (adjacent-information) or **E** (passive-in-Results) flags. §4 Limitations is the best
writing in the paper — particularly "A caveat we withdrew", which reports a pre-registered caveat
that turned out to be wrong. Keep that section intact.

---

## SUBMISSION GATE

```
Artifacts:    4 cited | 3 unreachable or missing              — BLOCK
Numbers:      31 verified | 0 estimates | 1 UNTRACED (0.555)  — BLOCK
              | 1 CONTRADICTED (FXa residual 0.5589)          — BLOCK
              | 6 presentation errors (34.3%, 45.8%, resamples,
                exhaustiveness, omitted CI, panel size)
Plausibility: possible YES | beats trivial baseline YES
              | units repurposed — 0.04 is both a measurement
                and a threshold; floor is n=1 with no interval — FIX
References:   /verify-refs run? NOT YET — refs [4] and [5] carry the title
Spine:        Clear — 1 section (§5) without a spine connection
Front matter: 4 claims overstating or contradicting the body   — BLOCK
Corrections:  every instance struck at the point of claim? YES
Output:       markdown is the source of truth; no PDF built yet
Prose:        4 flags (A×1, B×2, D×1) — none blocking

READY TO SUBMIT: No
```

### Blockers, in order

1. **"Public git repository" is false.** Either publish the pre-registration artifacts with
   history intact (after the patent-disclosure check), or drop "public" and state what a reader
   can actually verify. The paper's whole contribution rests on this.
2. **Supplementary material is cited and does not exist.** Produce it — it is a short table of
   arm → pre-registration file → commit hash → date — or cut the sentence.
3. **The residual band 0.507–0.555 is contradicted by a measured 0.5589** and mixes two
   residualisation procedures. Recompute under one procedure, or narrow the claim.
4. **"Docking loses on both targets" fails the paper's own reading rule** and pairs numbers from
   two different Factor Xa panels.

### Not blocking, but fix before submission

5. Bootstrap the resolution floor and give it an interval; state its receptor and n.
6. Add the 25,000-inactive comparability sentence to §3.3.
7. Correct 34.2%, 45.7%, per-arm resample counts, per-arm exhaustiveness, and print the
   search-effort CI.
8. Insert the three objection answers.
9. Either compute a Factor Xa Boltz-2 residual or state that the residual analysis is
   single-target.
