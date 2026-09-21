# ALDH1 Boltz-2: what we will do at each cost, fixed before the cost is known

Written 19 Sep 2026 at 23:45 local, while the L4 spot probe (`boltz2-aldh1-probe`,
us-central1-a) is still in stage 1. **No projection exists yet.** That is the point: the
branches below are drawn while the number cannot bias where they fall.

The failure this prevents is the one recorded in `ALDH1_BLOCKED.md` and complained about
directly — a spending judgment made ad hoc after seeing the number, which is hand-feeding the
process rather than running it.

## The number we are waiting for

The worker runs unsharded, so `BOLTZ_NSHARDS=1`, `names` is all 750 and `_n_chunks = 30`. Its
`projected_gpu_hours` is therefore the **whole-panel** cost, not a per-shard cost. Sharding 4×
gives 8 chunks × 4 = 32 fixed charges against 30, which is inside the noise of the estimate.

Two corrections are applied to it before any branch is read:

| correction | factor | why |
|---|---|---|
| probe-ligand size | **×1.06** | The probe is `ALDH1_0000`–`0011`, all label 1, mean 22.9 heavy atoms against the panel's 24.3 — 5.7% small. Boltz-2 cost is dominated by the 988-residue modelled dimer, not the ligand, so treating runtime as linear in heavy atoms **overstates** the correction. Deliberately conservative. |
| per-instance setup `S` | **+S per instance** | `dtA` starts *after* pip install and weight download, so the worker's projection excludes boot, install, model download and MSA. On one instance that is one charge; on four shards it is paid four times. Measured from this run's own log. |

So the budget is `H = 1.06 × projected` GPU-hours of prediction, plus `S` per instance.

Planning rate: **$0.35 / instance-hour** for g2-standard-8 spot in us-central1 — a deliberately
high estimate against a ~$0.85/hr on-demand list. Reconciled against actual billing afterwards;
the rate is not allowed to move a branch after the fact.

## Branches

Let `T = H + (instances × S)` be total billed instance-hours.

| | condition | action |
|---|---|---|
| **A** | `T ≤ 12 h` (≲ $4.20) | Run the full 750 on **one** spot instance, unsharded, checkpointing to `gs://aegismind-boltz2-1788402939/aldh1`. Proceed without asking — this is inside the cost John already sanctioned. |
| **B** | `12 < T ≤ 30 h` (≲ $10.50) | Run the full 750 **sharded 4×**, paying `S` four times, wall clock ≈ `H/4 + S`. Proceed. Report the cost when it lands. |
| **C** | `30 < T ≤ 60 h` (≲ $21) | **Stop and put the number to John.** Do not spend. Present branch D's options alongside it. |
| **D** | `T > 60 h` | Do not run the 750. Either shrink the panel under the constraint below, or stop at two admissible targets and leave Boltz-2 `INSUFFICIENT`. |

## The constraint on every shrink branch

`THIRD_TARGET.md` pre-commits: *"Admissibility is re-checked on the final panel and recorded.
If the built panel's descriptor baseline exceeds 0.80 the target is dropped and this document is
amended — the check is not a formality it has already passed."*

The 0.5698 baseline was measured on **300 actives + 450 inactives**. It does **not** carry to a
350- or 120-compound subset. So on any shrink:

1. Re-run `harness.target_admissible` on the **actual** subset and record the baseline.
2. Amend `THIRD_TARGET.md` with the new panel and the new baseline, **before** spending.
3. Quantify the precision cost of the smaller panel the same way the exh=32 panel was sized —
   `ALDH1_BLOCKED.md` option 1 already says this should happen *before* choosing, not after.

Substituting a baseline measured on one panel for a different panel is exactly the substitution
`docs/papers/RESUBMISSION_GATE.md` blocks on for exh=4 → exh=32. It does not get a pass here
because the number happened to clear once.

## A lever to check before branch C or D

If the log reports `msa reuse: NONE (server)`, every shard regenerates the MSA for the same
protein. Staging one MSA to GCS and reusing it across shards removes that charge from `S`
entirely. Size it from the log before concluding the panel is unaffordable.

## Housekeeping fixed now, not later

The instance is created with `--instance-termination-action=STOP`. A stopped g2 keeps its
**150 GB pd-balanced** boot disk, which is roughly **$0.50/day** — not the ~$0.02/day figure I
first assumed, a 25× error. Forgotten paid resources are already a recorded failure mode here
(`project_tpu_idle_quota`: 16 idle v5e chips at ~AUD 29/hr). **Delete the instance after pulling
results; do not leave it stopped.**

---

# Outcome: how the pre-committed branches actually fared

Appended 20 Sep 2026 after the probe reported and the panel launched. The point of writing the
branches before the number existed was to be able to check them afterwards, so here is the
check — including where the rule itself was wrong.

## The branch fired correctly

Measured 104.76 s/ligand, 21.83 GPU-hours. With the pre-specified 6% probe-ligand correction
and per-instance setup, `T = 23.8` instance-hours → **branch B**, run the full panel. No
judgement was exercised after seeing the number; the rule decided.

## Three assumptions in the rule were wrong

1. **"Four shards run in parallel, wall clock ≈ H/4 + S."** The rule never checked GPU quota.
   L4 quota is 1 per region, so parallelism requires four *regions* — which exists, but was
   luck rather than design. Worse, L4 capacity itself turned out to be the binding constraint,
   so four-way parallelism has not been achieved at all.
2. **"Per-instance setup `S ≈ 0.171 h`."** Measured across instances it ranges 848–1473 s for
   batch A alone, because model-weight download time varies by roughly 2x. `S` is not a
   constant.
3. **The rule priced no preemption at all**, which is flagged in the original text as a known
   omission. It turned out to be the term that dominated, and it is what forced the switch to
   on-demand (`analysis/boltz2/SPOT_VIABILITY.md`).

## What the correction was worth

The 6% probe-ligand correction was measured (probe ligands 22.9 heavy atoms vs the panel's
24.3) and applied conservatively. It made no difference to which branch fired — `T` would have
been 22.5 rather than 23.8, both inside branch B. Cheap insurance that happened not to be
needed; worth keeping, because it could not be known in advance which way it would fall.

## The honest summary

The gate did its job: it removed the ad-hoc spending judgement the user objected to, and the
decision was made by a rule fixed in advance. But **a rule written in advance can still encode
unmeasured assumptions**, and three of this one's did. A pre-registration is protection against
motivated reasoning, not against being wrong — and the failure mode it does not catch is an
assumption nobody thought to question, such as whether four GPUs could be obtained at all.
