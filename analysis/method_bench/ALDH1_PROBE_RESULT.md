# ALDH1 Boltz-2 probe on L4: the measured number, and what it cost me to doubt it

20 Sep 2026. Measured on the hardware the work would actually run on, per
`ALDH1_COST_DECISION.md`, whose branches were committed before this number existed.

## The number

| | s/ligand | projected GPU-h for 750 |
|---|---|---|
| Mpro L4 probe (`boltz2_probe_L4.json`, 3 Sep) | 100.99 | 21.07 |
| **ALDH1 L4, measured** | **104.76** | **21.83** |

Batch A 4/4 in 848.5 s (MSA pre-staged), batch B 8/8 in 838.1 s, on an NVIDIA L4 23 GB,
g2-standard-8 spot, us-central1-b. `boltz2_probe_L4_aldh1.json`.

## What I got wrong

Twice in the run-up I argued that the Mpro figure could not transfer, because ALDH1 is modelled
as a 494 aa dimer A+B -- about 988 residues against Mpro's ~612 -- and Boltz-2 scales
superlinearly in sequence length. The two numbers are **3.7% apart**. Cost here is dominated by
the per-ligand diffusion and affinity work, not by protein length, so the extrapolation I
distrusted was sound.

`THIRD_TARGET.md`'s "~21 GPU hours" was therefore right all along. The Kaggle abort recorded in
`ALDH1_BLOCKED.md` was caused *purely* by the T4-for-L4 hardware substitution (421.5 vs 104.8
s/ligand, 4.0x), not by anything about the target.

The lesson survives inverted. "Measured on the wrong hardware" was a real defect and the abort
that caught it was the system working. But "this estimate cannot transfer" was itself an
unmeasured claim, and I asserted it as confidently as the original error. A prediction about
cost is not safer for being pessimistic.

## Two other measurements worth keeping

**MSA generation is ~14% of batch A, not most of it.** Batch A took 984 s building the MSA from
the server and 848 s with it pre-staged -- a 136 s difference. The pre-staging change
(354ecebd) is still worth having, because it is paid per instance and per preemption, but the
commit message justified it with "most of it on that MSA", which the next run disproved.

**The two-point marginal went negative again** (-2.6 s/ligand, batch B faster than batch A).
The `max(_proj_diff, _proj_avg)` guard did its job and fell back to the scaled-chunk
projection. Without it this would have been a false PASS.

## Spot reliability, measured at n=2

| attempt | zone | outcome |
|---|---|---|
| 1 | us-central1-a | **preempted at 30 min**, mid batch B |
| 2 | us-central1-b | completed in 31 min |

Between them, every create across all three us-central1 zones returned a stockout for roughly
nine hours. Two data points establish no rate, but they do establish that a 30-minute window is
not guaranteed, and that capacity can be absent region-wide for hours.

## Decision under the committed rule

`ALDH1_COST_DECISION.md` fixed the branches in advance:

- `H = 1.06 x 21.826 = 23.14` GPU-hours (the 6% is the probe-ligand size correction).
- Per-instance setup `S ~ 0.171 h`: ~3.1 min boot and install, plus ~429 s of model download
  inside batch A once its 4 ligands are charged at the marginal rate.
- Four shards: `T = 23.14 + 4 x 0.171 = 23.8` instance-hours.

`12 < T <= 30` is **branch B: run the full 750, sharded 4x**. At the deliberately high $0.35/hr
planning rate that is **~$8.3**; at the ~$0.22/hr g2-standard-8 spot rate, ~$5.2.

Sharding needs four concurrent L4 spots. The regional quota is 1, but it is 1 in each of
us-central1, us-east1, us-east4, us-west1, us-west4 and europe-west4, so four shards run across
four regions. Cross-region reads of a US-CENTRAL1 bucket cost about $0.02/GB on ~100 MB, which
is immaterial.

## What branch B did not price

Preemption. The committed rule charges `S` once per instance; attempt 1 shows it is really once
per *preemption*. Read from the code rather than estimated, a restart costs:

| | |
|---|---|
| boot + install | ~3 min |
| model weight download (inside batch A) | ~7 min |
| **re-running the 12 probe ligands** | **~21 min** |
| in-flight chunk lost (25 x 105 s) | up to 44 min |

so **~31 min fixed plus up to 44 min of lost work**, not the "~10 min" an earlier draft of this
file claimed. The 21 min is avoidable waste: `ckpt_pull` runs *after* batches A and B, and
`rest` starts at `names[12:]`, so a resumed shard re-scores 12 compounds it already holds.
Left in place deliberately — restructuring the probe path on a worker about to run unattended
for ~22 hours is a worse risk than the minutes it saves, and if preemptions turn out frequent
the fix can be made against data rather than against a guess.

The chunk loss is the dominant term and *was* worth acting on, because 44 min of exposure
against an observed ~30 min time-to-preemption means most chunks would never commit. `CH` is
now env-overridable (`BOLTZ_CHUNK`), default unchanged at 25 so nothing already measured moves.

Any full-panel supervisor needs a cumulative instance-hour cap, because the failure mode is
churn that spends without finishing.
