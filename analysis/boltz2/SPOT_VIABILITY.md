# When spot GPU is cheaper on paper and more expensive in practice

20 Sep 2026, measured during the ALDH1 Boltz-2 panel run. Written because the conclusion is
reusable and the arithmetic is not obvious in advance.

## The rule

Spot is viable only when

    time-to-preemption  >>  un-checkpointed prologue

Below that, instances are reclaimed before committing anything and the run churns forever,
spending money at full rate while making no progress. This is not a slowdown; it is a
**non-convergence**, and it looks identical to "still working" from outside.

## Measured on L4 spot, us-central1 / us-east1 / northamerica-northeast2

| instance | zone | time to preemption |
|---|---|---|
| probe attempt 1 | us-central1-a | 30 min |
| shard 0 first | us-central1-c | 106 min |
| shard 0 second | us-central1-a | **14 min** |
| shard 1 | us-east1-d | **18 min** |

Median ~24 min. Confirmed as genuine `compute.instances.preempted` operations, not
self-termination — worth checking, because a startup script that runs `shutdown -h now` on any
worker exit produces an identical external signature.

Boltz-2's prologue was ~28 min: batch A ~14 min (model download plus 4 ligands) and batch B
~14 min (8 ligands). Neither was checkpointed. So roughly half of all instances committed
nothing, and across 40 minutes with four shards running, **zero compounds were scored**.

## What fixed it, and what only looked like it would

**Worked:**
- `ckpt_pull` before the probe, so a resumed shard skips the ~28 min prologue entirely
  (commit 0d6ba9f9). Only helps once a shard has its first checkpoint.
- Chunk size matched to time-to-preemption. The in-flight chunk is what a preemption destroys.
- Staging the MSA to GCS so batch A does not rebuild it — though this turned out to be only
  ~14% of batch A (984 s vs 848 s), not "most of it" as first claimed.

**Did not help, because it addressed the wrong thing:** most of a night's engineering on
preemption resilience — zone rotation, adoption, duplicate pruning, spend caps. All of it is
correct and worth keeping. None of it addresses a prologue longer than the instance's life.

## The economics, once measured

| | spot | on-demand |
|---|---|---|
| useful work per ~24 min life | ~14 min (10 min is boot + weight download) | continuous |
| instance-hours for 750 compounds | ~38 | 23.8 |
| cost | ~$13 | ~$20 |
| wall clock | 10-14 h plus stockout gaps | ~6 h |
| risk | may stall at the spend cap | none |

**Spot saved about $7 and cost about 8 hours plus a real chance of never finishing.** Switched
to on-demand on John's authorisation.

## The diagnostic error worth remembering

For roughly ninety minutes the supervisors reported "no L4 spot capacity in any zone" and I
reported that upward as a market condition. It was not. `cleanup()` ran inside the zone-search
loop and deleted instances that `gcloud create` had returned non-zero for and then successfully
created moments later. The searcher was eating its own winnings and truthfully reporting that
it had found nothing.

Both of my explanations before the right one were confident and wrong: first "the market has no
capacity", then "gcloud exit codes are unreliable" asserted as fact when it was an untested
theory. The evidence that settled it was one `ps` (duplicate supervisors) and one
`gcloud compute operations list` (inserts with no deletes). Neither took more than a minute.

**The generalisable form:** when a searcher reports finding nothing, check whether it is
destroying what it finds before concluding there is nothing to find. An absence produced by
your own instrument is indistinguishable from an absence in the world — the same failure mode
as `feedback_silent_failure_looks_like_data`, one level up.
