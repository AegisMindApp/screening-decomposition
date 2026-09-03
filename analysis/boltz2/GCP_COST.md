# What Boltz-2 would cost on Google Cloud: about $5–15 per target

**3 September 2026.** Anchored on our own measurement (229.4 s/ligand marginal, Kaggle P100,
Mpro dimer 2×298 aa, bf16 AMP, torch 2.4.1+cu121) rather than on the paper's figure.

## Estimate for one target (751 compounds)

| card | est. s/ligand | GPU-hours | on-demand | spot | basis |
|---|---|---|---|---|---|
| **A100 40GB** | 20 | 4.2 | **$15** | **$5** | native bf16 + tensor cores |
| L4 | 50 | 10.4 | $7 | $3 | 121 TFLOPS bf16 vs A100's 312 |
| T4 | 100 | 20.8 | $11 | $4 | fp16 tensor cores, **no native bf16** |
| P100 (measured) | 229 | 47.9 | $26 | $10 | Pascal: neither tensor cores nor bf16 |

GCP us-central1 list rates: A100 40GB $3.67/hr on-demand, ~$1.10 spot; L4 ~$0.70; T4 ~$0.55
all-in. Storage and egress here are under $1 and ignored.

All three targets (Mpro 751 + Factor Xa 854 + PD-L1 417 = 2,022 compounds) on A100 spot:
**about $12.**

## The one number that is measured, and the ones that are not

**Measured:** 229.4 s/ligand on a P100. Everything else in that table is an *extrapolation*
from relative bf16 throughput.

There is one corroboration worth noting: dividing our P100 rate by 11.5 gives 20 s/ligand,
which is exactly what the Boltz-2 paper reports. That is consistent with the paper having run
on A100-class hardware, and it makes the A100 row the best-supported estimate rather than the
most optimistic one. It is still not a measurement.

**Before committing to a full run, re-run the same two-batch probe on the actual instance.** It
costs ~10 minutes and replaces the extrapolation with a number. That is precisely the step that
prevented this from consuming a whole week of Kaggle quota.

## Practical notes

- **Spot is safe here.** The worker already writes `boltz2_scores.json` after every 25-compound
  chunk, so a preemption loses at most one chunk, not the run.
- **The dimer is the expensive choice.** A monomer would be roughly 2–4× cheaper and
  biologically wrong for Mpro, whose S1 subsite is formed by the partner protomer. Keep the
  dimer; the saving is not worth the wrong construct.
- **Setup is not free in time**, though it is in money: a GCE instance needs CUDA drivers and
  the boltz install. The Kaggle worker ports over, minus the Pascal bootstrap — on sm_70+ the
  stock torch works and the `get_arch_list()` check should pass untouched.

## Before creating anything

Billing is enabled on `aegismind-tpu` and **the Billing Budget API is disabled**, so no budget
alert can exist. Enable `billingbudgets.googleapis.com` and set a budget first. Verified the
same day that nothing is currently running — no TPU VMs across 15 zones, no GCE instances in
any zone — so there is no spend accruing today, and the point is to keep it that way.

## Is it worth it

Cost was never the barrier: **$5 is not a budget decision.** The barrier was that Kaggle's free
tier hands out a P100, which cannot run this model at a usable rate. Whether to spend the $5
depends on whether the underlying question still matters — adjudicating the Boltzina /
reliability-paper disagreement on our benchmark — not on the price.
