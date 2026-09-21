# ALDH1 Boltz-2 is blocked on hardware, not on design

19 Sep 2026. The panel, sequence, pocket and bundle are all built and correct. The run cannot
proceed on Kaggle.

## What happened

Shard 1 aborted itself under the worker's pre-registered ceiling:

```
== accelerator
Tesla T4, 15360 MiB
MARGINAL 421.5 s/ligand   PROJECTED 22.01 GPU-hours for 188   CEILING 8.0
ABORT per pre-registration. Not spending the quota.
```

**Kaggle gives a Tesla T4. Boltz-2 costs 421.5 s/ligand there.** I budgeted 101 s/ligand from
`boltz2_probe_L4.json` — an **L4** probe from a GCP run — and described it as "a measured probe
on this exact model". It was measured on the right model and the **wrong hardware**, and the
figure is out by 4.2×.

That is the same error as the stage-2 exhaustiveness misestimate twelve hours earlier: a cost
taken from a convenient measurement rather than from the machine the job will actually run on.
"Measured" is not sufficient. It has to be measured *where the work will happen*.

## What it would cost

| panel | GPU-hours on T4 | |
|---|---|---|
| 750 (built, matches Mpro 748 / FXa 885) | **88** | infeasible |
| 200 | 23 | still beyond a weekly quota |
| 120 | 14 | fits two shards under 8 h each |
| 80 | 9.4 | fits comfortably |

## Why the abort is the system working

The ceiling cost ~1.5 h of probe time instead of ~88 GPU-hours. It was set at 8.0 h in this run
precisely because the worker's own default of 6.0 looked too tight against my (wrong) estimate —
so the mechanism caught an error in the number that was used to configure it.

## Options, none of them free

1. **Shrink the panel to ~120 compounds.** Fits Kaggle. But a 120-compound panel against Mpro's
   748 and FXa's 885 is a weak third leg, and the transferability verdict would rest on it.
   The precision cost can be quantified from stage-1 data the same way the exh=32 panel was
   sized — that should be done before choosing, not after.
2. **Paid compute.** The L4/A100 route in `analysis/boltz2/GCP_COST.md`. 88 T4-hours is roughly
   21 L4-hours; real money, and John's call.
3. **Stop at two targets.** The harness already refuses a transferability verdict on two, which
   is honest. Boltz-2 stays `INSUFFICIENT` and is not deployed against objects.

Nothing here invalidates the ALDH1 selection — it remains the cleanest admissible target
(descriptor baseline 0.5698) and the bundle is built and verified. The block is purely cost.
