# L4 probe: 101 s/ligand, 21.07 GPU-hours — aborts against a 20.0 ceiling by 5%

**3 September 2026.** Run 4 on a GCE L4 (on-demand, us-central1-b) cleared every install
checkpoint and produced a clean measurement. The pre-registered ceiling then fired.

    batch A   4/4 in 851 s   (model load, weights, MSA generation -- one-time)
    batch B   8/8 in 808 s   (reused MSA -- the marginal rate)

    MARGINAL     101.0 s/ligand
    PROJECTED    21.07 GPU-hours for 751 compounds
    CEILING      20.0  ->  ABANDONED

## Where this sits

| card | s/ligand | GPU-hours (751) | vs L4 |
|---|---|---|---|
| P100 (Kaggle, measured) | 229.4 | 47.9 | 2.3× slower |
| **L4 (measured)** | **101.0** | **21.07** | — |
| paper's figure | 20 | 4.2 | 5.1× faster |

**My L4 estimate was 50 s/ligand. The measurement is 101 — I was 2× optimistic.** That is the
entire reason the probe exists, and it is the second time this session that an extrapolated
throughput figure was wrong in the optimistic direction.

The paper's ~20 s/ligand remains 5× faster than a real L4. That is consistent with it having
been produced on A100-class hardware, which this project cannot obtain — `NVIDIA_A100_GPUS` and
`A2_CPUS` are both 0.

## The awkward part, stated plainly

21.07 against a ceiling of 20.0 is an overshoot of **5%**, or about **AUD 2** at L4 on-demand
(~AUD 27 total versus the ~AUD 25 the ceiling implied). Both figures sit inside the AUD 100
project budget.

**I am not raising the ceiling myself.** It was pre-registered at 20.0 before the L4 probe ran,
and moving a threshold immediately after seeing a number that just misses it is precisely the
practice this project's credibility rests on refusing — even when, as here, the threshold is a
resource limit rather than a scientific one, and even when the sum is trivial. The scientific
reading rules were never in question; the point is that a line moved after the fact is not a
line.

Three ways forward, for the person whose money it is:

1. **Authorise the overage** — full 751 compounds, ~21.1 GPU-hours, ~AUD 27. Cleanest
   scientifically: the benchmark stays exactly the panel every other result on this page uses.
2. **Random stratified 700 of 751** — 19.6 GPU-hours, inside the existing ceiling with no line
   moved. Power is barely affected (the interval widens ~4%), and the subset is drawn at random
   rather than by any property, so it avoids the selection bias that damaged arm 2 and PD-L1.
3. **Stop.** The measurement is itself a result: Boltz-2 costs 101 s/ligand on the best card
   available to us, and the hypothesis stays untested.

## Still untested

No Boltz-2 affinity value has been compared against a label. The Boltzina / reliability-paper
disagreement remains exactly as open as before.

## Spend so far

Roughly 1.2 hours of `g2-standard-8` across four boots, mixed spot and on-demand — on the order
of **AUD 1**. Every failed run uploaded its log and terminated itself; nothing idled.
