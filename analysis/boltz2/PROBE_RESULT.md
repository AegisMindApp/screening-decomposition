# Boltz-2 probe: 229 s/ligand on a P100, 47.9 GPU-hours projected — abandoned per pre-registration

**3 September 2026.** The pre-registered compute abort fired. Everything in the pipeline works;
the run is simply unaffordable on the hardware we can get.

## The measurement

    batch A   4/4 ligands in 1251 s   (includes model load, weight download, MSA generation)
    batch B   8/8 ligands in 1835 s   (MSA reused; this is the marginal rate)

    MARGINAL     229.4 s/ligand
    PROJECTED    47.86 GPU-hours for 751 compounds
    CEILING       6.0 GPU-hours       -> ABANDONED

Splitting the probe into two sub-batches was what made this trustworthy. Charging batch A's
one-time costs to the rate would have given ~313 s/ligand; taking batch B alone isolates the
marginal cost, and it still breaches the ceiling by 8×. The abort is not an artefact of
start-up overhead.

## Against the paper's figure

The Boltz-2 paper reports **~20 s/ligand**, which is where our 4.17 GPU-hour estimate came
from. Measured on an accessible free-tier card: **229.4 s/ligand, 11.5× slower**. Two causes,
not separated here: the P100 is several generations old, and we submit the **dimer**
(2 × 298 aa) because Mpro is only active as one and that is what we docked against — a monomer
would be cheaper and biologically wrong.

The paper's throughput claim is explicitly conditioned on "modern parallel high-performance
computing infrastructure". That is consistent with what we measured; it is not a contradiction
of their number. It does mean the method is out of reach at this scale without paid compute.

## Why we are not running a subset instead

The ceiling permits **94 of 751 compounds**. That does not rescue the experiment:

- Our measurement floor is ~0.04 AUROC at n = 751. At n ≈ 94 the interval is roughly 2.8×
  wider, so the floor rises to ~0.11 — while the pre-registered primary rule requires resolving
  a **+0.04** difference against descriptors. The subset provably cannot answer the question.
- Subsetting is the exact failure mode that has bitten this project twice already: arm 2's 149
  size-selected dropouts and PD-L1's 38 missing compounds at −14.3 pp. A random stratified
  sample would avoid *that* bias, but not the resolution problem above.

Running 94 compounds would produce a number we could not interpret, at the cost of the entire
week's quota. That is worse than no number.

## Status of the hypothesis

**Untested.** No Boltz-2 affinity value was ever compared against a label. The disagreement
between Boltzina (arXiv 2508.17555) and the reliability evaluation (arXiv 2603.05532) is
exactly as open as before we started. Nothing here is evidence about Boltz-2's accuracy in
either direction — only about its cost on a P100.

## What the run did establish

The whole path is now working and reproducible: Pascal detection, `torch==2.4.1` +
`torchvision==0.19.1` (cu121) pinned together with `sm_60` kernels verified via
`get_arch_list()`, boltz installed and its CLI validated, CCD data and both weight sets
downloaded, MSA generated once and correctly reused (184 kB `msa/*_0.csv`, NUL-validated), and
affinity JSON harvested for 12/12 ligands across two batches. Only the cost blocks it.

To actually run this, the options are paid GPU time (~48 GPU-hours on a P100-class card, far
less on an A100), or a target with a smaller protein where the per-ligand cost falls.
