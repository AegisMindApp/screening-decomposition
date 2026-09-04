# Pre-registration — Boltz-2 on Factor Xa and PD-L1

**Written 4 September 2026, before either target runs.** Mpro arm 1 (blind) is complete:
Boltz-2 0.7913, descriptors 0.7654, combination 0.8588 (+0.0934, CI excludes zero). Mpro arm 2
(pocket-conditioned) is running.

## Configuration

Pocket-conditioned from the outset — Mpro arm 1 established that blind is a handicap against a
pocket-constrained comparator, so there is no reason to repeat it. `diffusion_samples = 1`,
`potentials` off, `affinity_mw_correction` off, matching the Mpro arms so results are
comparable across targets.

Pockets derived from each receptor at the Vina box centre, 8.0 Å cut-off:

- **Factor Xa** — 886 compounds, 405 actives (46%), chain A, 230 tokens. 13 residues, led by
  Gly206/Trp205/Glu207/Ser204 — the S1 pocket. Box centre [7.96, 5.32, 21.77].
- **PD-L1** — 455 compounds, 359 actives (79%), chains A+B, 246 tokens. 14 residues spanning
  **both** chains (B105, B106, B104, A39, A98…), correct for PD-L1 inhibitors, which bind at
  the dimer interface. Box centre [4.95, 37.59, 21.94].

## Reading rules — per target, fixed now

Primary endpoint `affinity_probability_binary`, out-of-fold, same folds across arms.

| target | descriptor baseline | Boltz-2 must reach |
|---|---|---|
| Factor Xa | **0.7022** | **> 0.7422** with CI excluding zero |
| PD-L1 | **0.8740** | **> 0.9140** with CI excluding zero |

Plus, per target: combination (descriptors + Boltz-2) must beat descriptors alone by > +0.04
with CI excluding zero. Controls unchanged — positive control reproducing the baseline above
± 0.005, permutation inside [0.45, 0.55], dropout bounded rather than the degenerate rule
argued.

## PD-L1 is a compromised benchmark — recorded before it runs

**Its descriptor baseline is 0.8740.** Seven free physicochemical properties separate this
panel almost perfectly. That is not a sign of a good benchmark; it is the confound we already
documented — the panel is **79% active** and its inactives are **127 Da heavier** than its
actives. Vina scores 0.5661 on it with a confidence interval that includes chance.

So a Boltz-2 result here mostly tests **whether it can out-predict a molecular-weight
artefact**, not whether it predicts binding. Both outcomes are weak:

- Below 0.914 → uninformative; it fails against a baseline that is itself measuring size.
- Above 0.914 → impressive but hard to attribute, since the target signal and the size
  confound point the same way.

It is run for completeness and because a null on a confounded panel is still worth recording,
but **no strong conclusion should be drawn from PD-L1 in either direction**, and that is stated
here rather than after the number lands.

Factor Xa carries the real weight: 46% actives, a 0.7022 baseline, and the one target that
passed every admissibility gate in the docking work.

## Compute and budget

GPU quota is 1, so these run serially after Mpro arm 2. Estimated at Mpro's 101 s/ligand:
FXa 24.9 GPU-h, PD-L1 12.8 GPU-h. Both proteins are far smaller than Mpro's 596-token dimer
(230 and 246 tokens), and Boltz-2's pairwise stack scales superlinearly with tokens, so the
true cost should be materially lower — the probe will measure it.

**Ceilings set to keep cumulative spend inside the AUD 100 project budget**: Factor Xa
**18 GPU-hours**, PD-L1 **10 GPU-hours**. Spend to date is ~AUD 27 (Mpro arm 1) with ~AUD 27
running (arm 2), leaving ~AUD 46; these ceilings cap the worst case at ~AUD 36. If a probe
projects above its ceiling the run abandons itself, as it did on the Kaggle P100.
