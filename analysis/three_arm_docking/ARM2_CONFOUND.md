# Arm 2 as run measures the wrong thing: DiffDock blind-docks away from the pocket

**Found 31 August 2026, before any arm-2 AUROC was computed.** Caught by a distribution check
that should be routine and nearly was not run.

## The signal

| | arm 1 (Vina poses) | arm 2 (DiffDock poses) |
|---|---|---|
| positive (clashing) Vina scores | **0%** | **63%** |
| mean | −7.58 | **+6.88** |
| max | −3.71 | **+67.43** |

A Vina affinity of +67 kcal/mol is a severe steric clash, not a weak binder. Arm 1's scores span
−10.4 to −3.7 and are all favourable; arm 2's are mostly not physical.

## The cause

DiffDock performs **blind docking** — it searches the whole protein rather than a specified site.
Our benchmark is defined by a 22 Å box centred on [9.05, 8.9, −1.51], and 7 of 12 sampled rank-1
poses land outside it, at a mean 12.7 Å from the centre and as far as 24.8 Å:

```
CHEMBL1188627   1.5 A   (in the pocket)
CHEMBL1874417   2.3 A
CHEMBL121556    2.5 A
...
CHEMBL1431798  17.4 A   (elsewhere on the protein)
CHEMBL1770312  24.3 A
CHEMBL1356238  24.8 A
```

7VU6 is a homodimer, so a blind search has a second copy of the site plus the dimer interface to
find. Nothing here is DiffDock malfunctioning — it is doing what it was built to do, on a
question our benchmark did not ask.

## Why this cannot be analysed as arm 2

The pre-registered contrast was **pose sampling**: does a learned pose generator, scored the same
way, change the ranking? What this run actually measures is **whether a blind search rediscovers
the pocket we defined** — a different question, and one whose answer contaminates the first. A
ligand placed 24 Å away scores badly because it is in the wrong place, not because sampling is or
is not the bottleneck.

Filtering to in-box poses would leave roughly 37% of compounds, far below the **90% coverage
floor** declared in the pre-registration, so the surviving subset would not be the same benchmark.

## What it does establish, and it is worth keeping

On a pocket-defined retrospective benchmark, **DiffDock-L's top-ranked pose is outside the
intended site more often than in it**. That is consistent with published results putting its
confidence score at AUROC 0.538–0.56 for virtual screening: a method that ranks its own blind
poses cannot be assumed to have found the site the benchmark is about.

Shards 0 and 2 are being allowed to finish so this can be quantified across all 751 compounds
rather than the 12 sampled here.

## The fix for arm 2 proper

Constrain pose generation to the pocket by giving DiffDock a receptor truncated to residues near
the box centre, while **scoring against the full protein** as arm 1 does. That keeps the
comparison single-variable — pose source — and removes the site-finding question from it.

Until that runs, arm 2 is **not reportable**, and the hypothesis's sampling clause stays open.
