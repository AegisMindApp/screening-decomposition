# Every receptor we have ever docked against has zero hydrogen-bond donors

**Found 31 August 2026**, while checking whether a generic docking-methodology checklist had
anything we had missed. Its second item — *"protonation states… can flip docking poses
entirely"* — pointed at something real.

## What the receptors actually contain

| receptor | HD (polar H) | N (donor) | NA (acceptor) | OA |
|---|---|---|---|---|
| Mpro 7VU6 | **0** | **0** | 754 | 852 |
| Factor Xa (CHEMBL244) | **0** | **0** | 313 | 346 |
| PD-L1 (CHEMBL3580522) | **0** | **0** | 315 | 360 |
| CHEMBL612545 | **0** | **0** | 644 | 722 |

Not one polar hydrogen, and not one donor nitrogen, in any receptor in the programme.

## The cause, in one line of source

`analysis/docking_gates/redock_control.py`:

```python
AD4_TYPES = {"C": "C", "N": "NA", "O": "OA", "S": "SA", ...}
```

**Every nitrogen is mapped to `NA` — hydrogen-bond acceptor — unconditionally**, and no
hydrogens are added at any point. AutoDock Vina identifies a donor by the polar hydrogen (`HD`)
attached to it; with no `HD` atoms, the donor set is empty by construction.

The function is honest about being minimal — *"Deliberately the same convention as production
rather than a better one: the gate should test the receptor the pipeline actually docks
against."* The gate was faithful to production. Production was wrong.

## What it costs, quantified on Mpro

All **754** nitrogens are donor-capable and all are typed as acceptors:

| | count |
|---|---|
| backbone amide N — a donor, typed as acceptor | **597** |
| side-chain donors (NZ, NE, NH1/2, ND1/2, NE1/2) | **157** |
| nitrogens correctly available as donors | **0** |

So the scoring function could never form a **protein-donor → ligand-acceptor** hydrogen bond —
roughly half of hydrogen-bonding physics — and simultaneously counted 597 backbone amide NH
groups as *acceptors*, which is not merely missing but inverted. Waters were dropped too
(`HOH/WAT/DOD` filtered), so conserved bridging waters are gone as well.

## Why this is a candidate explanation, and why it is not yet an answer

It is mechanistic, systemic, and of the right sign: a scoring function blind to half its
hydrogen bonds — and actively wrong about the other half — should rank poorly. It fits what arm
3 already measured, that **swapping only the scoring function moves Mpro 0.399 → 0.539 while the
atom coordinates stay identical**. If the geometry is adequate and the scorer is the problem, a
scorer given inverted H-bond typing is a concrete reason why.

But it does **not** explain the pattern on its own. All four receptors were prepared the same
way, yet Factor Xa reaches **0.6657** and PD-L1 **0.5948** while Mpro sits at **0.427**. A
defect present in every run cannot by itself account for results that differ across runs. The
honest reading is that this is a real defect of unknown effect size, not a diagnosis.

## The experiment

Re-prepare the Mpro receptor with polar hydrogens added at pH 7.4 and correct donor/acceptor
typing, change nothing else, and re-dock the same 753 compounds. Everything needed already
exists — the compound set, the box, the labels, the harness, and a 0.408 baseline at
exhaustiveness 4 to compare against.

Three outcomes, all informative:

| result | reading |
|---|---|
| AUROC rises materially | preparation was a substantial part of the failure, and **every docking number this project has published is affected** |
| AUROC unchanged | Vina's H-bond term is not what is wrong here, and the defect is real but not the cause |
| AUROC falls | the below-chance result is not a preparation artefact at all |

Until that runs, no docking result here should be strengthened *or* dismissed on the strength of
this finding.

## What else the checklist surfaced

- **Consensus rescoring** — "use one scoring function to generate poses, another to re-score" is
  already confirmed on our own data by arm 3 (+0.140 AUROC).
- **Y-randomization** — scramble the labels and confirm the pipeline does not "work" on noise.
  Cheap, never run here.
- **Conserved waters** — dropped by the same converter; untested.
- Its recommendation to use DiffDock for screening is contradicted by 2025–26 benchmarks putting
  DiffDock's confidence score at AUROC 0.538–0.56 for that purpose.
- It omits the two controls that have been most informative for us: a **property/descriptor
  baseline** (0.763 on Mpro, which no docking configuration has beaten) and **pre-registering
  the reading rule** before the numbers exist.
