# DiffDock's top-ranked pose lands outside the intended binding site about two thirds of the time

**Run 31 August 2026**, all three shards complete. `arm2_blind_scores.json`,
`arm2_blind_placement.json`.

This was not the experiment we set out to run. Arm 2 was meant to isolate pose **sampling** —
does a learned pose generator, scored identically, change the ranking? A routine distribution
check before computing any AUROC showed it could not answer that, and the reason is a result in
its own right.

## Measured on all 746 scored compounds

| | arm 1 (Vina poses) | arm 2 (DiffDock, blind) |
|---|---|---|
| poses scoring **positive** (steric clash) | **0%** | **62.9%** |
| mean | −7.58 | **+6.92** |
| median | −7.60 | +1.91 |
| worst | −3.71 | **+158.59** |

A Vina affinity of +158 kcal/mol is not a weak binder. Arm 1's scores are entirely favourable and
span a plausible 6.6 kcal/mol range; two thirds of arm 2's are not physical at all.

## Where the poses actually went (90 measured)

| | |
|---|---|
| within the 11 Å box half-width | **35.6%** |
| **outside the docking box** | **64.4%** |
| median distance from box centre | **14.1 Å** |
| maximum | **26.0 Å** |

The 62.9% clash rate and the 64.4% out-of-box rate agree closely, which is the expected
relationship: a ligand placed outside the pocket has nowhere good to sit.

## This is DiffDock working correctly

DiffDock performs **blind docking** — it searches the entire protein rather than a nominated
site. 7VU6 is a homodimer, so a blind search has the second copy of the active site and the dimer
interface available, and it takes them. Nothing here is a malfunction; it is a method answering
the question it was built for, on a benchmark that asked a different one.

Our benchmark is defined by a 22 Å box on one site. A method that does not aim at that site
cannot be compared, on this benchmark, to one that does.

## What it means for using DiffDock in screening

Published benchmarks put DiffDock's own confidence score at **AUROC 0.538–0.56** for virtual
screening — near chance. This adds the mechanism underneath that number: on a pocket-defined
target, its top-ranked pose is more often somewhere else on the protein than in the site. A
confidence score cannot rank binding at a site the pose is not in.

The practical consequence is that DiffDock is not a drop-in replacement for site-directed docking
in a retrospective screen, and any pipeline using it that way needs a placement check — the exact
check that stopped us reporting an AUROC here.

## What we did *not* conclude

We did **not** compute an arm-2 AUROC and report it. Filtering to in-box poses would have left
about 36% of compounds, far below the 90% coverage floor declared in the pre-registration, so the
survivors would not have been the same benchmark. The sampling clause of the hypothesis remains
**open**, and is being re-run with pose generation constrained to a pocket-truncated receptor
while scoring stays against the full protein.

Nothing here says DiffDock predicts poses badly. It says a blind method on a homodimer does not
reliably choose the site a retrospective benchmark is about, which is a statement about
experimental design, not about pose accuracy.
