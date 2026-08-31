# Constraining DiffDock to the pocket fixes *where* the ligand goes and not *what shape it is in*

**31 August 2026.** Pocket-constrained shards 0 and 1 complete (501 compounds); shard 2 running.

## The constraint worked, completely, for placement

| | blind | pocket-constrained |
|---|---|---|
| pose centroid within the 11 Å box | 35.6% | **100%** (60/60) |
| median distance from box centre | 14.1 Å | **3.6 Å** |
| maximum distance | 26.0 Å | 10.4 Å |

Truncating the receptor to 121 residues around the site removed the homodimer's second pocket and
the dimer interface, and DiffDock then put every measured pose inside the box. As a fix for the
blind-docking confound this is unambiguous.

## And the poses are still unscoreable four times in ten

| | arm 1 (Vina) | blind DiffDock | pocket DiffDock |
|---|---|---|---|
| poses scoring **positive** (steric clash) | **0%** | 62.9% | **38.9%** |
| mean | −7.58 | +6.92 | +3.72 |
| median | −7.60 | +1.91 | **−1.06** |
| worst | −3.71 | +158.6 | **+304.5** |

Every one of these poses is now in the right pocket, and 39% of them still score as steric
clashes — one at +304 kcal/mol. So the residual is **not** a placement problem.

## Two failure modes, now separated

This is the useful outcome of the arm, and it took both runs to see:

1. **Site selection.** DiffDock blind-docks; on a homodimer it takes the other site or the
   interface about two thirds of the time. **Fixed** by constraining the receptor.
2. **Local geometry.** A diffusion model generates coordinates directly; it is not minimising a
   force field, so bond lengths, angles and contacts need not satisfy one. Vina's scoring is a
   physics-based function that penalises exactly those violations. **Not fixed**, and not
   fixable by placement.

   > **Corrected 31 Aug by PoseBusters (`POSEBUSTERS_RESULT.md`).** Only *contacts* is right.
   > Bond lengths, bond angles, ring flatness and internal energy all pass at **100%** — the
   > internal chemistry is essentially perfect. The sole failure is **minimum distance to
   > protein, 47% pass**: the ligand is placed too close. This sentence overstated the case
   > against generative pose geometry.

The second is a known property of generative docking rather than a defect we introduced, and it
is why published pipelines that use DiffDock poses productively **rescore after a local
minimisation** rather than scoring the raw output.

## Arm 2 is still not reportable, and now for a stated reason

Filtering to non-clashing poses would leave ~61% of compounds — better than the blind run's 36%,
still below the **90% coverage floor** in the pre-registration. Reporting an AUROC over the
61% that happen to score favourably would select on the outcome variable.

**The protocol amendment that would settle it:** rescore with local minimisation
(`gnina --minimize`) so each DiffDock pose is relaxed inside the pocket before scoring. This makes
the comparison *fairer*, not looser: Vina's own poses in arm 1 are already locally optimal under
Vina's function by construction, so minimising DiffDock's poses removes an asymmetry rather than
adding a licence. Recorded here as an amendment before it runs.

## What this already establishes independent of arm 2

**A learned pose generator cannot be dropped into a physics-scored screening pipeline unchanged.**
It needs its site specified and its poses relaxed. Neither requirement is exotic, but a pipeline
that omits them produces scores that are 39–63% not physical — and an AUROC computed over those
would have looked like a finding about sampling.
