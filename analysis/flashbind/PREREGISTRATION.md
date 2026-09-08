# Pre-registration: does a third scoring approach clear the docking residual ceiling?

**Written 8 September 2026, before any FlashBind score exists.**

## The question, taken verbatim from the paper

§3.4 of `docs/papers/screening_decomposition_preprint.md`: *"on these two targets the ceiling
appears to be a property of the scoring approach rather than of the task. A third scoring
approach clearing it would be the test that generalises this."*

Six docking scores residualise to **0.494–0.573**. Boltz-2 residualises to **0.6567
[0.6439, 0.6648]**. One point outside a band is not a generalisation.

FlashBind (bioRxiv 10.64898/2025.12.22.695983, MIT, `AIDD-Lab/FlashBind`) is the sharpest
available third case because it is *architecturally between the two families*: a fast docking
model for placement, a learned EGNN head for scoring. It claims early enrichment competitive
with Boltz-2 at ~50× lower cost. Surfaced by `analysis/recency_monitor` on its first production
scan.

## Reading rules — fixed now

Primary quantity: **residual AUROC** on Mpro — FlashBind's score with the seven descriptors
regressed out, out-of-fold gradient-boosted, mean over 12 fold seeds with the full range, exactly
as in `analysis/docking_value/RESIDUAL_BAND.md`.

- **CEILING IS THE SCORING APPROACH (claim generalises)** — residual > **0.573**, the top of the
  docking band, with the 12-seed range clear of it. Two non-docking approaches beat the ceiling
  and it is a property of the approach.
- **ENRICHMENT AND ORTHOGONAL SIGNAL COME APART** — residual within **[0.494, 0.573]** while raw
  AUROC is competitive with Boltz-2's 0.7913. A method can rank actives well and still carry no
  more target-specific signal than docking. This would *weaken* our §3.4 claim and is the more
  interesting outcome.
- **BELOW THE BAND** — residual < 0.494. Reported as such.

Secondary, reported regardless: raw AUROC; marginal value over descriptors alone
(paired bootstrap, 10,000 resamples) against Boltz-2's +0.093 [+0.062, +0.125]; and the same on
Factor Xa if Mpro completes.

## Controls — the arm is void without these

1. **Positive control.** The descriptor baseline on the same compounds must reproduce **0.7654**
   to ±0.005. If the panel is not the panel, no comparison is readable.
2. **Coverage floor.** ≥ 90% of 751 compounds scored, with failures reported by reason.
3. **Permutation.** Shuffled labels must give AUROC in [0.45, 0.55].

## Leakage — checked BEFORE the residual is read

FlashBind's binary task ships against MF-PCBA (PubChem-derived). Our panel is ChEMBL-derived
SARS-CoV-2 Mpro. Mpro is a heavily used benchmark target and these sets may intersect.

Before any number is interpreted, establish whether Mpro, or our specific compounds, appear in
FlashBind's training data. If the training set cannot be established, **the result is reported
with leakage unresolved and is not used to support the generalisation claim** — a model trained
on the target it is being tested on would clear the ceiling for the wrong reason.

## Abort conditions

- ESM3 representations cannot be generated (gated weights) → report blocked, do not substitute a
  different encoder and call it FlashBind.
- Checkpoints unavailable or licence-incompatible → report blocked.

## Pre-committed disclosure

Reported either way, including the outcome that weakens §3.4. The manuscript is not yet submitted
to a journal; if this lands before submission it goes in, and if it contradicts the current text
the text changes.

## Pocket agreement — reported regardless of the residual

**Added 8 September 2026, still before any FlashBind score exists.**

The six docking scores that define the band [0.494, 0.573] were all produced inside the same
22 Å box centred on **[9.05, 8.90, −1.51]** (`bundle_CHEMBL4523582_7VU6/shard*/manifest.json`,
receptor 7VU6). Boltz-2 — the one method that cleared the ceiling at 0.6567 — chose its own
binding mode. FABind+ also predicts its own pocket.

So if FlashBind clears the ceiling, the two methods that cleared it are also the two methods
that were not box-constrained, and "the ceiling is a property of the scoring approach" has a
live alternative: **the ceiling is a property of being box-constrained**. Same shape as the
blind-vs-pocket-conditioned error that voided the DiffDock arm.

`pocket_indices.lmdb` is an output of the FABind+ stage, so the check is free. Report:

- the centroid of FABind+'s selected pocket residues per compound, in receptor coordinates
- its distance to [9.05, 8.90, −1.51], as a distribution over compounds
- the fraction of compounds whose pocket centroid falls inside the 22 Å box

Reading:

- **Pockets agree** (median centroid distance < 11 Å, i.e. inside the box) → FlashBind is a
  clean third *scoring* approach and the reading rules above stand unmodified.
- **Pockets disagree** → the residual is reported with the pocket confound named, and it does
  not by itself support the generalisation claim. The finding is then about the box, not the
  score, and that is the more useful result.

This is a **reported quantity, not a gate**. It cannot void the arm; it determines which of two
claims the number is allowed to support.

## Second leakage surface: FABind+, not just the scorer

`LEAKAGE_CHECK.md` cleared the *affinity head* — 0 protease assays among MF-PCBA's 105 AIDs.
That says nothing about the pose generator. FABind+ ships a PDBbind-trained checkpoint, and
PDBbind contains SARS-CoV-2 Mpro complexes. This advantages **pose placement** on our target,
in the direction of better performance.

Disclosed here, not treated as a blocker: the arm asks whether a third scoring approach carries
orthogonal signal, and a well-placed pose is the input to that question rather than the answer.
It is stated in the write-up either way.

## Abort decision, taken before the result exists

**8 September 2026.** FABind+ requires the real `torchdrug` (`inference_mol_utils.py:81`),
which declares `python<3.11` against Kaggle's 3.12 image. v5 installs past the metadata cap and
makes FABind+'s exact call as the test.

**If torchdrug fails to run on 3.12, the arm is reported BLOCKED.** No Python 3.10 environment,
no substituted featuriser. Recorded here before the v5 result is known, so it is a decision
rather than a rationalisation of one.

Two reasons it is the right call. A 3.10 environment would break the Pascal torch pin the whole
run rests on, so it is not a small change. And the affinity module's 56-dim featuriser is the
obvious substitute and the wrong one — it emits nodes only, no `edge_list` or `edge_weight`, and
swapping an encoder while still calling the output FlashBind is what the abort conditions above
already forbid.

A blocked arm leaves §3.4's question open with the reason stated, which is an honest outcome.
A patched-together arm answers it wrongly.

**Condition tested, 8 September 2026 — not met.** torchdrug imports fine on Python 3.12. It
failed on `rdkit.Chem.Draw.mplCanvas`, removed from current rdkit but present in every cp312
build up to **2024.3.6**. torchdrug is unmaintained (v0.2.1, July 2023; that import untouched
since 2021), so the fix belongs on the rdkit side: a supported version pin, not a patch to
either library and not a substituted component. The arm continues. The abort rule above stands
unchanged for a failure that actually meets it.
