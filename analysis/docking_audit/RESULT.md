# Every docking panel this project has run fails at least one bias control

**6 September 2026.** `sweep.py`, results in `SWEEP.json`. Rules fixed in `PREREGISTRATION.md`
(`e5c305113`) before any unaudited panel was scored; Amendment 1 corrects a control copied from
the wrong panel. Gate **PASSES** on all four positive controls.

This is not an attempt to rescue anything. It applies three diagnostics developed this week —
descriptor baseline, marginal value, residual AUROC — to every panel on disk, to establish which
historical docking claims were readable at all.

> **`CHEMBL612545` is not a target.** Added 21 September 2026, after the identity was
> established. The receptor used with this panel is genuine PD-L1 (PDB 5J89), but the ligand set
> is not: `CHEMBL612545` is a ChEMBL **`UNCHECKED`** record — `target_type` UNCHECKED, no
> organism, **zero target components**, 2,317,536 heterogeneous activities. It is a catch-all
> bucket for bioactivity with no validated target assignment. Human PD-L1 is `CHEMBL3580522`.
> Of twelve sampled "actives", one has any recorded activity against it; the remainder are
> largely HCN1 channel blockers and IL-6 release inhibitors.
>
> **This changes how the row below should be read.** Its "no target-specific signal" is not a
> finding about docking — it is the expected result for a set assembled from unrelated assays,
> which has no shared pharmacophore for any method to recover. The row is retained because
> removing it would hide a real error, but it must not be counted as a docking panel.
>
> Full documentation of the identity error is in the companion benchmark study —
> **https://github.com/AegisMindApp/retrospective-benchmark** (`README.md`,
> `PREREGISTRATION.md`) — where the same identifier is the subject rather than an aside.

## Result: 0 of 3 docking panels survive, plus one non-target control

Three real panels, then `CHEMBL612545`, which is listed separately below because it is
not a target and is not counted in the three.

| target | n | act% | descriptors | docking | marginal value | residual | verdict |
|---|---|---|---|---|---|---|---|
| Mpro | 751 | 34.2 | 0.7673 | 0.4530 | −0.0017 [−0.0038, +0.0003] | 0.5172 | confounded; no marginal value |
| Factor Xa | 854 | 45.4 | 0.7107 | 0.6775 | +0.0172 [−0.0032, +0.0378] | 0.5750 | confounded; no marginal value |
| **PD-L1** | 320 | 78.1 | **0.9145** | 0.5356 | +0.0296 [+0.0021, +0.0638] | 0.5556 | confounded; no marginal value; imbalanced |
| CHEMBL612545 | 414 | 17.9 | 0.6921 | 0.5780 | −0.0103 [−0.0229, +0.0012] | 0.4865 | no marginal value; **no target-specific signal**; imbalanced |

Descriptor baseline: seven free properties, out-of-fold logistic regression. Marginal value:
AUROC(descriptors + docking) − AUROC(descriptors), paired bootstrap, 10,000 resamples. Residual:
docking score with descriptors regressed out, out-of-fold gradient-boosted, mean of 12 fold seeds.

## What each row means

**PD-L1 is the finding.** Seven free descriptors reach **0.9145** on this panel. It is very nearly
separable from ligand properties alone, so no docking conclusion drawn on PD-L1 is evidence about
binding — it is evidence that actives and inactives differ in size and lipophilicity. The panel
was already flagged as confounded in this project's records; this quantifies it. Docking's own
AUROC is 0.5356, barely above chance, and its marginal value (+0.0296) sits below the 0.0393
protocol floor even though its interval excludes zero.

**`CHEMBL612545`'s residual of 0.4865 [0.4571, 0.5026] is an UNINTENDED positive control for the
residual metric** — distinct from the four designed positive controls the gate passes, noted at the
top of this file. It was not built to be one. This is the only row in the table with no target behind it, so the
correct answer for a residual metric is "nothing here" — and it is the only row whose interval
contains 0.50.

Read as a statement about docking the row would be vacuous: there are no target-specific actives
to rank, so failing to rank them is not informative. Read as a statement about the *instrument*
it is the most useful row present, because it demonstrates the residual metric can return chance
when chance is the truth. Every other diagnostic here reports a shortfall; without a case like
this one, a reader cannot tell a metric that detects absence from a metric that only ever
reports absence.

This mirrors the companion study's finding that the debias gate ranked the same non-target as its
most admissible library: a set assembled from unrelated assays has no coherent property structure,
which is invisible to a gate asking only whether actives separate from decoys.

**Mpro and Factor Xa reproduce their published numbers** and are here as controls, not findings.
Both are property-confounded by the 0.70 criterion; neither shows resolvable marginal value.

## Limitations, stated rather than discovered

- **CHEMBL612545 is audited on a subset.** The full panel is 2,249 compounds at 3.3% active; only
  414 carry SMILES on disk, and that subset runs at 17.9% active — roughly five-fold enriched. It
  is **not** the panel the original work used, and the verdict applies only to what was auditable.
- **PD-L1 at 78% active** is severely imbalanced. AUROC remains interpretable, but the descriptor
  baseline and any enrichment metric both inflate on such a panel.
- The 0.70 descriptor threshold is a convention fixed in advance, not a calibrated boundary. Two
  panels sit near it (Factor Xa 0.711, CHEMBL612545 0.692) and would flip under a small change.
  Their marginal-value verdicts do not depend on the threshold.
- Four panels is the entire population on disk, not a sample. No selection, but equally no
  generalisation beyond this project's own docking.

## What this changes

Any standing claim resting on a docking-score difference from these panels should be read as
unsupported until re-derived. The two panels in the screening-decomposition manuscript (Mpro,
Factor Xa) are already reported there with exactly this framing, so the manuscript needs no
change. PD-L1 was previously **deferred**; it should now be described as **confounded and
closed**, because there is no docking result to recover from it.

The expected direction held: the audit reduced the number of standing claims and rescued none.
