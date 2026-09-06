# Every docking panel this project has run fails at least one bias control

**6 September 2026.** `sweep.py`, results in `SWEEP.json`. Rules fixed in `PREREGISTRATION.md`
(`e5c305113`) before any unaudited panel was scored; Amendment 1 corrects a control copied from
the wrong panel. Gate **PASSES** on all four positive controls.

This is not an attempt to rescue anything. It applies three diagnostics developed this week —
descriptor baseline, marginal value, residual AUROC — to every panel on disk, to establish which
historical docking claims were readable at all.

## Result: 0 of 4 panels survive

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

**CHEMBL612545 is the only panel whose residual range includes 0.50** ([0.4571, 0.5026]). The part
of the docking score orthogonal to ligand properties does not rank actives at all.

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
