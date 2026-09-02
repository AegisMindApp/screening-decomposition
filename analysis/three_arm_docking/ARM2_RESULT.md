# Arm 2 result: pose sampling is not the failure — supported

**3 September 2026.** Reading rule pre-registered 31 August, before any DiffDock run:
*substituting DiffDock-L pose generation while retaining Vina scoring will change AUROC by
less than 0.05*; refutation requires the ΔAUROC(2−1) **CI entirely above +0.05**.

## Both controls that failed the previous attempt now pass

- **Box-size control**: 40 arm-1 poses scored at both 22 Å and 30 Å, **max |diff| = 0.0000**.
  Enlarging the box does not move scores for ligands that already fitted, so the enlargement
  needed to stop dropouts did not itself become a confound.
- **Dropouts**: 753/753 scored, **0 failures**. The previous run lost 149 compounds, and that
  loss was size-selected (43.6% actives among failures vs 31.8% among survivors).

## Result

| arm | pose source | scoring | AUROC (minimised) | AUROC (raw) |
|---|---|---|---|---|
| 1 | Vina | Vina | 0.4024 | 0.4038 |
| 2 | **DiffDock-L** | Vina | **0.4197** | **0.4599** |

| field | ΔAUROC (2−1) | 95% CI | CI entirely above +0.05? | verdict |
|---|---|---|---|---|
| minimised | **+0.0173** | [−0.0269, +0.0636] | no | **supported** |
| raw | +0.0561 | [−0.0064, +0.1177] | no | **supported** |

Minimised is the primary comparison: raw DiffDock poses carry clashes (39.7% positive before
minimisation, 1.6% after), so raw scores are contaminated by clash artefacts.

**Both fields support the hypothesis.** Note the raw point estimate (+0.0561) exceeds 0.05
while its interval spans zero — under the rule as written, refutation needs the *interval*
above +0.05, not the point estimate, so this is not a refutation. Recording this explicitly
because a stricter reading would have flipped the verdict, and the pre-registered wording is
the authority.

## The number that carries the most weight

**Pose-source score correlation r = 0.139** (minimised; −0.088 raw). The two pose sets are
essentially unrelated — DiffDock and Vina put the ligands in different places — and they still
produce the same chance-level ranking (0.4024 vs 0.4197). Different poses, same answer. That
is stronger evidence that sampling is not the bottleneck than the ΔAUROC alone.

Both arms remain far below the descriptor bar of **0.7653**.
