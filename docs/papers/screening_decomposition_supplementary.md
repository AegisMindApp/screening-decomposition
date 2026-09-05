# Supplementary material

**S1. Pre-registration index.** Every reading rule — endpoint, threshold, control, abort
condition — was committed before the corresponding data existed. Each row gives the file, the
commit that introduced it, and that commit's date. Amendments are appended to the file rather
than editing the original text, with the reason recorded in place.

| arm | pre-registration file | commit | committed |
|---|---|---|---|
| Scoring function (gnina rescore of Vina poses) | `analysis/three_arm_docking/PREREGISTRATION.md` | `cd8f55639` | 2026-08-31 08:46 |
| Receptor preparation repair | `analysis/receptor_prep/PREREGISTRATION.md` | `1f68b20c4` | 2026-08-31 14:26 |
| Pose ensemble vs top pose | `analysis/pose_ensemble/PREREGISTRATION.md` | `448d954e5` | 2026-09-02 22:52 |
| Pose generator (DiffDock-L) | `analysis/three_arm_docking/PREREGISTRATION.md` | `cd8f55639` | 2026-08-31 08:46 |
| Boltz-2, Mpro, arm 1 (blind) | `analysis/boltz2/PREREGISTRATION.md` | `9e3083d02` | 2026-09-03 08:27 |
| Boltz-2 analysis script, written before the data | `analysis/boltz2/analyse_boltz2.py` | `5afb5bec8` | 2026-09-03 13:51 |
| Boltz-2, Mpro, arm 2 (pocket-conditioned) | `analysis/boltz2/PREREGISTRATION_ARM2.md` | `e96d4fea1` | 2026-09-04 11:28 |
| Boltz-2, Factor Xa | `analysis/boltz2/PREREGISTRATION_TARGETS.md` | `68bd9b94a` | 2026-09-04 11:40 |

Commit dates above are git metadata carried through from the project's own history; they are not
a third-party timestamp, and a reader should treat them as the authors' record rather than as
independent verification of ordering.

Bars citing a "≈0.04 measurement floor" — the pose-ensemble arm and both Boltz-2 arms — were fixed
before the floor correction described in §3.1 and S3. They are read as written; see §3.2.

The search-effort comparison (exhaustiveness 4 vs 32) and the LIT-PCBA descriptor baseline were
not separately pre-registered; both are reported here as such. The search-effort arm reuses two
runs that pre-dated the pre-registration programme, and the LIT-PCBA baseline applies a method
already fixed for the Mpro and Factor Xa panels to a third-party dataset.

**S2. Per-arm experimental conditions.** The six interventions are single-factor comparisons
against differing baselines, not a partition of one pipeline. Four were run before the receptor
defect described in §3.2 was found and repaired.

| intervention | baseline AUROC | receptor | exhaustiveness | n | poses held fixed? |
|---|---|---|---|---|---|
| Scoring function (gnina CNN) | 0.3992 | unprotonated | 4 | 743 | **yes** |
| Receptor repair | 0.4079 | unprotonated → protonated | 4 | 751 | receptor coordinates fixed; ligands re-docked |
| Pose ensemble | — | unprotonated | 4 | 751 | **yes** (same run, 5 modes) |
| Pose generator (DiffDock-L) | — | unprotonated | 4 (Vina arm) | 751 | no |
| Pocket conditioning (Boltz-2) | 0.7913 | none — co-folding from sequence | n/a | 750 | n/a |
| Eight-fold search effort | 0.4079 (exh 4) | unprotonated | 4 → 32 | 753 | no |

**S3. Resolution floor, both variants.** `analysis/docking_value/floor_interval.py`,
`FLOOR_INTERVAL.json`. Paired bootstrap over compounds, 10,000 resamples, unprotonated receptor.

| | value | 95% CI | n | r between the two score vectors |
|---|---|---|---|---|
| Protocol floor — exhaustiveness-32 Vina cache vs gnina `--score_only` on the exhaustiveness-4 poses | 0.0393 | [+0.0235, +0.0557] | 745 | 0.9372 |
| Scoring floor — exhaustiveness-4 Vina cache vs gnina `--score_only` on the same exhaustiveness-4 poses | 0.0201 | [+0.0112, +0.0286] | 743 | 0.9700 |

Neither score vector contains a compound differing from its counterpart by more than 2 kcal/mol.
An earlier description of the protocol floor as a comparison of "identical poses" was incorrect
and is corrected in `analysis/docking_value/FLOOR_CORRECTION.md`.

**S4. LIT-PCBA per-target results.** `analysis/litpcba/descriptor_baseline_full.json` and
`descriptor_baseline_ave.json`.

| target | full: actives | inactives | AUROC | AVE: actives | inactives | AUROC |
|---|---|---|---|---|---|---|
| ADRB2 | 17 | 25,000 | 0.7008 | 4 | 25,000 | 0.6016 |
| ALDH1 | 7,168 | 25,000 | 0.6169 | 1,343 | 25,000 | 0.6123 |
| ESR1_ago | 13 | 5,583 | 0.7260 | 3 | 908 | 0.8737 |
| ESR1_ant | 102 | 4,948 | 0.7525 | 25 | 794 | 0.6497 |
| FEN1 | 369 | 25,000 | 0.7596 | 91 | 25,000 | 0.7471 |
| GBA | 166 | 25,000 | 0.7456 | 41 | 25,000 | 0.7390 |
| IDH1 | 39 | 25,000 | 0.7921 | 9 | 25,000 | 0.7838 |
| KAT2A | 194 | 25,000 | 0.6656 | 48 | 25,000 | 0.6032 |
| MAPK1 | 308 | 25,000 | 0.7285 | 77 | 15,250 | 0.7129 |
| MTORC1 | 97 | 25,000 | 0.5945 | 24 | 8,243 | 0.5686 |
| OPRK1 | 24 | 25,000 | 0.8640 | 6 | 25,000 | 0.7672 |
| PKM2 | 546 | 25,000 | 0.6743 | 136 | 25,000 | 0.6068 |
| PPARG | 27 | 5,211 | 0.8459 | 6 | 861 | 0.7195 |
| TP53 | 79 | 4,168 | 0.5769 | 16 | 795 | 0.5050 |
| VDR | 884 | 25,000 | 0.6918 | 165 | 25,000 | 0.6510 |
| **median** | | | **0.7260** | | | **0.6510** |

**S5. Panels.** Mpro: 751 compounds, 257 active (34.2%), PDB 7VU6. Factor Xa: 886 compounds, 405
active (45.7%). The Factor Xa marginal-value comparison in §3.4 uses an earlier 854-compound
freeze of the same panel; where a number from that freeze is quoted it is labelled.
