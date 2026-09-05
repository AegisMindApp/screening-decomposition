# screening-decomposition

Pre-registration, analysis code and raw outputs for **"What actually moves virtual screening
performance: two measured resolution limits and a pre-registered decomposition of six
interventions"**.

The manuscript is [`docs/papers/screening_decomposition_preprint.md`](docs/papers/screening_decomposition_preprint.md);
supplementary tables S1–S5 are [`docs/papers/screening_decomposition_supplementary.md`](docs/papers/screening_decomposition_supplementary.md).

## Why this repository exists

The paper's claim is that each reading rule — endpoint, threshold, control, abort condition —
was fixed **before** the corresponding data existed. That claim is only worth anything if you
can check it, so the history is here rather than a snapshot.

| arm | pre-registration | commit | date |
|---|---|---|---|
| Scoring function (gnina rescore) | `analysis/three_arm_docking/PREREGISTRATION.md` | `f5e42306a` | 2026-08-31 08:46 |
| Receptor preparation repair | `analysis/receptor_prep/PREREGISTRATION.md` | `05d9f78ec` | 2026-08-31 14:26 |
| Pose ensemble vs top pose | `analysis/pose_ensemble/PREREGISTRATION.md` | `d5d08c1ac` | 2026-09-02 22:52 |
| Pose generator (DiffDock-L) | `analysis/three_arm_docking/PREREGISTRATION.md` | `f5e42306a` | 2026-08-31 08:46 |
| Boltz-2 Mpro arm 1 (blind) | `analysis/boltz2/PREREGISTRATION.md` | `8df29ba21` | 2026-09-03 08:27 |
| Boltz-2 analysis script | `analysis/boltz2/analyse_boltz2.py` | `42f94969b` | 2026-09-03 13:51 |
| Boltz-2 Mpro arm 2 (pocket) | `analysis/boltz2/PREREGISTRATION_ARM2.md` | `064969caf` | 2026-09-04 11:28 |
| Boltz-2 Factor Xa | `analysis/boltz2/PREREGISTRATION_TARGETS.md` | `00baa6db9` | 2026-09-04 11:40 |

Check any of them:

```
git log -1 --format=%ad --date=iso f5e42306a
```

**What this does and does not establish.** These are commit dates carried through from the
repository this was exported from — git metadata, not a third-party timestamp. They are our
record. A reader can confirm internal consistency and ordering; a reader cannot, from this
repository alone, rule out that the metadata was authored. Two comparisons — the search-effort
arm and the LIT-PCBA descriptor baseline — were **not** pre-registered and say so in S1.

## Reproducing the headline result

The resolution floor is the paper's lead claim and runs in a few seconds with no dependencies
beyond the standard library:

```
python3 analysis/docking_value/floor_interval.py
```

It asserts a positive control before reporting anything — it must reproduce the previously
published numbers (r = 0.9369, AUROC 0.4182 / 0.3793) or it exits. That control is not
decorative: it **failed** on 6 September 2026 and caught a real defect, written up in
[`analysis/docking_value/FLOOR_CORRECTION.md`](analysis/docking_value/FLOOR_CORRECTION.md).
The floor had been reported as a comparison of identical poses; it was not.

The LIT-PCBA descriptor baseline needs RDKit, scikit-learn and the LIT-PCBA distribution
(not redistributed here — see [litpcba.drugdesign.fr](https://drugdesign.unistra.fr/LIT-PCBA)):

```
python3 analysis/litpcba/descriptor_baseline.py
```

## Read this too

[`docs/papers/AUDIT_screening_decomposition.md`](docs/papers/AUDIT_screening_decomposition.md)
is the pre-submission audit of our own manuscript. It records five blockers found against the
draft, including a published number of ours that was wrong (we cited AutoDock Vina on LIT-PCBA
at 0.61; the primary source says **0.581** — an error that, corrected, widens our own margin).
It is here because an audit that only appears when it flatters the authors is not an audit.

## Layout

```
analysis/docking_value/     resolution floors, marginal value of docking, floor correction
analysis/three_arm_docking/ scoring function and pose-source arms
analysis/receptor_prep/     the hydrogen-bond donor defect and its repair
analysis/pose_ensemble/     pose-ensemble vs top-pose scoring
analysis/boltz2/            Boltz-2 co-folding arms, both targets, and the size-bias analysis
analysis/litpcba/           descriptor baseline on LIT-PCBA, full and AVE-debiased
docs/papers/                manuscript, supplementary, pre-submission audit
```

## Licence

Code is MIT (`LICENSE`). The manuscript, supplementary material and result data are
CC BY 4.0. LIT-PCBA and ChEMBL data carry their own upstream terms and are not redistributed
here beyond the derived scores needed to reproduce the analyses.

## Citation

Goodman J (2026). *What actually moves virtual screening performance: two measured resolution
limits and a pre-registered decomposition of six interventions.* Preprint.
