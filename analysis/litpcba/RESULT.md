# Seven free descriptors on LIT-PCBA: 0.726 median, against published Vina at 0.61

**5 September 2026.** Our finding that docking loses to seven physicochemical descriptors
rested on two panels we assembled ourselves, either of which could be confounded — PD-L1
demonstrably is. LIT-PCBA is the field's standard "unbiased" benchmark, built explicitly to
keep actives and inactives within similar property ranges. This tests the claim on someone
else's data against someone else's docking numbers.

Descriptors: MW, clogP, HBD, HBA, rotatable bonds, TPSA, formal charge. Logistic regression,
out-of-fold on the full set; on the AVE variant the authors' own train/validation split is
honoured rather than re-folded. Inactives capped at 25,000 per target (seeded).

## Result

| | median AUROC | targets |
|---|---|---|
| **descriptors, LIT-PCBA full** | **0.7260** | 15 |
| descriptors, full, ≥50 actives | 0.6830 | 10 |
| **descriptors, AVE-debiased** | **0.6510** | 15 |
| descriptors, AVE, ≥50 actives | 0.6510 | 5 |
| *published* **AutoDock Vina**, LIT-PCBA | *0.61* | 15 |
| *published* **GNINA**, LIT-PCBA | *0.61–0.62* | 15 |

Published figures from arXiv 2605.01681 (Vina, line 135) and Sunseri & Koes 2021 (GNINA).

**On the full set, like for like, seven free descriptors beat published Vina by 0.116** — or
0.073 restricting to the ten targets with ≥50 actives. That replicates our Mpro and Factor Xa
result on the field's own benchmark, using docking numbers we did not produce, so it cannot be
attributed to our receptor preparation, box definition or panel construction.

**It also validates our panels.** Our descriptor baselines — Mpro 0.7654, Factor Xa 0.7041 —
sit inside LIT-PCBA's 0.577–0.864 range. They are typical, not unusually confounded.

## The AVE result cuts partly against the headline

AVE debiasing drops descriptors from 0.726 to 0.651. **So the full set does retain property
bias, and roughly 0.075 of the headline number is that bias.** The AVE procedure removes some
of it — and not all: 0.651 on a deliberately debiased set is still well clear of chance, from
seven trivial properties.

The median is not an artefact of small validation sets. Restricting to targets with ≥20 or
≥50 actives moves the AVE median by at most 0.002 (0.6510 → 0.6497 → 0.6510), even though
several AVE validation sets are tiny (ESR1_ago 3 actives, ADRB2 4, PPARG 6).

## The comparison we cannot make

**There is no published Vina number on the AVE-debiased set.** Comparing our AVE descriptor
result (0.651) to a full-set Vina number (0.61) would not be like for like, and we do not make
that comparison. Vina would presumably also fall under AVE debiasing; by how much is unknown
and out of reach here — docking 15 targets against hundreds of thousands of compounds is not
affordable. The honest position is that the descriptors-beat-docking comparison stands on the
full set only.

## Bug found and fixed

The first AVE run reported an empty table and **exited 0**. LIT-PCBA ships two layouts —
`actives.smi`/`inactives.smi` for the full set, `active_T/V.smi` + `inactive_T/V.smi` for AVE —
and the script recognised only the first, skipped every target, and looked like a clean run.
It now detects the layout, honours the train/validation split when present, and exits with a
fatal error naming the layouts it understands rather than printing nothing.
