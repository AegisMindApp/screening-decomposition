# Pre-registration — Boltz-2 arm 2: pocket-conditioned

**Written 4 September 2026, before arm 2 runs.** Arm 1 (blind) is complete and reported in
`RESULT.md`: Boltz-2 0.7913, descriptors 0.7654, desc+Boltz-2 0.8588.

## What changes, and what does not

**Exactly one variable changes: pocket conditioning.** Arm 1 ran blind — sequence and SMILES
only — while the Vina comparator docked into a 22 Å box centred on the known site. Arm 2
supplies that same site to Boltz-2 as a pocket constraint.

`diffusion_samples` stays at 1, `potentials` stays off, and `affinity_mw_correction` stays
off. Changing more than one thing would make any improvement unattributable, which is the
error this project has spent the week documenting in other people's benchmarks.

`affinity_mw_correction` is deliberately **not** enabled here. The question it addresses —
whether Boltz-2's affinity head carries the ligand-size bias we documented in Vina — can be
answered from arm 1's existing scores at zero compute cost by regressing them on molecular
weight. That analysis is separate and does not need a run.

## Pocket definition — derived, not guessed

Residues with any atom within 8.0 Å of the Vina box centre `[9.05, 8.9, -1.51]`, chain A:

    A164 HIS  A145 CYS  A41 HIS  A189 GLN  A165 MET  A143 GLY
    A142 ASN  A49 MET   A144 SER A166 GLU  A141 LEU  A26 THR

**Sanity check passed:** the set contains the catalytic dyad **Cys145 (4.41 Å)** and
**His41 (5.02 Å)**, and the canonical substrate-binding residues His163/His164/Glu166/Met165.
Chain B's equivalents are 41–49 Å away, confirming the site sits on the chain-A protomer only.
`max_distance: 5.0`.

## Reading rules — fixed now

Primary endpoint remains `affinity_probability_binary`. Same compounds, same folds, same
descriptor baseline (0.7654), out-of-fold only.

1. **Does conditioning beat the free baseline?** ΔAUROC(arm2 − descriptors) **> +0.04** with a
   paired-bootstrap 95% CI excluding zero. Arm 1 scored +0.0259 and failed this.
2. **Does conditioning beat blind?** ΔAUROC(arm2 − arm1) **> +0.04** with CI excluding zero.
   This is the question the arm exists to answer. Anything less is reported as **pocket
   conditioning did not help by a resolvable margin**, not as a trend.
3. **Combination.** desc+arm2 versus descriptors alone, same margin. Arm 1 scored +0.0934 and
   passed; arm 2 must be compared against arm 1's combination, not just against descriptors.

Controls unchanged and all must pass before any Boltz-2 number is read: positive control
reproducing 0.7654 ± 0.005, permutation inside [0.45, 0.55], dropout active-fraction within
10pp — with the standing note that the dropout rule is **degenerate below ~5 failures** and
impact should be bounded rather than the rule argued.

## Cost

~21 GPU-hours, ~AUD 27, inside the 25.0 GPU-hour ceiling authorised on 3 Sep and the AUD 100
project budget. Cumulative Boltz-2 spend after this arm: ~AUD 54.

## What would make this arm a null

Rule 2 failing. If supplying the correct binding site does not move the answer by more than
the measurement floor, that is a substantive finding in itself: it would put pocket
conditioning alongside pose source, search effort and receptor preparation in the set of
interventions that do not survive contact with this benchmark.
