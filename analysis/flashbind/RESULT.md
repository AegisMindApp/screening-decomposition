# FlashBind on the Mpro panel: the third scoring approach does not clear the ceiling

**Run 8 September 2026.** Rules fixed in `PREREGISTRATION.md` (c6b4504f9, amended 24bc26573
and 453e2031f) before any score existed. Kernel `oceansparx/flashbind-probe` v15, n = 750.

## Pre-registered controls — all three pass

| control | required | measured (n = 750) | verdict |
|---|---|---|---|
| 1. descriptor baseline reproduces the panel | 0.7654 ± 0.005 | **0.7664** [0.7578, 0.7714], 12 seeds | **PASS** |
| 2. coverage of 751 | ≥ 90% | **750/751 = 0.999** | **PASS** |
| 3. shuffled labels | AUROC in [0.45, 0.55] | **0.4978** [0.4583, 0.5416] | **PASS** |

Control 1 failed on first attempt at 0.8982 because I ran gradient boosting under `KFold`; the
published baseline is `StandardScaler + LogisticRegression` under `StratifiedKFold(5)`. The
control was wrong, not the panel. A positive control has to reproduce the *procedure* as well
as the data, or it tests nothing.

## Pocket agreement — the confound is ruled out, on the full panel

Full panel, n = 750: median centroid distance **8.00 Å** from the docking box centre
[9.05, 8.90, −1.51], range 7.53–16.67, **`frac_in_box` 0.999** — one compound of 750 places its
pocket centroid outside the 22 Å box. Reproduced at 8.00 / 8.01 / 8.00 Å across v13, v14 and
v15. (The 48-compound probe gave 7.88–7.90 Å with max 8.62; the probe is a subset of the panel
and is not quoted as independent confirmation.)

FABind+ picks its own pocket and lands in the same site the six docking scores were confined
to, so the residual reads as a property of the scoring approach rather than of being
unconstrained by the box.

## Result

| quantity | FlashBind | Boltz-2 | docking comparators |
|---|---|---|---|
| raw AUROC (Mpro arms) | **0.4185** [0.3698, 0.4674] | 0.7913 | 0.3791 – 0.5389 |
| residual AUROC, GB, 12 seeds (pre-registered band, incl. FXa) | **0.5384** [0.5287, 0.5456] | 0.6567 [0.6439, 0.6648] | 0.4941 – 0.5733 |
| — the same band, Mpro arms only | 0.5384 | 0.6567 | 0.4941 – 0.5426 |
| marginal value over 7 descriptors, 12 seeds | **+0.0247** [+0.0217, +0.0271] | +0.093 [+0.062, +0.125] | — |

**The residual sits inside the pre-registered docking band**, with the full 12-seed range
inside it. The arm does not support generalising "the ceiling is a property of the scoring
approach". Boltz-2 remains the only point outside the band, and one point is not a
generalisation.

The pre-registered bar of 0.573 is set by the Factor Xa arm, so it is a cross-target band. On
an Mpro-only band (0.4941–0.5426) FlashBind's 0.5384 is still inside. The verdict is the same
under either definition; the band was not re-opened after seeing the result.

**The marginal value is real but small.** +0.0247 with a 12-seed range that does not straddle
zero, and whose upper end (+0.0271) sits clear below the lower end of Boltz-2's +0.093
[+0.062, +0.125]. Reported over all 12 fold seeds: quoting one seed would have measured the
increment against whichever baseline draw suited it, and fold-seed movement on this panel is
0.01–0.03 — the size of the effect itself.

## The raw AUROC, in context

FlashBind's 0.4185 is below chance with 0.5 outside the bootstrap interval, and both
checkpoints agree independently (0.4445, 0.4067). **But it is mid-pack among the comparators on
this panel, not an outlier:** Vina repaired 0.4530, Vina raw 0.4079, gnina Vina term 0.3791,
gnina CNNscore 0.5041, gnina CNNaffinity 0.5389. Three of the five docking scores are also
below chance here. FlashBind fails on this panel the way the docking family fails on it, which
is the manuscript's own thesis rather than a distinguishing defect of this method.

## Why it inverts — and the limit of that explanation

The score is anti-correlated with size: rotatable bonds **−0.557**, molecular weight −0.444,
TPSA −0.346, HBD −0.278 (logP +0.270). On this panel the actives are the *larger* compounds
(Spearman(label, MolWt) = +0.197), since Mpro inhibitors are largely peptidomimetic.

That is not the whole story, and the first draft of this file over-claimed it. Inverse
molecular weight *alone* scores **0.3801** and inverse rotatable-bond count **0.3741** — both
further below chance than FlashBind's 0.4185. So the model is not simply a size proxy: it ranks
somewhat better than raw size inversion, and after descriptors are regressed out its residual
is 0.5384, above 0.5. There is non-size signal in it. There is just no more of it than docking
carries.

## What this does and does not say

On **one target**, against a panel whose actives skew large, FlashBind's binary probability
does not rank actives and carries no more target-specific signal than docking once descriptors
are regressed out.

It does **not** say the method is poor at its own task. The binary head is trained on MF-PCBA
HTS data, where "active" is defined differently from a ChEMBL potency threshold, and the
leakage check found no protease assay in any of those 105 AIDs — so this panel is out of domain
in both target and label definition. A single target cannot support a claim about the method in
general, and none is made here.

Two disclosed asymmetries: FABind+'s pose generator is PDBbind-trained and PDBbind contains
Mpro complexes, which advantages placement here; and Mpro is a homodimer, so the protein
representation had to be built per chain and concatenated (`RECEPTOR_MISMATCH.md`).

## The open question stays open

§3.4 asked for a third scoring approach clearing the ceiling. This one does not. Whether the
ceiling is a property of docking-family scoring or a property of this panel is still answered
by exactly one point.
