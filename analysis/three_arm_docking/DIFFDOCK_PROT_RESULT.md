# The DiffDock pose-source arm is not evaluable on this panel — and neither was the published version

**7 September 2026.** `rescore_diffdock_prot.py`, `analyse_diffdock_prot.py`. Rules in
`PREREGISTRATION_DIFFDOCK_PROT.md` (`ce636ea0a`) with Amendment 1. **The arm is VOID by its own
pre-registered dropout control. No effect was read.**

## What ran

The 750 pocket-conditioned DiffDock rank-1 poses were recovered from the original 31 August
Kaggle kernel outputs — no GPU re-run was needed, since protonation does not change DiffDock's
pose generation. Each was minimised and scored against the repaired receptor
(`md5 28657daa…`, `vina_md5` asserted equal to the recorded protocol).

**Coverage: 604/750 = 80.5%**, above the amended 75% floor. That is almost exactly the published
arm's own 604/751 = 80.4% — an unplanned but close replication of the original run's yield.

The 146 failures are not conversion errors. Vina rejects them with *"The ligand is outside the
grid box"*: even under pocket conditioning, DiffDock's top-ranked pose lands outside the 22 Å box
for about a fifth of compounds.

## Why it is void

| dropout control | result | |
|---|---|---|
| active fraction among excluded vs panel | **0.445 vs 0.342** (\|diff\| **0.103**, tolerance 0.10) | **FAIL** |
| mean molecular weight, excluded vs included | **474.6 vs 397.5 Da** | size-selected |

The exclusion is not random. **DiffDock systematically fails to place larger compounds in the
pocket, and larger compounds are disproportionately actives on this panel.** Dropping them
enriches the surviving subset in inactives, which biases any ranking comparison computed on it.

The active-fraction criterion fails narrowly — 0.103 against a 0.10 tolerance — and we are not
going to argue that 0.103 rounds to 0.10. The 77 Da mean-weight gap corroborates the same
conclusion independently, and the threshold was fixed before the numbers existed precisely so
that a marginal failure could not be talked away.

## This invalidates the published number too, not just the re-run

The published arm reported **+0.0173 [−0.0269, +0.0636]** on 604 of 751 compounds. That is the
same 80% subset, produced by the same mechanism, and therefore carries the same size- and
label-biased exclusion. The original run did not test for it.

So §3.2's pose-source row is not a conditional null awaiting a better receptor. **It is not
evaluable on this panel at all**, on either receptor, until DiffDock places the full compound set
inside the box — which is a property of DiffDock on this target, not of the analysis.

## Consequence for the manuscript

§3.2 currently reports the arm as **+0.017 [−0.027, +0.064], inside the floor**, and §4 lists it
as the last remaining conditional null. Both must change: the row becomes **not evaluable**, with
the reason stated. The paper loses a data point and gains an honest one.

This also removes the last of the four defective-receptor conditionals, though not in the way
intended — three arms were re-run or resolved, and the fourth turns out to have been
uninterpretable from the start.

## What the control was for

Amendment 1 added this dropout check because excluding a fifth of a panel non-randomly is a
failure mode this project has hit before. It fired on the first arm it was applied to. Had it not
been registered, the arm would have produced a plausible small effect on a biased subset — which
is exactly what the published version did.
