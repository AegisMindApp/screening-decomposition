# Our clash-rate measurement reproduces a published benchmark — and there is a proper tool for it

**31 August 2026.** Assessment of a suggested reading list on modern docking, checked rather than
accepted.

## The one item that matters: PoseBusters

We found that pocket-constrained DiffDock poses score as steric clashes **38.9%** of the time
(63% blind) against **0%** for Vina, and improvised the check — a positive Vina affinity implies
a clash.

That is a crude version of an established test. **PoseBusters** ([Chemical Science
2024](https://pubs.rsc.org/sc/article/15/9/3130/827511/PoseBusters-AI-based-docking-methods-fail-to))
does it properly: chemical validity, intramolecular geometry (bond lengths, angles, aromatic ring
planarity) and intermolecular checks (protein–ligand clash, volume overlap).

Its published finding is titled *"AI-based docking methods fail to generate physically valid
poses"*, and it reports deep-learning methods producing **up to 95% invalid poses despite RMSD
< 2 Å** — DiffDock among the five tested, against AutoDock Vina and GOLD as the classical
controls.

**So our 39–63% is not a symptom of our setup.** It independently reproduces a benchmarked
property of generative docking, at the milder end of the published range. That materially raises
confidence in the arm-2 decision: withholding the AUROC was the right call, not excessive caution.

**Action:** run PoseBusters over the pocket-constrained poses rather than relying on the
positive-score proxy. It distinguishes *why* a pose fails — clash versus bond geometry versus
stereochemistry — which the proxy cannot, and that distinction decides whether local minimisation
can rescue the arm or not.

## What else is worth taking

**Co-folding is a genuinely different class**, not another docking method: AlphaFold 3, Chai-1,
Boltz-1/2 predict the complex from sequence. **Boltz-2 predicts affinity explicitly**, which is
the quantity our benchmark measures — unlike DiffDock's confidence score, which is not an affinity
estimate and benchmarks at AUROC 0.538–0.56 for screening. If anything in this list deserves a
future arm, it is Boltz-2, and on the affinity head rather than the pose head.

**GNINA** is recommended there as Vina search plus CNN rescoring. We have already measured exactly
that: **+0.140 AUROC** on identical poses. The suggestion is sound and is already banked.

## What to treat carefully

The claim that co-folding methods "generally outperform" docking rests on **PoseBench**, which
scores **pose prediction** from apo structures. Our failure is in **ranking actives above
inactives**, which is a different quantity — a method can place a ligand beautifully and still
rank a library badly. The email's own closing point is the honest one and matches what we
measured: classical methods lead on physical validity, generative ones on speed and flexibility.

DynamicBind, FABFlex, ArtiDock and the specialised models are real but do not address our
bottleneck. We have now measured, on one benchmark with pre-registered rules, that the **scoring
function** dominates (+0.140), **receptor preparation** is second (+0.045) and **search effort**
does not matter (−0.019). Every method in that list is a better *sampler*. None of them is a
better *scorer*, which is the thing that moved our number most.

## Net

One concrete tool adopted (PoseBusters), one candidate for a future arm (Boltz-2 affinity), one
recommendation confirmed as already done (GNINA), and a useful external check that today's clash
rate is a known property rather than a mistake of ours.
