# gnina on Factor Xa: the paper's largest intervention reverses sign on a second target

24 September 2026. Answers the gnina arm of `PREREGISTRATION.md` (amendment of the same day,
written before any Factor Xa intervention data existed).
Re-derive: `python analysis/replication_fxa/gnina_fxa_verdict.py` → `GNINA_FXA_RESULT.json`.

## Verdict: DOES NOT REPLICATE

| | Mpro | **Factor Xa** |
|---|---|---|
| Vina poses + Vina score | 0.3992 | 0.6872 |
| **same poses** + gnina CNN | 0.5389 | 0.6289 |
| ΔAUROC | **+0.1397** [+0.0913, +0.1867] | **−0.0583** [−0.0930, −0.0235] |
| floor charged | Mpro exh-4 | Factor Xa exh-4, 0.00718 at n=849 |

Both gates passed before any AUROC was read:

- **coverage 95.8%** (849 of 886; 37 gnina failures) against a pre-registered floor of 90%
- **positive control r = 0.9934** — gnina's own Vina-style affinity for each pose against our
  cached Vina score. Mpro gave +0.970 against the same 0.90 threshold. gnina was demonstrably
  scoring the poses we handed it.

The poses themselves are the ones already inside the paper's Factor Xa floor: all 886 Vina scores
are **bit-identical** to the committed Kaggle seed-1 run, because Vina with a fixed seed and
`--cpu 1` is deterministic.

## Why this is a sign reversal and not the attenuation we predicted

The amendment recorded, before the run, that a **smaller** gnina gain was expected on Factor Xa:
Vina sits at 0.6872 there and 0.3992 on Mpro, so there is less headroom above the baseline. It
then fixed the test as **sign plus floor-clearing**, explicitly so that a shrunken effect could
not be waved through as "expected attenuation".

The effect does not shrink. It **reverses**, by −0.058, with its interval clear of the floor and
nowhere near Mpro's. The pre-registered rule fails on the sign, which is the term that was fixed
in advance precisely because it could not be rationalised afterwards.

## What it suggests, stated as a hypothesis rather than a finding

On Mpro, Vina ranks **worse than random** (0.3992) and gnina's CNN lifts it across chance to
0.5389. On Factor Xa, where Vina already ranks usefully (0.6872), the same rescoring degrades it.

The reading consistent with both numbers is that **+0.140 measured a recovery from a broken
baseline rather than a property of the scoring function**. A method that replaces a
worse-than-random ordering with a near-random one gains a lot of AUROC while producing nothing a
screen could use; applied where the ordering already carries signal, it destroys some.

That is a hypothesis this experiment does not settle — two targets cannot separate "gnina helps
bad baselines" from "gnina helps Mpro". It is offered because it explains the sign reversal
without either result being wrong, and because it is testable on a third target.

## What this does to the manuscript

§3.3 reports gnina as its cleanest SURVIVES, at +0.140 with an interval far clear of the floor.
That measurement stands: it is correct for Mpro and nothing here disputes it. What fails is the
implied generality — the row is presented as "changing the scoring function on identical poses"
helping, and on a second target it hurts.

The manuscript must not describe this intervention as replicated, and should not report the Mpro
row without the Factor Xa row beside it. **Deciding how §3.3 and the abstract are reworded is a
judgement for the author**, not something to fold in silently, so the wording is left unchanged
here pending that decision.

Note also that this cuts against the paper's own framing in a useful direction. The paper argues
resolution limits are underappreciated; this is a case where an effect **five times the floor**,
resolved beyond any doubt on its own panel, still fails to generalise. Clearing the resolution
bar is necessary and plainly not sufficient, and the paper is stronger for saying so about its
own headline number.

## Scope

One seed of Factor Xa poses (seed 1), one receptor, gnina v1.3 `--score_only --cnn_scoring
rescore --seed 42` — identical flags and binary version to the Mpro arm. Pose-to-pose noise is
shared between the two arms by construction, since both score the same files.
