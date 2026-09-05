# Boltz-2 on Mpro: the first method tested that adds real information over free descriptors

**4 September 2026.** 751 compounds dispatched, 750 scored, 1 failed. 20.7 GPU-hours on a
GCE L4 (projected 21.5), ~AUD 27. Analysis code committed `5afb5bec8` **before the data
existed**; reading rules in `PREREGISTRATION.md` (`9e3083d02` + amendments).

## Controls — all pass

| control | result | |
|---|---|---|
| **Positive** — reproduce the known descriptor baseline | **0.7654** vs 0.7653 expected | **PASS** |
| **Permutation** — shuffled labels | **0.4952**, inside [0.45, 0.55] | **PASS** |
| **Dropout** — 1 inactive of 751 | bounded 0.7897–0.7917 vs 0.7913 | **immaterial** |

The dropout rule as written fires VOID, but it is **degenerate at n = 1**: with a single
failure the active fraction among failures is 0% or 100% by construction and can never land
within 10pp of 34.3%. Rather than argue the rule, the impact was bounded — imputing the
dropped inactive at the best and worst possible scores moves AUROC by **0.002**. Immaterial.

(An earlier run of the analysis reported "21 missing" because it compared against the
771-compound manifest superset rather than the 751 dispatched. 20 of those have no SMILES and
were never in the Vina benchmark either. The AUROCs were computed on the correct set and are
unaffected.)

## Result

| method | AUROC |
|---|---|
| Vina, repaired receptor | 0.4530 |
| seven free descriptors | 0.7654 |
| **Boltz-2** `affinity_probability_binary` | **0.7913** |
| **descriptors + Boltz-2** | **0.8588** |

| pre-registered rule | Δ | 95% CI | verdict |
|---|---|---|---|
| **1.** Boltz-2 − descriptors > +0.04 | +0.0259 | [−0.0227, +0.0742] | **NOT DEMONSTRATED** |
| **3.** desc+Boltz-2 − descriptors > +0.04 | **+0.0934** | **[+0.0616, +0.1253]** | **SUPPORTED** |

**Boltz-2 alone does not beat the descriptor baseline by the pre-registered margin.** +0.0259
sits under the 0.04 bar with an interval spanning zero.

**But in combination it clears the bar decisively.** +0.0934 with a CI excluding zero, more
than twice the ~0.04 measurement floor. That means Boltz-2 carries information the descriptors
do not have — which is precisely what docking never managed:

| method | adds over descriptors |
|---|---|
| Vina (Mpro) | −0.0023 |
| Vina (Factor Xa) | +0.0150 (n.s.) |
| pose ensemble | +0.0194 (n.s.) |
| residual docking | reverses sign between targets |
| **Boltz-2** | **+0.0934, significant** |

After five candidate improvements that all landed inside the noise floor, this is the first to
clear it.

## Secondary endpoint is inverted

`affinity_pred_value` scores **0.2516** — far below chance. Sign-flipped it is 0.7484. It
carries signal with the wrong sign, exactly the pathology Vina shows on this target. Reported
as secondary and **not** substituted for the primary; the endpoint was fixed in advance
precisely to prevent that swap.

## What this does not show

> **Superseded 5 Sep 2026 — read this section with `ARM2_RESULT.md` beside it.** The
> lower-bound framing below was pre-registered and is **withdrawn**. Arm 2 supplied the pocket
> and scored 0.7783 against this blind run's 0.7913 (r = 0.952): running blind cost Boltz-2
> nothing, so +0.0934 needs no allowance and is not a lower bound.

The run was **blind** — no pocket conditioning — while the Vina comparator docked into a
defined 22 Å box on the known site. It also used `diffusion_samples = 1`, no `potentials`, and
no `affinity_mw_correction`. All recorded in the pre-registration on 4 Sep, before this number
existed. So the result is a **lower bound** on a handicapped configuration: Boltz-2 achieved
+0.0934 over descriptors while doing the harder version of the task.

The obvious next arm is pocket-conditioned with `affinity_mw_correction` enabled — cheap,
since the MSA is computed and cached. That last flag matters given this project's finding that
docking scores are dominated by compound size; if Boltz-2 ships a ligand-size correction, the
affinity head may currently carry the same bias.

## Scope

One target, one benchmark, one configuration. Mpro is a rigid, well-characterised protease of
the kind Boltz-2's own documentation lists as favourable, and its scaffolds are likely
represented in training data. This does not generalise to Factor Xa or PD-L1 without running
them.
