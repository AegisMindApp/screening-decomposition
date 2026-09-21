# Stage 1 — the run-to-run resolution limit of the Mpro docking protocol

**AUROC seed SD = 0.00497 over 10 replicates; 95% upper bound = 0.00817, and that bound is the
floor.** 755 compounds (257 active / 498 inactive), exhaustiveness 4, everything fixed but
`--seed`. Design in [`PREREGISTRATION.md`](PREREGISTRATION.md), committed before any replicate
ran, with the comparison table fixed before this number existed.

The manuscript's two "limits" were **0.020** and **0.039** — 2.4× and 4.8× too large, and
neither measured what the decision rule needed.

## Three findings, in order of how much they matter

### 1. The paper's effects are differences between two single runs, and its CIs contain no seed variance

Its intervals bootstrap **compounds**. Each intervention compares two conditions run **once
each**, so the relevant noise is √2 × the single-run SD = **0.0116**, a 95% half-width of
0.0227. Adding that in quadrature — the minimum correction, not double counting:

| intervention | effect | published CI | with seed variance | verdict |
|---|---|---|---|---|
| Scoring function (gnina CNN) | +0.140 | [+0.091, +0.187] | [+0.087, +0.193] | **survives** |
| Pose ensemble | +0.055 | [+0.043, +0.068] | [+0.029, +0.081] | **survives** |
| Receptor preparation repair | +0.045 | [+0.024, +0.067] | [+0.014, +0.076] | **survives** |
| Binding site to Boltz-2 | −0.013 | [−0.026, +0.000] | [−0.039, +0.013] | no longer significant |
| Eight-fold search effort | −0.019 | [−0.032, −0.006] | **[−0.045, +0.007]** | **no longer significant** |

The three positive findings survive. The **eight-fold search-effort result does not** — and it
is the one the paper reports as a real, CI-excludes-zero effect.

### 1b. Which arms carry the noise, and why √2 is not uniform

The manuscript uses four methods, and they are not seeded alike:

| worker | method | seeded |
|---|---|---|
| `kaggle_shard_worker.py` | AutoDock Vina — **the baseline** | **no** |
| `minimise_vina_worker.py` | Vina minimisation | **no** |
| `colab_worker.py` | Vina (Colab variant) | **no** |
| `gnina_rescore_worker.py` | gnina CNN rescoring | `--seed 42` |
| `diffdock_worker.py` | DiffDock-L | `--seed 42` |
| `minimise_rescore_worker.py` | rescoring | `--seed 42` |

**Every Vina worker is unseeded; every non-Vina method is seeded.** All the run-to-run noise
therefore sits in the baseline each intervention is measured against. So the √2 factor applies
only where two *independent* Vina runs are compared; where gnina rescores the same poses, the
pose noise is common to both arms and largely cancels:

| intervention | pose noise | corrected CI | |
|---|---|---|---|
| Scoring (gnina CNN) | shared | [+0.089, +0.191] | survives |
| Pose ensemble | shared | [+0.035, +0.075] | survives |
| Receptor repair | independent | [+0.014, +0.076] | survives |
| Boltz-2 site | independent | [−0.039, +0.013] | not significant |
| Eight-fold search | independent | [−0.045, +0.007] | **not significant** |

My first pass applied √2 uniformly, which overstated the noise on the two shared-pose arms.
They survive either way, and the two that fail are both independent-run comparisons where √2 is
correct — so the conclusion is unchanged and the reasoning is now right rather than merely
conservative.

### 2. That effect is the size of the spread between replicates of one of its own conditions

Ten exh=4 replicates span **0.0168** AUROC. The claimed exh32 − exh4 difference is **0.0185**.
A ratio of **1.10**. Comparing two conditions measured once each is barely distinguishable from
comparing two seeds of the *same* condition.

### 3. Eight-fold more search moves the ranking LESS than a reseed does

Independent of AUROC, measured on the 742 compounds common to their exh=32 run, their exh=4 run
and all ten of my exh=4 seeds:

| comparison | Spearman ρ |
|---|---|
| their exh=32 vs their exh=4 — **8× more search** | **0.9593** |
| my exh=4 vs my exh=4 — **same protocol, different seed** | **0.9551** |
| their exh=4 vs my exh=4 | 0.9543 |
| their exh=32 vs my exh=4 | 0.9573 |

Changing exhaustiveness from 4 to 32 perturbs the compound ranking by **+0.0042 Spearman**
relative to simply re-seeding. An eight-fold increase in compute is a *smaller* perturbation
than running the identical protocol again.

This is a ranking statement, independent of the AUROC arithmetic above, and it explains the
result there rather than merely restating it: the intervention is smaller than the noise it was
measured against.

**ρ = 0.9525 across seeds is also a ceiling on any surrogate model** — nothing can agree with
Vina better than Vina agrees with itself. Worth stating in any paper that reports a surrogate's
correlation with docking scores, which conventionally does not.

### 4. The manuscript's own exh=4 run is an ordinary replicate

Its 0.4079 sits **inside** our ten-seed range (0.4074–0.4242), 1.17 upper-bound SDs from the
mean. Earlier in this work I reported it as "4.7σ" and anomalous; that used the **point** SD
at n=5, and it was wrong. Recorded because I wrote the rule against exactly that overstatement
two sections above where I then committed it.

## Two defects in the published pipeline, found while building this

- **No `--seed` at all.** `kaggle_shard_worker.py:46` builds the Vina command without one, so
  every result in the manuscript ran at a seed drawn from system entropy. Seed variance was not
  merely unmeasured — it was already inside both reported "limits".
- **Whole-number affinities discarded.** The score regex requires a decimal, and Vina prints an
  integral affinity as `-7`. Those compounds were recorded as *docking failures*. ~1 in 250;
  unlikely to bias AUROC, but it is silent data loss reported as something else.

## Scope, and what is still owed

Exhaustiveness **4**, not the 32 the headline results use. Stage 1 is an upper bound on stage 2
only if less search is more stochastic — **stage 2 tests that rather than assuming it** and is
running now (6 seeds × 2 shards).

Per-compound instability is far larger than the aggregate suggests: median per-ligand score SD
**0.096 eV**, median rank span **70 places out of 755**, max 549. Individual compound rankings
move a great deal between seeds; AUROC barely does. Anyone ranking compounds off a single Vina
run should know that.

This does **not** address the editor's third point, chemical-series-controlled splits.

Raw: `results/seed_exh4_s*.json` · Analysis: `seed_floor_exh4.json` · Self-test: `selftest.py`
