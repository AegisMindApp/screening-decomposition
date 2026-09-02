# Disagreement-as-signal: predicted 0.45 → 0.60+, measured 0.53 and barely outside control

**2 September 2026.** An external review proposed that ranking only compounds where two
scoring functions agree would lift Mpro "from 0.45 to 0.60+". Tested on 743 compounds with
Vina (repaired receptor) and gnina CNNaffinity. Vina/CNN rank agreement: Spearman r = 0.445.

Filtering changes the population, so every filtered AUROC is compared against **2000 random
subsets of identical size** — the control that can refute the claim. This is the same
selection effect that produced the 149 size-selected dropouts in arm 2.

| keep | n | active % | AUROC vina | AUROC cnn | random-subset 95% | verdict |
|---|---|---|---|---|---|---|
| most-agreeing 25% | 185 | 41.1% | **0.5321** | 0.5313 | [0.3651, 0.5257] | just outside |
| most-agreeing 50% | 371 | 34.0% | 0.4668 | 0.4813 | [0.3991, 0.4907] | inside — no info |
| most-agreeing 75% | 557 | 35.4% | 0.4623 | 0.5023 | [0.4205, 0.4708] | inside — no info |
| all | 743 | 34.6% | 0.4451 | 0.5389 | — | baseline |

**Not confirmed.** The best case reaches **0.5321**, not 0.60+, and only by discarding 75% of
the library. It clears the random-subset ceiling by 0.0064 — a margin smaller than the ~0.04
measurement floor established in `MARGINAL_VALUE.md`, so it is not resolvable on this
benchmark. Note also the active rate rises to 41.1% in the agreeing subset: part of the
movement is composition, not information.

For scale: descriptors alone reach **0.7653** on the same compounds, for free.

## Status of the other proposed directions

| direction | status |
|---|---|
| Residual docking (rank on 3D residual) | **Tested, reverses sign**: Mpro +0.077 [+0.035, +0.119], FXa −0.119 [−0.149, −0.090], both significant. See `MARGINAL_VALUE.md`. |
| "Large residual ⇒ 3D matters" (interpretability) | **Backwards on our data.** The residual helped the *broken* target (Mpro) and hurt the *working* one (FXa). |
| Report ΔAUROC over descriptors | **Already implemented and committed**: −0.0023 Mpro (sig.), +0.0150 FXa (n.s.). |
| Few-shot target-adaptive scoring | Partly pre-empted: our descriptor model *is* a target-specific fit (0.7653 / 0.7129) and adding the 3D term does not improve it. |
| Train on FEP instead of PDBbind | Not testable here — no FEP capability, and FEP's own prospective ranking accuracy is not established as a ceiling worth chasing. |
| **Pose-ensemble scoring** | **The one survivor.** Untested by us; poses are stored top-1 only, so it needs a re-dock saving all modes. CPU-only, no GPU quota needed. |

The PDBbind-trained GNN figure quoted in the review (0.417 [0.372, 0.461]) is ours and is
accurate — `docs/papers/benchmark_paper_draft.md:413`.
