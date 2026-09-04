# Boltz-2 is as property-driven as Vina — but its residual carries far more signal

**4 September 2026.** Computed from arm 1's existing scores at zero compute cost, which is why
`affinity_mw_correction` was deliberately left off arm 2: the question it addresses is
answerable without a run.

| | corr with MW | descriptors explain | raw AUROC | **residual AUROC** |
|---|---|---|---|---|
| Vina (repaired receptor) | 0.562 | 62.8% | 0.4527 | 0.5382 |
| **Boltz-2** `prob_binary` | 0.357 | 63.1% | **0.7913** | **0.6418** |

Residual AUROC = what survives after the seven descriptors are regressed out (out-of-fold,
gradient-boosted).

## The finding

**Descriptors explain essentially the same fraction of both scores — 63% each.** Boltz-2 is
not less property-driven than docking in aggregate; its correlation with molecular weight is
lower (0.357 vs 0.562) but the total descriptor-explainable share is identical.

**What differs is the residual.** Vina's target-specific remainder scores 0.5382. Boltz-2's
scores **0.6418**. For context, every docking score measured this week — Vina raw, Vina
repaired, gnina CNNaffinity, gnina CNNscore, across Mpro and Factor Xa — produced a residual
in the band **0.507–0.555**. Boltz-2 sits clearly outside it.

That is the mechanism behind arm 1's +0.0934 in combination with descriptors. Boltz-2 does not
add value by being a better molecular-property calculator; it adds value because the part of
its prediction that is *not* explainable by cheap descriptors actually ranks actives.

## Why this matters for the earlier conclusion

The session's docking work concluded that the target-specific content of a docking score is
worth ~0.51–0.55 AUROC, and that a derivative docking algorithm was not worth building because
the prize was ~0.011 AUROC over free descriptors. That conclusion stands **for docking**.

This result shows the ceiling was a property of the *scoring approach*, not of the task. A
co-folding affinity model raises the target-specific component to 0.64 on the same compounds,
same labels, same folds — while running blind, against a comparator that was given the binding
site.

## Caveat

One target, one configuration, arm 1's blind run. Mpro is rigid and well-represented in
training data. Arm 2 (pocket-conditioned) is running and will show whether supplying the site
raises the residual further.
