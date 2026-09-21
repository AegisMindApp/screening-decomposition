# Choosing a clean third target: ALDH1

Written 19 Sep 2026, before the panel is built or any Boltz-2 compute is spent.

## Why a third target at all

The harness refuses a transferability verdict below three admissible targets, because two give
one transfer direction each way and no replication. Mpro and Factor Xa are admissible. PD-L1
is **not** — descriptor baseline 0.874, 79% active, inactives 127 Da heavier than actives — and
was recorded as compromised in `analysis/boltz2/PREREGISTRATION_TARGETS.md` in September, before
it ran. It was nonetheless used in the first harness run, where it carried the largest residual
margin and both apparently-usable cold transfers. Removing it is what left us at two.

## Selection, measured rather than argued

Admissibility was re-checked **on the panel we would actually build** — 300 actives + 450
inactives, 40% active — because the published LIT-PCBA baselines are computed on the full
imbalanced sets (25,000 inactives) and balance changes AUROC.

| target | descriptor baseline on the built panel | class |
|---|---|---|
| **ALDH1** | **0.5698** | oxidoreductase |
| MTORC1 | 0.5922 | kinase complex |
| KAT2A | 0.6656 | acetyltransferase |
| VDR | 0.6854 | nuclear receptor |
| PKM2 | 0.7244 | kinase |
| MAPK1 | 0.7346 | kinase |

Three reasons for ALDH1:

1. **Cleanest baseline by a clear margin — 0.5698.** For comparison the two admissible targets
   we already have score 0.7654 (Mpro) and 0.7041 (Factor Xa) on descriptors. Lower is better
   here: it is the headroom in which a method can demonstrate it exceeds properties at all.
2. **A different protein class.** Mpro is a viral cysteine protease and Factor Xa a serine
   protease; two proteases is a narrow basis for claiming transferability. An oxidoreductase
   makes the transfer genuinely hard, and a failure would be informative rather than ambiguous.
3. **7,168 actives available**, so a 300-active subsample is not scraping the barrel — unlike
   IDH1 (39), ADRB2 (17) or ESR1_ago (13).

**Rejected:** PPARG (0.8459) and OPRK1 (0.8640) fail admissibility outright — properties already
separate those panels.

## Pre-registered before the run

- **Panel**: 300 actives + 450 inactives from `analysis/litpcba/full/ALDH1/`, drawn with
  `random.Random(20260919)`, fixed before selection. Compounds that RDKit cannot parse are
  skipped and counted.
- **Admissibility is re-checked on the final panel** and recorded. If the built panel's
  descriptor baseline exceeds 0.80 the target is dropped and this document is amended — the
  check is not a formality it has already passed.
- **Boltz-2 cost**: 101 s/ligand measured on L4 (`boltz2_probe_L4.json`), so ~21 GPU hours for
  750 compounds. Sharded three ways at ~7 h each, inside the session limit with margin. That
  estimate is from a *measured* probe on this exact model, not extrapolated from a convenient
  subset — the error that cost two stage-2 kernels.
- **Nothing is claimed from ALDH1 alone.** Its purpose is to take the harness from two
  admissible targets to three so a transferability verdict becomes possible at all.

## What this does not fix

A third target makes a verdict *possible*, not *strong*. Three targets give six transfer
directions with no replication within a pair. If Boltz-2 comes back `DEPLOYABLE` on three, that
is a starting position for pointing it at objects, not a finished validation.
