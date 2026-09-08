# Leakage check: is SARS-CoV-2 Mpro in FlashBind's training data?

**8 September 2026.** Required by `PREREGISTRATION.md` (`c6b4504f9`) before any FlashBind score
is interpreted. A model trained on the target it is being tested on would clear the residual
ceiling for the wrong reason.

## What the binary checkpoint is trained on

`predict_binary.sh` runs against `./data/mf-pcba/`. The dataset repo
(`clorf6/FlashBind`, MIT) ships `mf-pcba`, `lit-pcba`, `pubchem`, `antibiotic`, `ESIBank`,
`casp16`, `fep4`, `openfe`, and uses SAIR for affinity regression.

## MF-PCBA: no Mpro, no protease at all

MF-PCBA is 60 curated PubChem assay pairs. The AIDs are recoverable from the filenames in
`davidbuterez/mf-pcba/retrieve-scripts` (`{DR_AID}-{SD_AID}.sh`): **105 distinct AIDs**, all 105
resolved to assay names through the PubChem PUG REST summary endpoint.

| pattern | matches |
|---|---|
| `\bsars\b`, `\bcov-?2\b`, `\b3cl\b`, `\bmpro\b`, `main protease`, `\bcoronavirus\b`, `\bcovid\b` | **0** |
| `\bprotease\b` | **0** |

Not one assay in MF-PCBA targets a protease of any kind, let alone Mpro.

## LIT-PCBA: no Mpro either

Its 15 targets are ADRB2, ALDH1, ESR1_ago, ESR1_ant, FEN1, GBA, IDH1, KAT2A, MAPK1, MTORC1,
OPRK1, PKM2, PPARG, TP53, VDR — verified against our own per-target results in
`analysis/litpcba/descriptor_baseline_full.json`. Mpro is absent.

## Position

**Leakage from the two named binary-task sources is ruled out for the target.** Neither contains
Mpro, and MF-PCBA contains no protease at all, so the binary head has not been fit to this
protein class.

**Not ruled out:** compound-level overlap. `pubchem.tar.zst` is broad and our ChEMBL compounds
could appear in it under other assays. That is a weaker concern — a model that has seen a
compound in an unrelated assay has not learned this target's structure–activity relationship —
but it is unresolved, and any result is reported with that stated.

## A caution about how this check was run

The first pass matched substrings and reported **3 apparent SARS-CoV-2 hits**. All three were
artefacts: "cov" inside "Dis**cov**ery" and "Non-**Cov**alent". Word-boundary matching gives
zero. This is the third substring collision of the day — the others dropped every tetracycline
determinant by matching "tetr" inside "tetracycline", and inflated literature-scan scores by
counting nested phrases three times. Substring matching on scientific text produces confident
wrong answers in both directions.
