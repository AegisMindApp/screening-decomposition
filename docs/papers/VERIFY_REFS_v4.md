# Reference and claim-provenance check, `resolution_v4_draft.md` — 25 September 2026

## Reference metadata: 8 of 8 VERIFIED

| # | First author | Year | Venue | DOI | Status |
|---|---|---|---|---|---|
| 1 | Trott | 2010 | J Comput Chem 31:455–461 | 10.1002/jcc.21334 | VERIFIED |
| 2 | Eberhardt | 2021 | J Chem Inf Model 61:3891–3898 | 10.1021/acs.jcim.1c00203 | VERIFIED |
| 3 | McNutt | 2021 | J Cheminform 13:43 | 10.1186/s13321-021-00522-2 | VERIFIED |
| 4 | Butina | 1999 | J Chem Inf Comput Sci 39:747–750 | 10.1021/ci9803381 | VERIFIED |
| 5 | Lakens | 2017 | Soc Psychol Personal Sci 8:355–362 | 10.1177/1948550617697177 | VERIFIED |
| 6 | Tran-Nguyen | 2020 | J Chem Inf Model 60:4263–4273 | 10.1021/acs.jcim.0c00155 | VERIFIED |
| 7 | Sieg | 2019 | J Chem Inf Model 59:947–961 | 10.1021/acs.jcim.8b00712 | VERIFIED |
| 8 | Passaro | 2025 | bioRxiv | 10.1101/2025.06.14.659707 | VERIFIED |

Ref 1 returns 2009 from CrossRef — that is the online-first date; the issue is 2010, volume and
pages as cited. No correction needed.

**The draft carried no DOIs; all eight are now added.** That matters more than it sounds, because
a bibliographic search alone would have put five wrong objects in the reference list.

### The trap, hit again

Searching CrossRef by title — the obvious approach — returned, as top hit:

| ref | what the search returned | what it should be |
|---|---|---|
| 2 | ChemRxiv preprint `10.26434/chemrxiv.14774223.v1` | the JCIM article |
| 3 | ChemRxiv preprint `10.26434/chemrxiv.13578140.v1` | the J Cheminform article |
| 5 | PsyArXiv preprint `10.31234/osf.io/97gpc_v1` | the SPPS article |
| 6 | **supplementary-information record** `…0c00155.s002` | the article `…0c00155` |
| 7 | **supplementary-information record** `…8b00712.s001` | the article `…8b00712` |

Five of eight. Each looks plausible, each resolves, and a `.s002` suffix is easy to miss in a
reference list. The same failure was recorded against v3; it recurs because bibliographic search
ranks by string similarity and a preprint or an SI record matches the title as well as the
article does. **Resolve the canonical DOI directly and compare what comes back.**

## Claim provenance

| Ref | Claim in manuscript | Verdict |
|---|---|---|
| 1,2 | "AutoDock Vina [1,2], box centred on the co-crystallised ligand" — tool attribution only | SUPPORTED |
| 3 | "scoring function (gnina CNN [3]) on identical poses" — tool attribution | SUPPORTED |
| 4 | "Butina clustering [4] on Morgan fingerprints … Tanimoto 0.65" — method attribution | SUPPORTED |
| 5 | "the standard two-one-sided-tests construction at α = 0.05 [5]" | SUPPORTED |
| 6 | ~~"Panel construction follows the LIT-PCBA convention [6]"~~ | **UNSUPPORTED — corrected** |
| 7 | "threshold-defined actives carry a known risk that measured enrichment reflects property matching rather than molecular recognition [7]" | SUPPORTED |
| 8 | "correct binding site to a co-folding model [8]" | SUPPORTED |

### The one defect, and how it got there

v3 read "Panel construction follows the LIT-PCBA convention [6] **for the third target**" — ALDH1,
which genuinely is a LIT-PCBA target. v4 dropped ALDH1 from scope, and when the v4 audit found
refs 6 and 7 orphaned (listed but never cited), I re-attached them to §2.1 — which describes
**Mpro and Factor Xa**, neither of which is a LIT-PCBA target. Both are ChEMBL
pChEMBL-threshold panels.

So a fix for a trivial defect created a **false methodological claim about another group's
dataset**: that our panels were built to a standard they were not built to, and carry debiasing
they do not carry. It is the more serious of the two errors by a wide margin.

Corrected to state the opposite, which is both true and more useful to a reader:

> Both panels are built from ChEMBL activity data by pChEMBL threshold — **not** from a curated
> unbiased benchmark such as LIT-PCBA [6], and we do not claim the debiasing that set applies.

The Sieg citation [7] is retained and its characterisation checked: the paper's claim is about
benchmark bias letting a model exploit non-recognition features, which is that work's subject.
§3.5's descriptor comparison (0.714 against Vina's 0.416) is the paper's own evidence that the
risk is realised on this panel, and that is now said outright.

## Gate

```
References:       8 cited | 8 listed | 0 orphans | 8 VERIFIED with canonical DOIs
Claim provenance: 7 SUPPORTED | 1 UNSUPPORTED — FIXED | 0 CONTRADICTED
Quotations:       none in the manuscript
```

**verify-refs: PASS.** With the v4 audit's four fixes, no submission gate remains open.
