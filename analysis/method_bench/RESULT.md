# First harness run — and a target that should never have counted

**CORRECTED 19 Sep 2026, same day.** The first version of this file reported Vina as
`NOT_BEYOND_PROPERTIES` **across 3 targets**. One of those targets, PD-L1, is inadmissible: its
descriptor baseline is 0.874, so seven free physicochemical properties already separate the
panel and nothing measured on it can show a method exceeds properties. Our own
`analysis/boltz2/PREREGISTRATION_TARGETS.md` recorded it as compromised in September, before it
ran — *"no strong conclusion should be drawn from PD-L1 in either direction"* — and I used it
anyway.

With PD-L1 removed both methods have **2 admissible targets** and the harness returns
`INSUFFICIENT` for each, refusing a transferability verdict. The harness now tests target
admissibility itself rather than relying on someone remembering the caveat.

**What survives the correction:** the per-target residuals, which are the substance. On both
admissible targets Vina sits at or below a measured pure-descriptor null — Factor Xa −0.0128,
Mpro +0.0171, against a 0.00817 floor. That, plus the manuscript's residual table and the
independent AcrB decoy control, is what the generation gate rests on. It does **not** rest on a
transferability verdict, because we do not have one.

Run 19 Sep 2026 on data already on disk — no new compute. `harness.py`, self-tested by
`selftest.py` before being trusted.

## Vina — 2 admissible targets (PD-L1 excluded) — `INSUFFICIENT`

| target | in-target | residual | measured null | margin |
|---|---|---|---|---|
| Factor Xa | 0.6645 | 0.5283 | 0.5411 | **−0.0128** |
| Mpro | 0.4268 | 0.5222 | 0.5051 | +0.0171 |

PD-L1 is excluded (descriptor baseline 0.874). For the record it scored residual 0.5740 against
a null of 0.4726 — the largest margin of the three, and the reason excluding it matters: it was
carrying the verdict.

Worst margin **−0.013**, against a floor of 0.008. On Factor Xa — the target where Vina's
in-target AUROC looks best at 0.664 — its residual is *below* what a pure-descriptor score
achieves. **The in-target number is properties.**

Cold transfer is unusable in 4 of 6 directions (0.14 to 0.69), and the two that clear 0.60 are
both *into* PD-L1 or *from* it, which is the panel most likely to be idiosyncratic (see caveat).

This formalises what the decoy control found for AcrB and what the manuscript found for the
Mpro residuals. It is now a reproducible verdict from a harness rather than a one-off analysis.

## Boltz-2 — 2 targets — `INSUFFICIENT`

The harness refuses to judge transferability on two targets, which is correct: two targets give
one direction each way and no replication. What the two do show:

| target | in-target | residual | measured null | margin |
|---|---|---|---|---|
| Mpro | 0.7910 | 0.6406 | 0.5835 | **+0.0571** |
| Factor Xa | 0.7227 | 0.5474 | 0.5431 | **+0.0042** |

Against Boltz-2's own floor of 0.0201: **Mpro clears, Factor Xa does not.** So even setting
transferability aside, "beyond properties" holds on one target of two. That is a sharper
statement than the manuscript's "replicates in direction but weakly", and it is the same
conclusion reached by a rule fixed in advance rather than by reading the numbers.

Cold transfer: FXa→Mpro **0.680** (usable), Mpro→FXa **0.497** (chance). One direction of two.

## What is missing, and what it costs

**Boltz-2 on PD-L1.** That is the single job standing between here and a transferability verdict
for the method we would most like to deploy. It is a GPU run over 460 compounds, queued behind
the gnina determinism test.

Everything else is already in place: panels, labels, SMILES, descriptors, and a floor per
method.

## Caveats

- **PD-L1's panel is unusual** — 334 actives against 84 inactives (80% active), and its
  molecular-weight AUROC is **0.123**, i.e. strongly *anti*-correlated with activity. Verdicts
  involving PD-L1 should be read with that in mind; it is the reason two cold transfers look
  usable when the other four do not.
- **Vina's floor (0.00817) is measured at exhaustiveness 4** on Mpro. Applying it to Factor Xa
  and PD-L1 assumes it carries across targets, which is unverified — the same substitution the
  resubmission gate blocks on for exh=32. Flagged rather than fixed.
- Boltz-2's floor of 0.0201 comes from `TRANSFER.json` and was not re-derived here.

## Why the numbers can be trusted at all

`selftest.py` plants four methods with known answers and caught a real defect on first run: a
**pure-descriptor** method scored residual 0.548 and passed a 0.5-based "beyond properties"
check, because out-of-fold residualisation cannot fully remove descriptor signal. The null is
now measured per target — the residual of the method's own out-of-fold descriptor prediction —
and the planted case is correctly flagged. Without that fix, **Vina would have passed as
"beyond properties" on all three targets.**
