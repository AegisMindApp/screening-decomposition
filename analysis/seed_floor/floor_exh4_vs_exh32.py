#!/usr/bin/env python3
"""Is the exh=4 reproducibility floor an upper bound on the exh=32 floor?

`analysis/seed_floor/RESULT.md` measured the floor at exhaustiveness 4 and the manuscript's
headline results are at 32. Carrying the exh=4 number across assumes less search is MORE
stochastic, so that the cheap measurement bounds the expensive one. That assumption was never
tested -- `docs/papers/RESUBMISSION_GATE.md` blocks on exactly this substitution -- and the raw
comparison is confounded, because the two stages ran different panels (755 vs 343 compounds)
and different seed counts, and AUROC variance scales inversely with panel size.

So compare MATCHED: same compounds, same seeds, only exhaustiveness differing.

The variance test has to respect the design. The two arms share seeds, so their AUROCs are
correlated and an F-test (which assumes independent samples) is the wrong instrument. The
Pitman-Morgan test is the paired equivalent: for paired (x, y) the variances are equal exactly
when (x-y) and (x+y) are uncorrelated, so it tests that correlation.

    python analysis/seed_floor/floor_exh4_vs_exh32.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy import stats

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from analyse_seeds import auroc                        # noqa: E402
from paired_exhaustiveness import load                 # noqa: E402


def chi2_upper(sd: float, n: int, conf: float = 0.95) -> float:
    """One-sided upper confidence bound on a standard deviation from n observations.

    A point sd from a handful of seeds is itself a small-sample estimate. Reporting it as the
    floor understates the floor roughly half the time, so the reported figure is the bound.
    """
    return float(sd * np.sqrt((n - 1) / stats.chi2.ppf(1 - conf, n - 1)))


def main() -> int:
    e4, lab = load(4, str(HERE / "results"))
    e32, lab32 = load(32, str(HERE / "results"))
    lab.update(lab32)
    shared = sorted(set(e4) & set(e32))
    if len(shared) < 3:
        print(f"need >=3 seeds with both conditions; have {shared}")
        return 1
    common = sorted(set.intersection(*[set(e4[s]) for s in shared],
                                     *[set(e32[s]) for s in shared]))
    y = [lab[c] for c in common]
    a4 = np.array([auroc(y, [e4[s][c] for c in common]) for s in shared])
    a32 = np.array([auroc(y, [e32[s][c] for c in common]) for s in shared])
    n = len(shared)

    print(f"MATCHED on {len(common)} compounds and seeds {shared}\n")
    sd4, sd32 = a4.std(ddof=1), a32.std(ddof=1)
    for nm, a, sd in (("exh=4 ", a4, sd4), ("exh=32", a32, sd32)):
        print(f"  {nm}  mean {a.mean():.4f}  sd {sd:.5f}  "
              f"95% upper bound {chi2_upper(sd, n):.5f}")

    print(f"\n  ratio sd(exh32)/sd(exh4) = {sd32 / sd4:.2f}")

    # Pitman-Morgan: paired test of equal variances.
    r, p = stats.pearsonr(a32 - a4, a32 + a4)
    print(f"  Pitman-Morgan paired variance test: r={r:+.3f}  p={p:.4f}  (n={n})")
    if p < 0.05:
        verdict = ("the two floors DIFFER" if sd32 > sd4 else
                   "exh=32 is genuinely quieter")
    else:
        verdict = "the two floors are NOT distinguishable at this n"
    print(f"  -> {verdict}")

    # An independent check that the matched comparison is measuring what it claims.
    # Stage 1 measured this same exh=4 arm on 755 compounds across 10 seeds. If AUROC variance
    # really scales as 1/n_compounds, that number predicts the sd we just measured on a smaller
    # panel -- a prediction made from data this script did not use.
    SD_STAGE1, N_STAGE1 = 0.00497, 755
    pred = SD_STAGE1 * np.sqrt(N_STAGE1 / len(common))
    print(f"\n  Panel-size check (prediction from stage 1, not a fit):")
    print(f"    exh=4 sd was {SD_STAGE1:.5f} on {N_STAGE1} compounds; scaled to {len(common)} "
          f"that predicts {pred:.5f}")
    print(f"    measured here: {sd4:.5f}   ratio {sd4 / pred:.3f}")
    print(f"    => the floor is PANEL-SIZE SPECIFIC. A floor measured on {N_STAGE1} compounds is")
    print(f"       roughly {np.sqrt(N_STAGE1 / len(common)):.2f}x too small for a claim computed "
          f"on {len(common)}.")

    print("\n  What this licenses:")
    bigger = sd32 > sd4
    print(f"    exh=4 sd is {'SMALLER' if bigger else 'LARGER'} than exh=32 sd, so substituting "
          f"the exh=4 floor for exh=32\n    is {'NOT conservative' if bigger else 'conservative'}"
          f" at the point estimate.")
    print("    Either way the honest move is to use the exh=32 floor that is now measured,")
    print("    rather than argue about whether a substitution happens to be safe.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
