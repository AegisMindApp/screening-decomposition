#!/usr/bin/env python3
"""A formal equivalence framework for the five interventions — the editor's second remedy.

The manuscript's decision rule is binary: an effect either clears the resolution floor or it
does not. That conflates two very different situations, and the conflation is exactly what a
reviewer should object to:

  * an effect MEASURED to be negligible — its interval sits inside the resolution, so we can
    say the intervention does nothing worth having; and
  * an effect we simply could not resolve — its interval is wider than the resolution, so the
    data are silent.

Reporting both as "does not survive" claims knowledge in the second case that we do not have.

THE FRAMEWORK. With the measured floor as the equivalence margin d (the natural choice: it is
the smallest difference this protocol can resolve), each effect falls into exactly one of three
states:

  SURVIVES     lower bound of the 95% CI  >  +d        the effect exceeds the resolution
  NEGLIGIBLE   the 90% CI lies inside (-d, +d)         TOST: equivalent to zero within resolution
  INCONCLUSIVE otherwise                                underpowered; the data do not decide

The 90% CI is the standard TOST construction for two one-sided tests at alpha = 0.05; using the
95% CI there would be the common error that makes equivalence harder to declare than it should
be.

    python analysis/seed_floor/equivalence.py
"""
from __future__ import annotations

import glob
import json
import math
from pathlib import Path

import numpy as np
from scipy import stats

HERE = Path(__file__).resolve().parent
PANEL_REF = 755          # panel the manuscript's interventions are computed on


def floor_at(n_compounds: int) -> tuple[float, float]:
    """exh=32 floor, rescaled to a panel of n compounds.

    Floors are panel-size specific: AUROC variance goes as 1/n, verified to within 9% in
    floor_exh4_vs_exh32.py. Charging a floor measured on 343 compounds against a claim computed
    on 755 without rescaling would overstate it by sqrt(755/343) = 1.48x.
    """
    lad = json.loads((HERE.parent / "method_bench/method_ladder.json").read_text())
    f = lad["noise_floors"]["vina_exh32"]
    k = math.sqrt(f["panel_n_compounds"] / n_compounds)
    return f["auroc_sd"] * k, f["auroc_floor_95upper"] * k


def classify(lo95, hi95, lo90, hi90, d):
    if lo95 > d:
        return "SURVIVES"
    if lo90 > -d and hi90 < d:
        return "NEGLIGIBLE"
    return "INCONCLUSIVE"


def main() -> int:
    sd_f, d = floor_at(PANEL_REF)
    print(f"equivalence margin d = {d:.5f}  (exh=32 floor, 95% upper bound, rescaled to "
          f"n={PANEL_REF}; point sd {sd_f:.5f})\n")

    # --- the exhaustiveness effect, from the six paired seeds we actually measured
    import sys
    sys.path.insert(0, str(HERE))
    from paired_exhaustiveness import load
    from analyse_seeds import auroc
    e4, lab = load(4, str(HERE / "results"))
    e32, lab32 = load(32, str(HERE / "results"))
    lab.update(lab32)
    seeds = sorted(set(e4) & set(e32))
    common = sorted(set.intersection(*[set(e4[s]) for s in seeds], *[set(e32[s]) for s in seeds]))
    y = [lab[c] for c in common]
    dl = np.array([auroc(y, [e32[s][c] for c in common]) - auroc(y, [e4[s][c] for c in common])
                   for s in seeds])
    n = len(dl); se = dl.std(ddof=1) / math.sqrt(n)
    d_local = floor_at(len(common))[1]
    t95 = stats.t.ppf(0.975, n - 1); t90 = stats.t.ppf(0.95, n - 1)
    lo95, hi95 = dl.mean() - t95 * se, dl.mean() + t95 * se
    lo90, hi90 = dl.mean() - t90 * se, dl.mean() + t90 * se
    # TOST explicitly, so the p-values are on record rather than implied by the interval.
    p_lo = stats.t.sf((dl.mean() + d_local) / se, n - 1)
    p_hi = stats.t.cdf((dl.mean() - d_local) / se, n - 1)
    print(f"eight-fold search effort, from {n} paired seeds on {len(common)} compounds")
    print(f"  delta {dl.mean():+.5f}   95% CI [{lo95:+.5f}, {hi95:+.5f}]")
    print(f"                           90% CI [{lo90:+.5f}, {hi90:+.5f}]   d={d_local:.5f}")
    print(f"  TOST p (vs -d) {p_lo:.4f}   TOST p (vs +d) {p_hi:.4f}   "
          f"equivalent: {max(p_lo, p_hi) < 0.05}")
    print(f"  => {classify(lo95, hi95, lo90, hi90, d_local)}\n")

    # --- the other four, from their stored intervals
    lad = json.loads((HERE.parent / "method_bench/method_ladder.json").read_text())
    print(f"the remaining interventions, against d={d:.5f} at n={PANEL_REF}")
    for k, v in lad["interventions"].items():
        if k.startswith("_") or k == "eightfold_search_effort":
            continue
        lo, hi = v["ci_with_seed"]
        half90 = (hi - lo) / 2 * (stats.norm.ppf(0.95) / stats.norm.ppf(0.975))
        mid = (hi + lo) / 2
        print(f"  {k:32} {v['delta']:+.4f}  95% [{lo:+.4f},{hi:+.4f}]  "
              f"-> {classify(lo, hi, mid - half90, mid + half90, d)}")
    print()
    print("Reading: NEGLIGIBLE is a positive finding -- the intervention is measured to do")
    print("nothing worth having. INCONCLUSIVE is an admission that the data do not decide.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
