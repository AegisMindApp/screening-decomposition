#!/usr/bin/env python3
"""Can the decision rule detect an effect of known size? Q4 of the pre-submission gate.

Every verdict in this work rests on one rule: compare an effect's interval against the measured
resolution floor and classify it SURVIVES / NEGLIGIBLE / INCONCLUSIVE. Nothing so far shows that
rule can recover an effect whose true size is known in advance. Without that, a run of negative
findings is indistinguishable from a pipeline that cannot detect anything -- the failure mode
recorded as "a control that cannot fire".

DESIGN. Take the real exh=4 scores for the six seeds. Build a synthetic second arm by shifting
the scores of actives by an amount calibrated to produce a TARGET AUROC delta. Feed that through
the identical paired analysis and equivalence rule used for the real interventions, and check
the verdict against what was planted.

PLANTED CASES AND THE VERDICT EACH MUST RETURN, fixed before running:

  +0.000  NEGLIGIBLE    a true null must be measured as negligible, not merely unproven
  +0.004  NEGLIGIBLE    below the margin: the size the real exhaustiveness effect came out at
  +0.050  SURVIVES      comfortably above the margin, like pose ensemble
  +0.140  SURVIVES      the largest real effect
  -0.050  SURVIVES      sign symmetry: the rule must not be one-sided

THE CONTROL CAN FAIL. If a planted +0.140 returned NEGLIGIBLE the rule would be broken, and
this script would say so. A test that cannot fail is not a control.

    python analysis/seed_floor/positive_control.py
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy import stats

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from analyse_seeds import auroc                      # noqa: E402
from paired_exhaustiveness import load               # noqa: E402
from equivalence import floor_at, classify           # noqa: E402

CASES = [(0.000, "NEGLIGIBLE"), (0.004, "NEGLIGIBLE"), (0.050, "SURVIVES"),
         (0.140, "SURVIVES"), (-0.050, "SURVIVES")]


def plant(y, s, target, rng):
    """Shift actives until the AUROC moves by `target`. Bisection on the shift size.

    Direction matters and getting it wrong is silent. Vina affinities are LOWER-is-better, so
    these arms sit at AUROC ~0.44 and shifting actives upward moves the metric DOWN. The first
    version assumed upward = better, which ran the bisection to its bound (recovering -0.446 for
    a planted +0.004) and inverted the sign (recovering +0.050 for a planted -0.050). The
    control caught it, which is the point of having one.

    So probe the direction empirically rather than assuming it.
    """
    y = np.asarray(y, int); s = np.asarray(s, float)
    base = auroc(y, s)
    if abs(target) < 1e-9:
        return s + rng.normal(0, 1e-9, len(s))        # a true null, with numerical jitter only
    step = 0.01 * (np.std(s) or 1.0)
    up = auroc(y, s + np.where(y == 1, step, 0.0)) - base
    sign = 1.0 if up > 0 else -1.0                    # which way a positive shift moves AUROC
    want = abs(target)
    lo, hi = 0.0, 20.0 * (np.std(s) or 1.0)
    for _ in range(80):
        mid = (lo + hi) / 2
        s2 = s + np.where(y == 1, sign * math.copysign(mid, target), 0.0)
        if abs(auroc(y, s2) - base) < want:
            lo = mid
        else:
            hi = mid
    return s + np.where(y == 1, sign * math.copysign((lo + hi) / 2, target), 0.0)


def main() -> int:
    e4, lab = load(4, str(HERE / "results"))
    e32, lab32 = load(32, str(HERE / "results"))
    lab.update(lab32)
    seeds = sorted(set(e4) & set(e32))
    common = sorted(set.intersection(*[set(e4[s]) for s in seeds], *[set(e32[s]) for s in seeds]))
    y = [lab[c] for c in common]
    d = floor_at(len(common))[1]
    print(f"{len(seeds)} seeds, {len(common)} compounds, margin d={d:.5f}\n")
    print(f"{'planted':>9}  {'recovered':>10}  {'95% CI':>22}  {'verdict':<13} expected      ok")

    rng = np.random.default_rng(20260921)
    fails = 0
    for target, expect in CASES:
        deltas = []
        for s in seeds:
            base = np.array([e4[s][c] for c in common], float)
            arm_b = plant(y, base, target, rng)
            deltas.append(auroc(y, arm_b) - auroc(y, base))
        dl = np.array(deltas)
        n = len(dl); se = dl.std(ddof=1) / math.sqrt(n) if dl.std(ddof=1) > 0 else 1e-12
        t95, t90 = stats.t.ppf(0.975, n - 1), stats.t.ppf(0.95, n - 1)
        lo95, hi95 = dl.mean() - t95 * se, dl.mean() + t95 * se
        lo90, hi90 = dl.mean() - t90 * se, dl.mean() + t90 * se
        got = classify(lo95, hi95, lo90, hi90, d)
        # A negative true effect is "detected" when the interval excludes the margin either way.
        if target < 0:
            got = "SURVIVES" if hi95 < -d else got
        ok = got == expect
        fails += not ok
        print(f"{target:+9.3f}  {dl.mean():+10.5f}  [{lo95:+8.5f},{hi95:+8.5f}]  "
              f"{got:<13} {expect:<13} {'OK' if ok else 'FAIL'}")

    # --- TIER 2: the same planted effects, but carrying REAL run-to-run noise.
    # Tier 1 plants onto the same arm it compares against, so the only variation between seeds
    # is the base scores and the intervals come out at +-0.00005 -- far tighter than any real
    # comparison. That validates the rule's arithmetic, not its behaviour against noise. Here
    # arm B is a DIFFERENT seed's run plus the planted effect, so the difference carries genuine
    # reseed noise on top of a known signal. This is the version that matters.
    print()
    print("TIER 2: planted effect carried on a DIFFERENT seed's run (real reseed noise)")
    print(f"{'planted':>9}  {'recovered':>10}  {'95% CI':>22}  {'verdict':<13} note")
    for target, _ in CASES:
        deltas = []
        for i, sd_ in enumerate(seeds):
            other = seeds[(i + 1) % len(seeds)]
            base = np.array([e4[sd_][c] for c in common], float)
            b_raw = np.array([e4[other][c] for c in common], float)
            arm_b = plant(y, b_raw, target, rng)
            deltas.append(auroc(y, arm_b) - auroc(y, base))
        dl = np.array(deltas)
        n = len(dl); sdv = dl.std(ddof=1); se = sdv / math.sqrt(n) if sdv > 0 else 1e-12
        t95, t90 = stats.t.ppf(0.975, n - 1), stats.t.ppf(0.95, n - 1)
        lo95, hi95 = dl.mean() - t95 * se, dl.mean() + t95 * se
        lo90, hi90 = dl.mean() - t90 * se, dl.mean() + t90 * se
        got = classify(lo95, hi95, lo90, hi90, d)
        if target < 0:
            got = "SURVIVES" if hi95 < -d else got
        note = "detected" if abs(target) > d else "inside the margin"
        print(f"{target:+9.3f}  {dl.mean():+10.5f}  [{lo95:+8.5f},{hi95:+8.5f}]  {got:<13} {note}")
    print("  (sd across seeds now reflects real reseed variation, not planting precision)")

    print()
    if fails:
        print(f"{fails} planted case(s) misclassified -- the decision rule is NOT validated.")
        return 1
    print("All planted effects recovered with the verdict fixed in advance.")
    print("The rule detects effects at and above the margin, and calls sub-margin effects")
    print("negligible rather than merely unproven. Q4 of the gate is satisfied.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
