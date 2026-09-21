#!/usr/bin/env python3
"""Can the harness tell a good method from a bad one? Three planted methods with known answers.

A harness nobody has tested against a known answer is decoration. Each case below has an
answer fixed by construction, and the middle one is the case that matters -- a method that is
PURE DESCRIPTORS must come out at chance on the residual, or the residualisation is broken and
every 'beyond properties' claim the harness makes is worthless.

    python analysis/method_bench/selftest.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from harness import evaluate  # noqa: E402


def make(rng, n=400, n_desc=7):
    D = rng.normal(size=(n, n_desc))
    # Noise 1.8, not 0.6: at 0.6 the descriptor baseline is 0.914 and the target is
    # INADMISSIBLE under the harness's own rule -- the planted panels would be as compromised
    # as PD-L1 and the self-test would exercise the wrong code path.
    y = (D[:, 0] + rng.normal(0, 1.8, n) > 0).astype(int)   # activity driven by descriptor 0
    return D, y


def build(kind, rng, n=400):
    D, y = make(rng, n)
    if kind == "oracle":
        score = y + rng.normal(0, 0.01, n)          # near-perfect
    elif kind == "descriptors_only":
        score = D[:, 0] * 2.0 + D[:, 1] * 0.5       # a pure function of descriptors
    elif kind == "random":
        score = rng.normal(size=n)
    elif kind == "genuine":
        score = y * 1.5 + D[:, 0] * 0.3 + rng.normal(0, 0.8, n)  # signal beyond properties
    return {"y": y.tolist(), "score": np.asarray(score).tolist(), "D": D.tolist()}


def main() -> int:
    rng = np.random.default_rng(0)
    fails = []
    for kind, want in (("oracle", "high"), ("descriptors_only", "chance_residual"),
                       ("random", "chance"), ("genuine", "high")):
        targets = {f"T{i}": build(kind, rng) for i in range(3)}
        r = evaluate(kind, targets, floor=0.02)
        it = np.mean(list(r["in_target"].values()))
        res = np.mean(list(r["residual"].values()))
        print(f"{kind:18} in-target {it:.3f}  residual {res:.3f}  -> {r['verdict'][:52]}")
        if want == "high" and it < 0.75:
            fails.append(f"{kind}: in-target {it:.3f} should be high")
        if want == "chance" and abs(it - 0.5) > 0.08:
            fails.append(f"{kind}: in-target {it:.3f} should be ~0.5")
        # THE load-bearing check: a pure-descriptor method must not look 'beyond properties'
        if want == "chance_residual":
            if abs(res - 0.5) > 0.10:
                fails.append(f"descriptors_only: residual {res:.3f} should be ~0.5 -- "
                             f"residualisation is not removing descriptor signal")
            if "NOT_BEYOND_PROPERTIES" not in r["verdict"]:
                fails.append("descriptors_only was NOT flagged NOT_BEYOND_PROPERTIES")
    for f in fails:
        print("FAIL:", f)
    print("\nharness self-test:", "PASS" if not fails else f"{len(fails)} FAILURES")
    return 1 if fails else 0


if __name__ == "__main__":
    raise SystemExit(main())
