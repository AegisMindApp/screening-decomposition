#!/usr/bin/env python3
"""The exhaustiveness effect, estimated PAIRED on seed.

Written before stage 2 finished, so the estimator is fixed in advance.

The manuscript estimates the eight-fold search-effort effect from ONE exh=32 run minus ONE
exh=4 run, both unseeded. That difference carries the run-to-run noise of both arms, and
`RESULT.md` shows it does not survive being charged for it.

Stage 1 ran seeds 1-10 at exh=4 and stage 2 runs seeds 1-6 at exh=32 — the same seed numbers —
so the same comparison can be made **within seed**:

    delta_s = AUROC(exh32, seed s) - AUROC(exh4, seed s)

Pairing cancels whatever the two runs share. It is imperfect here: changing exhaustiveness
changes how many Monte-Carlo runs consume the RNG stream, so the streams diverge after the
start rather than staying locked. Whether pairing actually buys anything is therefore something
to MEASURE, not assume — the script reports the paired SD against the unpaired expectation
sqrt(2) x sd(single run), and says which is smaller.

Cross-program comparisons (gnina, Boltz-2, DiffDock) get no benefit from matching seed NUMBERS:
seed 42 in two different programs indexes unrelated streams. What pairs those is reusing the
same poses.

    python analysis/seed_floor/paired_exhaustiveness.py
"""
from __future__ import annotations

import argparse
import glob
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy import stats

HERE = Path(__file__).resolve().parent
from analyse_seeds import auroc  # noqa: E402


def load(exh: int, d: str):
    """seed -> merged scores, only for seeds whose every shard is present and complete."""
    # Key on (seed, n_shards, panel), not seed alone. A superseded run with a different shard
    # count and a different panel sat in the same directory and was merged into the same seed,
    # so four complete 350-panel shards were dragged down by two dead 771-panel partials.
    files = defaultdict(list)
    for f in sorted(glob.glob(f"{d}/seed_exh{exh}_s*.json")):
        j = json.loads(Path(f).read_text())
        files[(j["seed"], j.get("n_shards", 1), j.get("panel"))].append(f)
    out, labels = {}, {}
    for (seed, _nsh, _panel), fs in files.items():
        docs = [json.loads(Path(f).read_text()) for f in fs]
        n_sh = max(x.get("n_shards", 1) for x in docs)
        if any(not x.get("complete") for x in docs) or \
           {x.get("shard", 0) for x in docs} != set(range(n_sh)):
            continue
        m = {}
        for x in docs:
            m.update(x["scores"])
            labels.update(x["labels"])
        out[seed] = m
    return out, labels


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default=str(HERE / "results"))
    a = ap.parse_args()

    e4, lab = load(4, a.dir)
    e32, lab32 = load(32, a.dir)
    lab.update(lab32)
    shared = sorted(set(e4) & set(e32))
    if not shared:
        print(f"no seed has BOTH exh=4 and exh=32 complete yet "
              f"(exh4: {sorted(e4)}, exh32: {sorted(e32)})")
        return 0

    # One compound set for every arm and every seed, or the differences are partly panel.
    common = sorted(set.intersection(*[set(e4[s]) for s in shared],
                                     *[set(e32[s]) for s in shared]))
    y = [lab[c] for c in common]
    print(f"{len(shared)} seeds with both conditions, {len(common)} compounds common to all")

    rows = []
    for s in shared:
        a4 = auroc(y, [e4[s][c] for c in common])
        a32 = auroc(y, [e32[s][c] for c in common])
        rows.append((s, a4, a32, a32 - a4))
        print(f"  seed {s:2d}: exh4 {a4:.4f}  exh32 {a32:.4f}  delta {a32-a4:+.4f}")

    d = np.array([r[3] for r in rows])
    sd4 = float(np.std([r[1] for r in rows], ddof=1)) if len(rows) > 1 else float("nan")
    print()
    print(f"PAIRED delta = {d.mean():+.5f}  sd {d.std(ddof=1):.5f}  "
          f"n={len(d)}  same sign {int((np.sign(d) == np.sign(d.mean())).sum())}/{len(d)}")
    if len(d) > 1:
        t = stats.ttest_1samp(d, 0.0)
        w = stats.wilcoxon(d) if len(d) >= 6 else None
        ci = d.mean() + np.array([-1, 1]) * stats.t.ppf(0.975, len(d) - 1) * d.std(ddof=1) / np.sqrt(len(d))
        print(f"  95% CI [{ci[0]:+.5f}, {ci[1]:+.5f}]   t-test p={t.pvalue:.4f}"
              + (f"   Wilcoxon p={w.pvalue:.4f}" if w else ""))
        print(f"  excludes zero: {ci[0] * ci[1] > 0}")
        # did pairing actually help? unpaired would cost sqrt(2) x the single-run sd
        print()
        print(f"  paired sd of the difference      : {d.std(ddof=1):.5f}")
        print(f"  unpaired expectation sqrt(2)xsd  : {np.sqrt(2) * sd4:.5f}")
        print(f"  -> pairing {'HELPS' if d.std(ddof=1) < np.sqrt(2) * sd4 else 'does NOT help'}"
              f" (measured, not assumed)")
    print()
    print(f"manuscript's single unpaired estimate: +0.0185")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
