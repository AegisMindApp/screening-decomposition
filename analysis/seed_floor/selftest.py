#!/usr/bin/env python3
"""Self-test for the seed-floor analysis, run before the real replicates landed.

Checks the analysis on input whose answer is known, because every number the JCIM repair
rests on comes out of `analyse_seeds.py`. Validating it against real data is impossible --
there is no ground-truth resolution limit to compare against -- so it is validated here
against injected values instead.

    python analysis/seed_floor/selftest.py
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
from scipy import stats

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from analyse_seeds import auroc  # noqa: E402


def main() -> int:
    fails = []

    # 1. Orientation. Vina scores are ENERGIES: more negative binds better. A sign slip
    #    inverts every AUROC silently and nothing downstream would look wrong.
    lab = [1, 1, 1, 0, 0, 0]
    if auroc(lab, [-9, -8, -7, -3, -2, -1]) != 1.0:
        fails.append("perfect ranking did not give AUROC 1.0 -- score sign is inverted")
    if auroc(lab, [-1, -2, -3, -7, -8, -9]) != 0.0:
        fails.append("inverted ranking did not give AUROC 0.0")

    # 2. End to end, with a known per-ligand jitter.
    tmp = tempfile.mkdtemp()
    rng = np.random.default_rng(0)
    labels = {f"C{i}": int(i < 150) for i in range(400)}
    base = {k: (-8.0 if v else -7.0) + rng.normal(0, 1.0) for k, v in labels.items()}
    JIT = 0.05
    for s in range(1, 9):
        r = np.random.default_rng(100 + s)
        json.dump({"exh": 4, "seed": s, "shard": 0, "n_shards": 1, "protocol": {},
                   "labels": labels, "failures": {}, "complete": True,
                   "scores": {k: v + r.normal(0, JIT) for k, v in base.items()}},
                  open(f"{tmp}/seed_exh4_s{s}_shard0of1.json", "w"))
    subprocess.run([sys.executable, str(HERE / "analyse_seeds.py"), "--exh", "4",
                    "--dir", tmp, "--out", f"{tmp}/a.json"], capture_output=True, text=True)
    d = json.load(open(f"{tmp}/a.json"))
    if not 0.8 * JIT < d["per_ligand_score_sd_median"] < 1.2 * JIT:
        fails.append(f"injected jitter {JIT} not recovered "
                     f"({d['per_ligand_score_sd_median']})")

    # 3. The chi-square upper bound IS the reported resolution limit, so it must be right.
    n = d["n_seeds"]
    want = d["auroc_sd"] * np.sqrt((n - 1) / stats.chi2.ppf(0.05, n - 1))
    if not np.isclose(want, d["RESOLUTION_LIMIT_auroc_sd_upper95"]):
        fails.append("95% upper bound does not match the chi-square closed form")

    # 4. Intersection. A compound missing from ONE seed must drop out of ALL of them, or the
    #    per-seed AUROCs cover different panels and the contrast is partly panel composition.
    os.remove(f"{tmp}/seed_exh4_s8_shard0of1.json")
    d7 = json.load(open(f"{tmp}/seed_exh4_s7_shard0of1.json"))
    d7["scores"].pop("C0")
    json.dump(d7, open(f"{tmp}/seed_exh4_s7_shard0of1.json", "w"))
    subprocess.run([sys.executable, str(HERE / "analyse_seeds.py"), "--exh", "4",
                    "--dir", tmp, "--out", f"{tmp}/b.json"], capture_output=True, text=True)
    d2 = json.load(open(f"{tmp}/b.json"))
    if d2["n_compounds_common"] != d["n_compounds_common"] - 1:
        fails.append(f"intersection did not drop the missing compound "
                     f"({d2['n_compounds_common']} vs {d['n_compounds_common'] - 1})")

    for f in fails:
        print("FAIL:", f)
    print("seed-floor analysis self-test:", "PASS" if not fails else f"{len(fails)} FAILURES")
    return 1 if fails else 0


if __name__ == "__main__":
    raise SystemExit(main())


def test_score_regex_accepts_a_whole_number_affinity():
    """Vina prints an integral affinity WITHOUT a decimal point.

    Observed live on CHEMBL4975056: mode 1 printed as "   1           -7", and the original
    pattern -- `(-?\\d+\\.\\d+)`, which is still in kaggle_shard_worker.py:32 and so was used
    for the manuscript's own runs -- discarded it as "rc=0 no-score". Vina had succeeded.
    Silent data loss, dressed as a docking failure.
    """
    import importlib
    import run_seeds
    importlib.reload(run_seeds)
    cases = {
        "   1           -7          0          0": "-7",      # the one that broke it
        "   1       -7.362          0          0": "-7.362",
        "   1        7.5           0          0": "7.5",
        " 1  -10": "-10",
    }
    bad = []
    for line, want in cases.items():
        m = run_seeds.SCORE.search(line + "\n")
        if not m or m.group(1) != want:
            bad.append(f"{line!r} -> {m.group(1) if m else None}, wanted {want}")
    for b in bad:
        print("FAIL:", b)
    assert not bad, "score regex dropped a valid affinity"
