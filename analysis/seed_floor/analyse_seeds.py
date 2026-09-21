#!/usr/bin/env python3
"""Turn the seed replicates into a resolution limit. Implements PREREGISTRATION.md.

    python analysis/seed_floor/analyse_seeds.py --exh 4
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


def auroc(labels, scores):
    """Vina scores are energies: MORE NEGATIVE is a better predicted binder, so the ranking
    variable is the negated score. Getting this backwards silently inverts every AUROC."""
    y = np.asarray(labels, int)
    x = -np.asarray(scores, float)
    pos, neg = x[y == 1], x[y == 0]
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    r = stats.rankdata(np.concatenate([pos, neg]))
    return float((r[:len(pos)].sum() - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg)))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--exh", type=int, required=True)
    ap.add_argument("--dir", default=str(HERE / "results"))
    ap.add_argument("--out", default=None)
    ap.add_argument("--include-partial", action="store_true",
                    help="include seeds still running. Off by default: every seed is scored "
                         "on the INTERSECTION, so one partial seed truncates the panel for "
                         "all of them. Measured live -- adding a seed with 150 docks done "
                         "cut the common panel from ~760 compounds to 149, and AUROC variance "
                         "scales inversely with panel size, so the floor came out inflated.")
    a = ap.parse_args()

    # Group by SEED before judging completeness. Skipping incomplete FILES is wrong once a
    # seed is sharded: a seed whose shard 1 of 3 never finished would still be merged from
    # shards 0 and 2 and treated as complete, so a partial seed silently shrinks the
    # intersection for every other seed. Stage 1 was unsharded so this never bit; stage 2 is.
    files = defaultdict(list)
    for f in sorted(glob.glob(f"{a.dir}/seed_exh{a.exh}_s*.json")):
        files[json.loads(Path(f).read_text())["seed"]].append(f)

    by_seed, labels, incomplete, skipped = defaultdict(dict), {}, [], []
    for seed, fs in sorted(files.items()):
        docs = [json.loads(Path(f).read_text()) for f in fs]
        n_sh = max(d.get("n_shards", 1) for d in docs)
        have = {d.get("shard", 0) for d in docs}
        bad = [Path(f).name for f, d in zip(fs, docs) if not d.get("complete")]
        missing = sorted(set(range(n_sh)) - have)
        if bad or missing:
            incomplete += bad
            reason = []
            if bad:
                reason.append(f"{len(bad)} shard(s) not complete")
            if missing:
                reason.append(f"shard(s) {missing} absent of {n_sh}")
            if not a.include_partial:
                skipped.append({"seed": seed, "reason": "; ".join(reason),
                                "scored_so_far": sum(len(d["scores"]) for d in docs)})
                continue
        for d in docs:
            by_seed[seed].update(d["scores"])
            labels.update(d["labels"])
    if not by_seed:
        raise SystemExit(f"no seed files for exh={a.exh} in {a.dir}")

    seeds = sorted(by_seed)
    # EVERY seed must be scored on the SAME compounds. If seed A times out on a ligand that
    # seed B scores, their AUROCs are computed over different panels and the difference between
    # them is partly panel composition, not seed. That would reintroduce exactly the kind of
    # contamination this measurement exists to remove.
    common = set.intersection(*(set(by_seed[s]) for s in seeds))
    names = sorted(common)
    y = [labels[n] for n in names]
    dropped = {s: len(set(by_seed[s]) - common) for s in seeds}

    aurocs = {s: auroc(y, [by_seed[s][n] for n in names]) for s in seeds}
    vals = np.array([aurocs[s] for s in seeds], float)
    n = len(vals)
    sd = float(vals.std(ddof=1)) if n > 1 else float("nan")
    # 95% one-sided UPPER bound on sigma from chi-square on n-1 df. The upper bound is the
    # resolution limit, not the point estimate: an SD from ~10 replicates carries ~24% relative
    # standard error, and quoting the point estimate as a threshold would overstate the
    # precision of the floor itself -- the same species of error that got the paper declined.
    sd_hi = (float(sd * np.sqrt((n - 1) / stats.chi2.ppf(0.05, n - 1)))
             if n > 1 else float("nan"))

    M = np.array([[by_seed[s][nm] for s in seeds] for nm in names], float)
    per_lig_sd = M.std(axis=1, ddof=1) if n > 1 else np.zeros(len(names))
    ranks = np.array([stats.rankdata(-M[:, i]) for i in range(n)])
    rank_span = (ranks.max(axis=0) - ranks.min(axis=0)) if n > 1 else np.zeros(len(names))

    out = {
        "exhaustiveness": a.exh, "n_seeds": n, "seeds": seeds,
        "n_compounds_common": len(names),
        "n_actives": int(sum(y)), "n_inactives": int(len(y) - sum(y)),
        "dropped_per_seed_vs_common": dropped,
        "incomplete_files": incomplete,
        "skipped_partial_seeds": skipped,
        "auroc_per_seed": {str(s): aurocs[s] for s in seeds},
        "auroc_mean": float(vals.mean()), "auroc_min": float(vals.min()),
        "auroc_max": float(vals.max()), "auroc_range": float(vals.max() - vals.min()),
        "auroc_sd": sd,
        "RESOLUTION_LIMIT_auroc_sd_upper95": sd_hi,
        "per_ligand_score_sd_median": float(np.median(per_lig_sd)),
        "per_ligand_score_sd_max": float(per_lig_sd.max()) if len(per_lig_sd) else None,
        "n_ligands_score_identical_across_seeds": int((per_lig_sd == 0).sum()),
        "rank_span_median": float(np.median(rank_span)),
        "rank_span_max": float(rank_span.max()) if len(rank_span) else None,
    }
    if incomplete:
        out["WARNING"] = (f"{len(incomplete)} shard file(s) not marked complete; this is a "
                          f"partial result and the floor may move")

    dest = Path(a.out or HERE / f"seed_floor_exh{a.exh}.json")
    dest.write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))
    print(f"\n-> {dest}")
    if n > 1:
        print(f"\nRESOLUTION LIMIT (exh={a.exh}): AUROC sd = {sd:.5f}, "
              f"95% upper bound = {sd_hi:.5f}, over {n} seeds on {len(names)} compounds")
        print("Compare against the manuscript's claimed effects; per PREREGISTRATION.md an "
              "effect below this floor does not survive.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
