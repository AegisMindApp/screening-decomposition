#!/usr/bin/env python3
"""Launch stage 2 (exhaustiveness 32) as Kaggle capacity frees up.

Stage 2 is the protocol the manuscript's headline results actually use. Stage 1 (exh=4) is an
upper bound on it only if less search really is more stochastic; this run TESTS that rather
than assuming it.

Sharded two ways per seed. Measured cost is ~200 s of CPU per ligand, so one full seed is
~43 core-hours -- more than a 4-core Kaggle session can reliably finish. A session that runs
out leaves a partial seed, and a partial seed shrinks the cross-seed intersection for every
OTHER seed, so the cheap failure must be a shard that does less rather than a seed that
half-finishes.

    python analysis/seed_floor/push_stage2.py --seeds 1 2 3 4 5 --n-shards 2
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import requests

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from kaggle_jobs import API, USER, _hdr, body, slug  # noqa: E402

CAP = "session count"


def push_one(seed, shard, n_shards, exh, panel=0):
    r = requests.post(f"{API}/kernels/push", headers=_hdr(),
                      json=body(seed, exh, shard, n_shards, panel), timeout=300)
    try:
        j = r.json()
    except Exception:
        return False, f"HTTP {r.status_code}"
    err = j.get("error") or j.get("errorNullable") or ""
    return (r.status_code == 200 and not err), (err or f"HTTP {r.status_code}")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", type=int, nargs="+", required=True)
    ap.add_argument("--n-shards", type=int, default=2)
    ap.add_argument("--exh", type=int, default=32)
    ap.add_argument("--panel", type=int, default=0)
    ap.add_argument("--interval", type=int, default=600)
    ap.add_argument("--max-hours", type=float, default=72.0)
    a = ap.parse_args()

    todo = [(s, sh) for s in a.seeds for sh in range(a.n_shards)]
    print(f"stage 2: exh={a.exh}, {len(todo)} kernels "
          f"({len(a.seeds)} seeds x {a.n_shards} shards)", flush=True)
    t0 = time.time()
    while todo and (time.time() - t0) < a.max_hours * 3600:
        seed, shard = todo[0]
        ok, why = push_one(seed, shard, a.n_shards, a.exh, a.panel)
        if ok:
            todo.pop(0)
            print(f"[{time.strftime('%H:%M')}] exh{a.exh} seed {seed} shard {shard} LAUNCHED; "
                  f"{len(todo)} left", flush=True)
            time.sleep(60)          # let Kaggle register the session before testing the next
            continue
        if CAP in str(why):
            print(f"[{time.strftime('%H:%M')}] waiting — at session cap "
                  f"({len(todo)} kernels left)", flush=True)
        else:
            print(f"[{time.strftime('%H:%M')}] seed {seed} shard {shard} REFUSED, "
                  f"non-capacity: {why} — dropped", flush=True)
            todo.pop(0)
            continue
        time.sleep(a.interval)

    if todo:
        print(f"GAVE UP with {len(todo)} unlaunched: {todo}", flush=True)
        return 1
    print("stage 2 fully launched", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
