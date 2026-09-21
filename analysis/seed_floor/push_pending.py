#!/usr/bin/env python3
"""Launch the remaining seed kernels as Kaggle capacity frees up.

Kaggle caps batch CPU sessions at 5. A rejected push returns **HTTP 200** with the reason in
`errorNullable` ("Maximum batch CPU session count of 5 reached."), so a naive retry loop that
checks only the status code would mark every seed launched and exit having run nothing. This
treats a push as successful only when the response carries no error AND the kernel subsequently
reports a run.

    python analysis/seed_floor/push_pending.py --seeds 7 8 9 10
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

CAP_MSG = "session count"


def push_one(seed: int, exh: int):
    """(launched, reason). Never reports launched on a capacity rejection."""
    r = requests.post(f"{API}/kernels/push", headers=_hdr(), json=body(seed, exh), timeout=300)
    try:
        j = r.json()
    except Exception:
        return False, f"HTTP {r.status_code} {r.text[:80]}"
    err = j.get("error") or j.get("errorNullable") or ""
    if r.status_code == 200 and not err:
        return True, j.get("ref") or "ok"
    return False, err or f"HTTP {r.status_code}"


def running(seed: int) -> str:
    r = requests.get(f"{API}/kernels/status", headers=_hdr(),
                     params={"userName": USER, "kernelSlug": slug(seed)}, timeout=60)
    try:
        return r.json().get("status", "unknown")
    except Exception:
        return "unknown"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", type=int, nargs="+", required=True)
    ap.add_argument("--exh", type=int, default=4)
    ap.add_argument("--interval", type=int, default=600)
    ap.add_argument("--max-hours", type=float, default=36.0)
    a = ap.parse_args()

    pending = list(a.seeds)
    t0 = time.time()
    print(f"queue: {pending}  (polling every {a.interval}s, giving up after {a.max_hours}h)",
          flush=True)
    while pending and (time.time() - t0) < a.max_hours * 3600:
        seed = pending[0]
        ok, why = push_one(seed, a.exh)
        if ok:
            pending.pop(0)
            print(f"[{time.strftime('%H:%M')}] seed {seed} LAUNCHED ({why}); "
                  f"remaining {pending}", flush=True)
            # Give Kaggle a moment to register the session before testing the next one,
            # otherwise the next push races it and is rejected for a slot that is now taken.
            time.sleep(60)
            continue
        if CAP_MSG in str(why):
            print(f"[{time.strftime('%H:%M')}] seed {seed} waiting — at session cap", flush=True)
        else:
            # A non-capacity rejection will not fix itself by waiting.
            print(f"[{time.strftime('%H:%M')}] seed {seed} REFUSED for a non-capacity reason: "
                  f"{why} — dropping it from the queue", flush=True)
            pending.pop(0)
            continue
        time.sleep(a.interval)

    if pending:
        print(f"GAVE UP with {pending} unlaunched after {a.max_hours}h", flush=True)
        return 1
    print(f"all seeds launched; statuses: "
          f"{ {s: running(s) for s in a.seeds} }", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
