#!/usr/bin/env python3
"""Run seed replicates on Kaggle CPU kernels, one seed per kernel.

One seed per kernel deliberately. A kernel that runs out of session time leaves a partial seed,
and a partial seed shrinks the cross-seed intersection for EVERY other seed — so the cheap
failure mode is a kernel that does less, not a chain that half-finishes. Each kernel is
independent and resumable; Kaggle queues whatever it cannot run concurrently.

Reads the bundle from the private dataset `oceansparx/mpro-7vu6-docking-bundle`, copies it to
/kaggle/working so the vina binary can be made executable (/kaggle/input is read-only), and
re-checks the md5s there — the whole experiment rests on the protocol being bit-identical.

    python analysis/seed_floor/kaggle_jobs.py --push --seeds 2 3 4 5
    python analysis/seed_floor/kaggle_jobs.py --status --seeds 2 3 4 5
    python analysis/seed_floor/kaggle_jobs.py --fetch --seeds 2 3 4 5
"""
from __future__ import annotations

import argparse
import base64
import json
from pathlib import Path

import requests

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
USER = "oceansparx"
API = "https://www.kaggle.com/api/v1"
DATASET = f"{USER}/mpro-7vu6-docking-bundle"
OUT = HERE / "results"


def _token() -> str:
    for line in (REPO / ".env").read_text().splitlines():
        if line.startswith("KAGGLE_API_TOKEN="):
            return line.split("=", 1)[1].strip().strip("\"'")
    raise SystemExit("KAGGLE_API_TOKEN not found in .env")


def _hdr():
    return {"Authorization": f"Bearer {_token()}", "Content-Type": "application/json"}


def slug(seed: int, exh: int = 4, shard: int = 0, n_shards: int = 1) -> str:
    base = f"aegismind-seedfloor-exh{exh}-s{seed}"
    return base if n_shards == 1 else f"{base}-sh{shard}of{n_shards}"


TEMPLATE = r'''
import base64, os, shutil, subprocess, sys, glob, zipfile
SEED = __SEED__
EXH = __EXH__
SHARD = __SHARD__
NSHARDS = __NSHARDS__
PANEL = __PANEL__

# Do not assume the mount path. The probe kernel found /kaggle/input containing only
# ['datasets'], not the flat /kaggle/input/<slug> layout, so the bundle is located by SEARCHING
# for its marker directory. Printing what was actually mounted is what made this diagnosable at
# all -- the first probe died on a bare FileNotFoundError that named only the path it wanted.
MARKER = "bundle_CHEMBL4523582_7VU6"
src = None
for root, dirs, files in os.walk("/kaggle/input"):
    if MARKER in dirs:
        src = root
        break
    if os.path.basename(root) == MARKER:
        src = os.path.dirname(root)
        break
if src is None:
    for root, dirs, files in os.walk("/kaggle/input"):
        print("  ", root, dirs[:5], files[:5], flush=True)
    raise SystemExit("bundle not found anywhere under /kaggle/input")
print("bundle root:", src, flush=True)

# NOT /kaggle/working: everything there becomes kernel output, and ~900 bundle
# files then bury the one file that matters behind a paginated listing.
work = "/tmp/bundle"
os.makedirs(work, exist_ok=True)
for d in os.listdir(src):
    sp = os.path.join(src, d)
    if os.path.isdir(sp) and d.startswith("bundle_"):
        shutil.copytree(sp, os.path.join(work, d), dirs_exist_ok=True)
print("bundle contents:", sorted(os.listdir(work))[:6], flush=True)

# /kaggle/input is read-only, so vina cannot be chmod'ed there.
for v in glob.glob(f"{work}/**/vina", recursive=True):
    os.chmod(v, 0o755)

open("/kaggle/working/run_seeds.py", "wb").write(base64.b64decode("__RUNNER_B64__"))
rc = subprocess.run([sys.executable, "/kaggle/working/run_seeds.py",
                     "--exh", str(EXH), "--seeds", str(SEED), "--jobs", "4",
                     "--shard", str(SHARD), "--n-shards", str(NSHARDS),
                     "--panel", str(PANEL),
                     "--bundle", work, "--out-dir", "/kaggle/working"]).returncode
print("runner rc", rc, flush=True)
sys.exit(rc)
'''


def source(seed: int, exh: int, shard: int = 0, n_shards: int = 1, panel: int = 0) -> str:
    s = (TEMPLATE.replace("__SEED__", str(seed)).replace("__EXH__", str(exh))
         .replace("__SHARD__", str(shard)).replace("__NSHARDS__", str(n_shards))
         .replace("__PANEL__", str(panel)))
    s = s.replace("__RUNNER_B64__",
                  base64.b64encode((HERE / "run_seeds.py").read_bytes()).decode())
    assert "__" not in s.replace("__main__", ""), "placeholder left unreplaced"
    return s


def body(seed: int, exh: int, shard: int = 0, n_shards: int = 1, panel: int = 0):
    sl = slug(seed, exh, shard, n_shards)
    return {"slug": f"{USER}/{sl}", "newTitle": sl,
            "text": source(seed, exh, shard, n_shards, panel),
            "language": "python", "kernelType": "script",
            "isPrivate": True, "enableGpu": False, "enableInternet": False,
            "datasetDataSources": [DATASET], "competitionDataSources": [],
            "kernelDataSources": [], "modelDataSources": [], "categoryIds": []}


def verify(a) -> int:
    """Classify every shard by what the status endpoint actually reports.

    The useful distinction is one the raw --status output loses in unparsed error text: 403 means
    the slug was never created, 404 means it exists but has never produced a run. A shard in
    either state contributes nothing to the intersection, and a missing shard does not error --
    it narrows the panel and drops the seed out of the paired estimate. Absence is the dangerous
    case because it reads as nothing.

    WHAT THIS DOES NOT TELL YOU -- corrected 19 Sep 2026, same night, before anything acted on it.
    The first version asserted that PUSHED_NEVER_RAN meant `push_stage2.py` had popped the shard
    after a 200 and forgotten it, and printed "RE-PUSH REQUIRED". That was wrong. Seed 4 shard 1
    was reported in this state at 23:35 and the LIVE PUSHER LAUNCHED IT AT 23:40 -- it had never
    left the pusher's queue. Acting on the verdict would have double-pushed it.

    Of the candidate mechanisms, one is refuted and one is now observed directly.

    Refuted: a capacity-refused push does NOT create the kernel. Seed 4 shard 2 was attempted
    after 23:40, refused at the session cap, and still answered 403.

    Observed: the state is mostly REGISTRATION LAG behind a successful push. Minutes later the
    pusher launched that same seed 4 shard 2 and it immediately read PUSHED_NEVER_RAN -- kernel
    created, run not yet visible to /kernels/status. So the common cause of this state is a shard
    that was just launched perfectly, which is the opposite of lost.

    That leaves no measured case of a genuinely lost shard, and the original reading of seed 4
    shard 1 is recorded as UNESTABLISHED rather than explained away.

    So PUSHED_NEVER_RAN is an observation to reconcile against any running pusher, NOT an
    instruction to re-push. Check `pgrep -af push_stage2` first; this function cannot see it.
    """
    order = ["PUSHED_NEVER_RAN", "NOT_PUSHED", "error", "cancelAcknowledged", "queued",
             "running", "complete"]
    seen = {}
    for s, sh in [(x, y) for x in a.seeds for y in range(a.n_shards)]:
        r = requests.get(f"{API}/kernels/status", headers=_hdr(),
                         params={"userName": USER,
                                 "kernelSlug": slug(s, a.exh, sh, a.n_shards)}, timeout=60)
        if r.status_code == 200:
            st = (r.json() or {}).get("status") or "error"
        elif r.status_code == 404:
            # The kernel exists but has never produced a run. This is the silent gap.
            st = "PUSHED_NEVER_RAN"
        elif r.status_code == 403:
            # Kaggle answers 403, not 404, for a slug that was never created at all.
            st = "NOT_PUSHED"
        else:
            st = "error"
        seen.setdefault(st, []).append((s, sh))
    for st in order + [k for k in seen if k not in order]:
        if st in seen:
            print(f"{st:18} {len(seen[st]):3}  {seen[st]}")
    stuck = seen.get("PUSHED_NEVER_RAN", [])
    if stuck:
        print(f"\nRECONCILE {len(stuck)} shard(s) against any running pusher: {stuck}")
        print("  The kernel exists and has never run. That is LOST only if nothing is still "
              "holding it -- a live pusher retrying at the session cap looks identical from "
              "here. Run `pgrep -af push_stage2` before re-pushing; a shard still in a "
              "pusher's queue will be double-pushed otherwise.")
    return 2 if stuck else 0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", type=int, nargs="+", required=True)
    ap.add_argument("--exh", type=int, default=4)
    # Stage 2 measures ~200 s of CPU per ligand, so one full seed is ~43 core-hours -- more
    # than a 4-core Kaggle session can reliably finish, and a session that runs out leaves a
    # partial seed that shrinks the intersection for every other seed. Two shards per seed is
    # ~5.4 h each, comfortably inside the limit.
    ap.add_argument("--n-shards", type=int, default=1)
    for f in ("push", "status", "fetch", "verify"):
        ap.add_argument(f"--{f}", action="store_true")
    a = ap.parse_args()

    if a.verify:
        return verify(a)

    if a.push:
        # Kaggle returns HTTP 200 for a REJECTED push, putting the reason in errorNullable --
        # e.g. "Maximum batch CPU session count of 5 reached." Reporting the status code alone
        # counts a rejected kernel as launched, and the seed then silently never runs.
        launched, rejected = [], []
        for s, sh in [(x, y) for x in a.seeds for y in range(a.n_shards)]:
            r = requests.post(f"{API}/kernels/push", headers=_hdr(),
                              json=body(s, a.exh, sh, a.n_shards), timeout=300)
            err = ""
            try:
                err = (r.json() or {}).get("error") or (r.json() or {}).get("errorNullable") or ""
            except Exception:
                err = r.text[:120]
            if r.status_code == 200 and not err:
                launched.append((s, sh))
                print(f"seed {s} shard {sh}: launched")
            else:
                rejected.append(((s, sh), err or f"HTTP {r.status_code}"))
                print(f"seed {s} shard {sh}: NOT LAUNCHED — {err or r.status_code}")
        print(f"\n{len(launched)} launched {launched}, {len(rejected)} rejected "
              f"{[x[0] for x in rejected]}")
        return 0 if not rejected else 2
    if a.status:
        for s, sh in [(x, y) for x in a.seeds for y in range(a.n_shards)]:
            r = requests.get(f"{API}/kernels/status", headers=_hdr(),
                             params={"userName": USER,
                                     "kernelSlug": slug(s, a.exh, sh, a.n_shards)}, timeout=60)
            st = r.json().get("status") if r.status_code == 200 else r.text[:60]
            print(f"seed {s} shard {sh}: {st}")
        return 0
    if a.fetch:
        OUT.mkdir(parents=True, exist_ok=True)
        # (s, shard) pairs, not bare seeds: the earlier shard-aware refactor converted the push
        # and status loops but left this one iterating seeds only, so `sh` was unbound and
        # every fetch raised UnboundLocalError.
        for s, sh in [(x, y) for x in a.seeds for y in range(a.n_shards)]:
            # The listing is paginated at 500 files and the kernel copies ~900 bundle files
            # into /kaggle/working, so the result JSON lands on a later page. Fetching only the
            # first page reported "no result file yet" for a run that had finished perfectly.
            files, tok, hit = [], None, None
            want = f"seed_exh{a.exh}_s{s}_shard{sh}of{a.n_shards}.json"
            for _ in range(12):
                params = {"userName": USER,
                          "kernelSlug": slug(s, a.exh, sh, a.n_shards)}
                if tok:
                    params["pageToken"] = tok
                r = requests.get(f"{API}/kernels/output", headers=_hdr(),
                                 params=params, timeout=300)
                if r.status_code != 200:
                    break
                j = r.json()
                page = j.get("files", [])
                files += page
                hit = next((f for f in page if f["fileName"].endswith(want)), None)
                if hit:
                    break
                tok = j.get("nextPageToken")
                if not tok:
                    break
            if hit is None:
                print(f"seed {s}: no result file after {len(files)} files listed")
                continue
            blob = requests.get(hit["url"], timeout=600).content
            d = json.loads(blob)
            # A partial seed is usable but shrinks the intersection for every other seed, so
            # say so loudly rather than letting it quietly narrow the panel.
            tag = "COMPLETE" if d.get("complete") else "PARTIAL"
            (OUT / want).write_bytes(blob)
            print(f"seed {s} shard {sh}: {tag}, {len(d['scores'])} scored, "
                      f"{len(d['failures'])} failed")
        return 0
    ap.print_help()
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
