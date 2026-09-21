#!/usr/bin/env python3
"""Run Boltz-2 on ALDH1 — the clean third target — sharded across Kaggle GPU.

Reuses analysis/boltz2/boltz2_worker.py unchanged so the third target goes through the
identical code path as Mpro and Factor Xa. Sharding is via the worker's own BOLTZ_SHARD /
BOLTZ_NSHARDS environment variables.

    python analysis/method_bench/run_aldh1_boltz2.py --zip --upload
    python analysis/method_bench/run_aldh1_boltz2.py --push
    python analysis/method_bench/run_aldh1_boltz2.py --status --fetch
"""
from __future__ import annotations

import argparse
import base64
import json
import zipfile
from pathlib import Path

import requests

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
BUNDLE = HERE / "bundle_aldh1"
WORKER = REPO / "analysis" / "boltz2" / "boltz2_worker.py"
ZIP = HERE / "aldh1_bundle.zip"
USER, SLUG = "oceansparx", "aldh1-boltz2-bundle"
API = "https://www.kaggle.com/api/v1"
NSHARDS = 4
# 750 compounds / 4 shards x 101 s/ligand measured = ~5.3 h per shard. The worker's default
# self-abort is 6.0 h, which would trip on a slightly slow shard; 8.0 leaves margin and still
# aborts well inside a GPU session.
CEILING_HOURS = 8.0


def _token():
    for line in (REPO / ".env").read_text().splitlines():
        if line.startswith("KAGGLE_API_TOKEN="):
            return line.split("=", 1)[1].strip().strip("\"'")
    raise SystemExit("KAGGLE_API_TOKEN not in .env")


def _hdr(ct=True):
    h = {"Authorization": f"Bearer {_token()}"}
    if ct:
        h["Content-Type"] = "application/json"
    return h


def slug(sh):
    return f"aegismind-boltz2-aldh1-sh{sh}of{NSHARDS}"


SRC = r'''
import base64, os, sys, glob
os.environ["BOLTZ_SHARD"] = "__SHARD__"
os.environ["BOLTZ_NSHARDS"] = "__NSHARDS__"
os.environ["BOLTZ_CEILING_HOURS"] = "__CEIL__"
os.environ["BOLTZ_OUT"] = "/kaggle/working"

# Do not assume the mount path; find the bundle by a file it must contain.
src = None
for root, dirs, files in os.walk("/kaggle/input"):
    if "sequence.txt" in files and "pocket.json" in files:
        src = root; break
if src is None:
    for root, dirs, files in os.walk("/kaggle/input"):
        print("  ", root, dirs[:4], files[:4], flush=True)
    raise SystemExit("ALDH1 bundle not found under /kaggle/input")
print("bundle:", src, sorted(os.listdir(src)), flush=True)

open("/kaggle/working/boltz2_worker.py", "wb").write(base64.b64decode("__WORKER_B64__"))
sys.argv = ["boltz2_worker.py", "--bundle", src]
exec(open("/kaggle/working/boltz2_worker.py").read())
'''


def main() -> int:
    ap = argparse.ArgumentParser()
    for f in ("zip", "upload", "push", "status", "fetch"):
        ap.add_argument(f"--{f}", action="store_true")
    a = ap.parse_args()

    if a.zip:
        with zipfile.ZipFile(ZIP, "w", zipfile.ZIP_DEFLATED) as z:
            for f in sorted(BUNDLE.iterdir()):
                z.write(f, f.name)
        print(f"{ZIP.name}: {ZIP.stat().st_size/1e6:.2f} MB, {len(list(BUNDLE.iterdir()))} files")

    if a.upload:
        blob = ZIP.read_bytes()
        r = requests.post(f"{API}/datasets/upload/file/{len(blob)}/0",
                          headers=_hdr(ct=False), data={"fileName": ZIP.name}, timeout=120)
        j = r.json()
        requests.put(j["createUrl"], data=blob, timeout=1800)
        meta = {"title": "ALDH1 Boltz-2 bundle", "slug": SLUG, "ownerSlug": USER,
                "licenseName": "other", "isPrivate": True,
                "subtitleNullable": "ALDH1 panel, sequence and pocket for Boltz-2",
                "descriptionNullable": "LIT-PCBA ALDH1: 750 compounds (300 active), 4wp7 "
                                       "sequence and pocket. Public-origin.",
                "files": [{"token": j["token"], "path": ZIP.name}], "isTitleUnique": False}
        c = requests.post(f"{API}/datasets/create/new", headers=_hdr(), json=meta, timeout=600)
        print("create:", c.status_code, c.text[:160])

    if a.push:
        wb = base64.b64encode(WORKER.read_bytes()).decode()
        launched, rejected = [], []
        for sh in range(NSHARDS):
            text = (SRC.replace("__SHARD__", str(sh)).replace("__NSHARDS__", str(NSHARDS))
                    .replace("__CEIL__", str(CEILING_HOURS)).replace("__WORKER_B64__", wb))
            body = {"slug": f"{USER}/{slug(sh)}", "newTitle": slug(sh), "text": text,
                    "language": "python", "kernelType": "script", "isPrivate": True,
                    "enableGpu": True, "enableInternet": True,
                    "datasetDataSources": [f"{USER}/{SLUG}"], "competitionDataSources": [],
                    "kernelDataSources": [], "modelDataSources": [], "categoryIds": []}
            r = requests.post(f"{API}/kernels/push", headers=_hdr(), json=body, timeout=300)
            err = (r.json() or {}).get("error") or (r.json() or {}).get("errorNullable") or ""
            (launched if (r.status_code == 200 and not err) else rejected).append((sh, err))
            print(f"  shard {sh}: {'launched' if not err else 'NOT LAUNCHED — ' + err}")
        print(f"{len(launched)} launched, {len(rejected)} rejected")

    if a.status:
        for sh in range(NSHARDS):
            r = requests.get(f"{API}/kernels/status", headers=_hdr(),
                             params={"userName": USER, "kernelSlug": slug(sh)}, timeout=60)
            print(f"  shard {sh}: {r.json().get('status') if r.status_code == 200 else r.text[:50]}")

    if a.fetch:
        out = HERE / "scores"
        out.mkdir(exist_ok=True)
        for sh in range(NSHARDS):
            r = requests.get(f"{API}/kernels/output", headers=_hdr(),
                             params={"userName": USER, "kernelSlug": slug(sh)}, timeout=180)
            files = r.json().get("files", []) if r.status_code == 200 else []
            hit = next((f for f in files if "boltz2_scores" in f["fileName"]), None)
            if not hit:
                print(f"  shard {sh}: no scores yet")
                continue
            blob = requests.get(hit["url"], timeout=300).content
            (out / f"boltz2_aldh1_s{sh}.json").write_bytes(blob)
            print(f"  shard {sh}: {len(json.loads(blob))} scored")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
