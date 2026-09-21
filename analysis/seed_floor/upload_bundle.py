#!/usr/bin/env python3
"""Put the Mpro docking bundle on Kaggle as a PRIVATE dataset, so seed replicates can run
across Kaggle CPU kernels instead of only the local box.

Everything uploaded is public-origin: ChEMBL-derived ligand pdbqt files, the PDB 7VU6 receptor,
and a public AutoDock Vina binary. The dataset is created private regardless.

The vina binary and receptor are what the run_seeds md5 assertions check, so the zip must
preserve them byte-for-byte — verified here before upload and again by the runner on Kaggle.

    python analysis/seed_floor/upload_bundle.py --zip      # build + verify the archive
    python analysis/seed_floor/upload_bundle.py --upload   # create the private dataset
    python analysis/seed_floor/upload_bundle.py --status
"""
from __future__ import annotations

import argparse
import hashlib
import json
import zipfile
from pathlib import Path

import requests

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
SRC = REPO / "analysis" / "retrospective_benchmark" / "kaggle"
DIRS = ["bundle_CHEMBL4523582_7VU6/shard0", "bundle_CHEMBL4523582_7VU6/shard1",
        "bundle_CHEMBL4523582_7VU6_overflow/shard0"]
ZIP = HERE / "mpro_bundle.zip"
USER = "oceansparx"
SLUG = "mpro-7vu6-docking-bundle"
API = "https://www.kaggle.com/api/v1"


def _token() -> str:
    for line in (REPO / ".env").read_text().splitlines():
        if line.startswith("KAGGLE_API_TOKEN="):
            return line.split("=", 1)[1].strip().strip("\"'")
    raise SystemExit("KAGGLE_API_TOKEN not found in .env")


def _hdr(ct=True):
    h = {"Authorization": f"Bearer {_token()}"}
    if ct:
        h["Content-Type"] = "application/json"
    return h


def md5(p: Path) -> str:
    return hashlib.md5(p.read_bytes()).hexdigest()


def build_zip() -> int:
    n = 0
    with zipfile.ZipFile(ZIP, "w", zipfile.ZIP_DEFLATED) as z:
        for d in DIRS:
            root = SRC / d
            if not root.exists():
                raise SystemExit(f"missing {root}")
            for f in sorted(root.rglob("*")):
                if f.is_file() and "__pycache__" not in str(f):
                    z.write(f, f"{d}/{f.relative_to(root)}")
                    n += 1
    # The md5 assertions are the point of the experiment: verify the round trip NOW, not on
    # Kaggle after an hour of docking.
    with zipfile.ZipFile(ZIP) as z:
        for d in DIRS:
            for name in ("vina", "receptor.pdbqt"):
                disk = md5(SRC / d / name)
                inzip = hashlib.md5(z.read(f"{d}/{name}")).hexdigest()
                assert disk == inzip, f"{d}/{name} changed in the zip ({disk} vs {inzip})"
    print(f"{ZIP.name}: {n} files, {ZIP.stat().st_size/1e6:.1f} MB, "
          f"vina+receptor md5 verified byte-identical through the archive")
    return 0


def upload() -> int:
    if not ZIP.exists():
        raise SystemExit("run --zip first")
    blob = ZIP.read_bytes()
    # fileName must be FORM-encoded, not JSON. Sent as JSON the server accepts the request and
    # returns a valid token, but the GCS object is created as the placeholder "filename.ext" --
    # so the upload succeeds, the PUT returns 200, and only the later create call fails, with
    # the unhelpful "Path must be non-null".
    r = requests.post(f"{API}/datasets/upload/file/{len(blob)}/0",
                      headers=_hdr(ct=False), data={"fileName": ZIP.name}, timeout=120)
    if r.status_code != 200:
        print("request-upload failed:", r.status_code, r.text[:400])
        return 1
    j = r.json()
    tok, url = j.get("token"), j.get("createUrl")
    if not url:
        print("no createUrl in response:", json.dumps(j)[:400])
        return 1
    up = requests.put(url, data=blob, timeout=1800)
    print("put:", up.status_code)
    if up.status_code not in (200, 201, 204):
        print(up.text[:300])
        return 1
    meta = {
        "title": "Mpro 7VU6 docking bundle",
        "slug": SLUG,
        "ownerSlug": USER,
        "licenseName": "other",
        "isPrivate": True,
        "subtitleNullable": "AutoDock Vina bundle for the Mpro retrospective panel",
        "descriptionNullable": (
            "Receptor (PDB 7VU6), AutoDock Vina binary and ChEMBL-derived ligand pdbqt files "
            "for the SARS-CoV-2 Mpro retrospective screening panel. Used to measure the "
            "run-to-run resolution limit of the docking protocol by replicating an unchanged "
            "protocol across independent --seed values. All contents are public-origin."),
        # "Path must be non-null": the create endpoint wants the filename on each file entry,
        # not just the upload token the PUT returned.
        "files": [{"token": tok, "path": ZIP.name, "description": "bundle archive"}],
        "isTitleUnique": False,
    }
    c = requests.post(f"{API}/datasets/create/new", headers=_hdr(), json=meta, timeout=600)
    print("create:", c.status_code, c.text[:400])
    return 0 if c.status_code == 200 else 1


def status() -> int:
    r = requests.get(f"{API}/datasets/list", headers=_hdr(),
                     params={"user": USER, "pageSize": 50}, timeout=60)
    rows = r.json() if r.status_code == 200 else []
    print(f"{len(rows)} dataset(s)")
    for d in rows:
        print(" ", d.get("ref"), "|", d.get("title"), "|", d.get("totalBytes"))
    return 0


def main() -> int:
    ap = argparse.ArgumentParser()
    for f in ("zip", "upload", "status"):
        ap.add_argument(f"--{f}", action="store_true")
    a = ap.parse_args()
    if a.zip:
        return build_zip()
    if a.upload:
        return upload()
    if a.status:
        return status()
    ap.print_help()
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
