"""Push the Boltz-2 Mpro run. GPU + internet required (weights download, MSA server)."""
import json, sys, tempfile, shutil, urllib.request, urllib.error
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parents[1] / "retrospective_benchmark" / "kaggle"))
from push_mpro import token, _current_version, wait_for_dataset, launcher_src   # noqa: E402

HERE = Path(__file__).parent
HANDLE = "oceansparx/mpro-boltz2"
# The push API rejects a bare kernel name with "Invalid slug"; it must be owner-qualified.
SLUG = "oceansparx/mpro-boltz2-run"

def main():
    tok = token()
    if "--kernels-only" not in sys.argv:
        before = _current_version(HANDLE, tok)
        import kagglehub
        from kagglehub.config import set_kaggle_api_token
        set_kaggle_api_token(tok)
        with tempfile.TemporaryDirectory() as tmp:
            for f in ("mpro_smiles.json", "manifest.json", "mpro_sequence.txt"):
                shutil.copy(HERE / "bundle" / f, Path(tmp) / f)
            kagglehub.dataset_upload(HANDLE, tmp,
                                     version_notes="mpro 751 + dimer seq for boltz2")
        want = (before or 0) + 1
        if not wait_for_dataset(HANDLE, tok, want_version=want):
            raise SystemExit(f"{HANDLE}: version {want} never became attachable - refusing to "
                             f"push a kernel that would read stale data")
    else:
        cur = _current_version(HANDLE, tok)
        if not cur:
            raise SystemExit(f"{HANDLE}: --kernels-only but no published version exists")
        print(f"  reusing {HANDLE} v{cur} (no re-upload)")
    src = launcher_src(HERE / "boltz2_worker.py", worker_args=())
    base = {"slug": SLUG, "newTitle": "mpro boltz2 run", "text": src, "language": "python",
            "kernelType": "script", "isPrivate": True,
            "enableGpu": True, "enableTpu": False, "enableInternet": True,
            "datasetDataSources": [HANDLE], "competitionDataSources": [],
            "kernelDataSources": [], "modelDataSources": [], "categoryIds": []}

    # v3 died on a compute-capability mismatch: modern torch has no kernels below sm_70 and the
    # session was given a P100-class card. Our notes say the push API cannot choose the
    # accelerator, but a rejected push costs no GPU time, so it is free to confirm rather than
    # assume. Try the T4 request first and fall back to the plain body.
    for attempt in ({**base, "acceleratorType": "nvidiaTeslaT4"}, base):
        tag = "with T4 request" if "acceleratorType" in attempt else "plain (card as allocated)"
        req = urllib.request.Request("https://www.kaggle.com/api/v1/kernels/push",
            data=json.dumps(attempt).encode(),
            headers={"Authorization": f"Bearer {tok}", "Content-Type": "application/json"})
        try:
            r = json.load(urllib.request.urlopen(req, timeout=180))
            if r.get("error"):
                print(f"  {tag}: rejected err={r.get('error')!r}")
                continue
            print(f"  kernel {SLUG}: v{r.get('versionNumber')} {tag} "
                  f"bad_datasets={r.get('invalidDatasetSources')}")
            return 0
        except urllib.error.HTTPError as e:
            print(f"  {tag}: HTTP {e.code} {e.read(500).decode()[:250]}")
    return 1

if __name__ == "__main__":
    sys.exit(main())
