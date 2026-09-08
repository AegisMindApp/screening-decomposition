"""Push the FlashBind staged probe to Kaggle GPU."""
import json, sys, tempfile, shutil, urllib.request, urllib.error
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parents[1] / "retrospective_benchmark" / "kaggle"))
from push_mpro import token, _current_version, wait_for_dataset, launcher_src  # noqa: E402

HERE = Path(__file__).parent
HANDLE = "oceansparx/flashbind-mpro"
SLUG = "oceansparx/flashbind-probe"

def main():
    tok = token()
    before = _current_version(HANDLE, tok)
    import kagglehub
    from kagglehub.config import set_kaggle_api_token
    set_kaggle_api_token(tok)
    with tempfile.TemporaryDirectory() as tmp:
        for f in HERE.glob("bundle/*"):
            shutil.copy(f, Path(tmp) / f.name)
        kagglehub.dataset_upload(HANDLE, tmp, version_notes="mpro 751 + probe subset for FlashBind")
    want = (before or 0) + 1
    if not wait_for_dataset(HANDLE, tok, want_version=want):
        raise SystemExit(f"{HANDLE}: v{want} never became attachable")
    full = "--full" in sys.argv
    text = launcher_src(HERE / "flashbind_worker.py")
    if full:
        # The launcher runs the worker as a subprocess, so the flag has to be in the
        # environment it inherits, not in argv the worker never reads.
        text = "import os\nos.environ['FLASHBIND_FULL'] = '1'\n" + text
    body = {"slug": SLUG, "newTitle": "flashbind probe", "text": text,
            "language": "python", "kernelType": "script", "isPrivate": True,
            "enableGpu": True, "enableTpu": False, "enableInternet": True,
            "datasetDataSources": [HANDLE], "competitionDataSources": [],
            "kernelDataSources": [], "modelDataSources": [], "categoryIds": []}
    req = urllib.request.Request("https://www.kaggle.com/api/v1/kernels/push",
        data=json.dumps(body).encode(),
        headers={"Authorization": f"Bearer {tok}", "Content-Type": "application/json"})
    r = json.load(urllib.request.urlopen(req, timeout=180))
    print(f"  kernel {SLUG}: v{r.get('versionNumber')} err={r.get('error')!r}")
    return 0 if not r.get("error") else 1

if __name__ == "__main__":
    sys.exit(main())
