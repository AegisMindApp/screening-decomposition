"""Push a Boltz-2 shard to Kaggle. P100-only; the worker's Pascal bootstrap handles the card.

Kaggle caps a session at ~12 h, so long panels are sharded (interleaved, so each shard keeps
the panel's label mix). The launcher sets BOLTZ_* env vars, which the stock launcher_src in
push_mpro.py cannot do.
"""
import argparse, base64, json, sys, tempfile, shutil, urllib.request, urllib.error
from pathlib import Path
HERE = Path(__file__).parent
sys.path.insert(0, str(HERE.parents[1] / "analysis" / "retrospective_benchmark" / "kaggle"))
sys.path.insert(0, str(HERE.parents[0] / "retrospective_benchmark" / "kaggle"))
from push_mpro import token, _current_version, wait_for_dataset          # noqa: E402

def launcher(worker: Path, env: dict) -> str:
    b64 = base64.b64encode(worker.read_bytes()).decode()
    return (
        "import base64, glob, os, subprocess, sys\n"
        f"open('/kaggle/working/_worker.py','wb').write(base64.b64decode('{b64}'))\n"
        + "".join(f"os.environ[{k!r}] = {v!r}\n" for k, v in env.items()) +
        "os.environ['BOLTZ_OUT'] = '/kaggle/working'\n"
        "cands = [d for d in sorted(glob.glob('/kaggle/input/*')) if os.path.isdir(d)]\n"
        "print('input dirs:', cands, flush=True)\n"
        "bundle = next((d for d in cands\n"
        "               if glob.glob(os.path.join(d, '**', 'manifest.json'), recursive=True)),\n"
        "              cands[0] if cands else None)\n"
        "print('bundle ->', bundle, 'env:', {k: v for k, v in os.environ.items()\n"
        "      if k.startswith('BOLTZ_')}, flush=True)\n"
        "if not bundle:\n"
        "    sys.exit('no dataset mounted under /kaggle/input')\n"
        "rc = subprocess.run([sys.executable, '/kaggle/working/_worker.py',\n"
        "                     '--bundle', bundle]).returncode\n"
        "sys.exit(rc)\n")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bundle", required=True, help="local bundle dir")
    ap.add_argument("--handle", required=True, help="kaggle dataset handle owner/name")
    ap.add_argument("--slug", required=True, help="kernel slug, owner-qualified")
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--nshards", type=int, default=1)
    ap.add_argument("--ceiling", type=float, required=True, help="GPU-hour abort ceiling")
    ap.add_argument("--kernels-only", action="store_true")
    a = ap.parse_args()
    tok = token()
    if not a.kernels_only:
        before = _current_version(a.handle, tok)
        import kagglehub
        from kagglehub.config import set_kaggle_api_token
        set_kaggle_api_token(tok)
        with tempfile.TemporaryDirectory() as tmp:
            for f in Path(a.bundle).iterdir():
                if f.is_file(): shutil.copy(f, Path(tmp) / f.name)
            kagglehub.dataset_upload(a.handle, tmp, version_notes=f"boltz2 {Path(a.bundle).name}")
        want = (before or 0) + 1
        if not wait_for_dataset(a.handle, tok, want_version=want):
            raise SystemExit(f"{a.handle}: v{want} never became attachable")
    env = {"BOLTZ_CEILING_HOURS": str(a.ceiling), "BOLTZ_POCKET": "1"}
    if a.nshards > 1:
        env["BOLTZ_NSHARDS"] = str(a.nshards); env["BOLTZ_SHARD"] = str(a.shard)
    body = {"slug": a.slug, "newTitle": a.slug.split("/")[-1].replace("-", " "),
            "text": launcher(HERE / "boltz2_worker.py", env), "language": "python",
            "kernelType": "script", "isPrivate": True,
            "enableGpu": True, "enableTpu": False, "enableInternet": True,
            "datasetDataSources": [a.handle], "competitionDataSources": [],
            "kernelDataSources": [], "modelDataSources": [], "categoryIds": []}
    req = urllib.request.Request("https://www.kaggle.com/api/v1/kernels/push",
        data=json.dumps(body).encode(),
        headers={"Authorization": f"Bearer {tok}", "Content-Type": "application/json"})
    try:
        r = json.load(urllib.request.urlopen(req, timeout=180))
        print(f"  {a.slug}: v{r.get('versionNumber')} err={r.get('error')!r} "
              f"shard={a.shard}/{a.nshards} ceiling={a.ceiling}h")
        return 0 if not r.get("error") else 1
    except urllib.error.HTTPError as e:
        print(f"  {a.slug}: HTTP {e.code} {e.read(400).decode()[:250]}"); return 1

if __name__ == "__main__":
    sys.exit(main())
