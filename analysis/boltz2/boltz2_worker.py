"""Boltz-2 on the Mpro benchmark. Reading rules: analysis/boltz2/PREREGISTRATION.md (9e3083d02).

Stage 1 is a timing probe. It runs 4 ligands, then 8 more, and takes the SECOND batch's rate as
the marginal per-ligand cost -- model load and MSA are one-time and would otherwise inflate the
projection and trigger a false abort. If the projected full-panel cost exceeds the
pre-registered 6.0 GPU-hour ceiling, the run is abandoned and reported as not attempted.

Every per-compound failure is recorded with its message. Nothing is swallowed: a compound that
produces no score appears in the failures file. A bare try/except that hid a parse error is what
made an earlier analysis silently report zero scored compounds.
"""
import json, os, glob, time, subprocess, sys, shutil
from pathlib import Path

BUNDLE = Path(sys.argv[sys.argv.index("--bundle") + 1]) if "--bundle" in sys.argv else Path(".")
OUT = Path(os.environ.get("BOLTZ_OUT", "/kaggle/working"))
# Pre-registered abort threshold. 6.0 on Kaggle (weekly quota); raised to 20.0 on paid GCP
# where the binding constraint is money -- see PREREGISTRATION.md Amendment 2, committed
# before any GCP run. No scientific reading rule changed.
CEILING_HOURS = float(os.environ.get("BOLTZ_CEILING_HOURS", "6.0"))
PROBE_A, PROBE_B = 4, 8
_TAG = ("_s" + os.environ["BOLTZ_SHARD"]) if os.environ.get("BOLTZ_NSHARDS","1") != "1" else ""

def sh(cmd, **kw):
    return subprocess.run(cmd, shell=isinstance(cmd, str), capture_output=True, text=True, **kw)

# Check the card FIRST. v3 spent 14 minutes installing boltz, downloading CCD data and both
# weight files, and generating four MSAs, only to die on:
#     CUDA error: no kernel image is available for execution on the device
# That is a compute-capability mismatch -- modern torch wheels ship no kernels for sm_60
# (P100). Kaggle's push API cannot select the accelerator, so the card is whatever the session
# was given, and the only fixes are a torch downgrade or changing the accelerator in the UI.
# Detecting it here costs seconds instead of a quarter-hour of a sub-8-hour weekly quota.
print("== accelerator", flush=True)
_smi = sh("nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader").stdout
print(_smi.strip() or "(nvidia-smi produced nothing)", flush=True)

# Kaggle hands out P100 (sm_60) and the push API cannot change that: acceleratorType and
# machineShape are both accepted and silently ignored, and manual UI selection has never
# produced a T4 either. Modern torch (2.10+cu128) ships no kernels below sm_70, which is how v3
# died after 14 minutes of install, weight download and MSA generation.
#
# The repo already solved this for the scaling-law work: install a cu121 torch that still
# carries sm_60 kernels, BEFORE boltz. boltz 2.2.1 requires only torch>=2.2, so 2.4.1 satisfies
# it and pip will not upgrade back to a build that cannot run this card.
PASCAL = any(g in _smi for g in ("P100", "P40", "K80"))
# torchvision must be downgraded WITH torch. v5 installed torch 2.4.1 alone and boltz then died
# on "operator torchvision::nms does not exist" -- the preinstalled torchvision was compiled
# against torch 2.10 and its C++ ops do not register on 2.4.1. torchmetrics imports torchvision
# at module level (functional.image.arniqa), so pytorch-lightning cannot even load without it.
# 2.4.1 pairs with torchvision 0.19.1.
PIN = ["torch==2.4.1", "torchvision==0.19.1"]
IDX = ["--index-url", "https://download.pytorch.org/whl/cu121"]
if PASCAL:
    print(f"  Pascal-or-older; installing {' '.join(PIN)} (cu121, ships sm_60)", flush=True)
    _r = sh([sys.executable, "-m", "pip", "install", "-q", *PIN, *IDX], timeout=3600)
    if _r.returncode:
        sys.exit(f"pascal torch bootstrap failed rc={_r.returncode}:\n{(_r.stderr or '')[-2000:]}")
    print("  torch/torchvision pinned", flush=True)

print("== install", flush=True)
# Pin the torch/torchvision pair that is ALREADY WORKING before installing boltz.
#
# boltz[cuda] upgraded the deep-learning image's torch 2.9.1+cu129 to 2.14.0+cu130, which
# orphaned the preinstalled torchvision (built against 2.9.1) and killed the run with
# "RuntimeError: operator torchvision::nms does not exist" -- torchmetrics imports torchvision
# at module level, so pytorch-lightning cannot load and boltz never starts. Same symptom we hit
# on Kaggle for the opposite reason (there we downgraded torch and left torchvision behind).
#
# boltz only requires torch>=2.2, so constraining it to whatever is installed satisfies the
# dependency and stops pip moving a stack that already works. On Pascal this runs after the
# 2.4.1/0.19.1 pin above, so the constraint simply locks that pair in instead.
_pv = sh([sys.executable, "-c",
          "import torch;print(torch.__version__.split('+')[0]);"
          "import torchvision;print(torchvision.__version__.split('+')[0])"])
_cons = None
if _pv.returncode == 0 and len(_pv.stdout.strip().splitlines()) == 2:
    _t, _tv = _pv.stdout.strip().splitlines()
    _cons = OUT.parent / "boltz_constraints.txt"
    _cons.parent.mkdir(parents=True, exist_ok=True)
    _cons.write_text(f"torch=={_t}\ntorchvision=={_tv}\n")
    print(f"  constraining torch=={_t} torchvision=={_tv}", flush=True)

# boltz dispatches to cuequivariance kernels on sm_70+ cards and imports cuequivariance_torch
# at PREDICTION time, well after MSA generation. Those packages live in the [cuda] extra.
# Omitting them cost a 518 s L4 run that died on ModuleNotFoundError after four MSA server
# calls. Pascal never takes that path, which is why the Kaggle P100 run reached prediction
# without them -- the first card that actually works is the one that needs the extra.
_base = [sys.executable, "-m", "pip", "install", "-q"]
_cflag = ["-c", str(_cons)] if _cons else []
r = sh([*_base, *_cflag, "boltz[cuda]"])
if r.returncode:
    print(f"  boltz[cuda] constrained failed; retrying unconstrained\n{(r.stderr or '')[-600:]}",
          flush=True)
    r = sh([*_base, "boltz[cuda]"])
    if r.returncode:
        r = sh([*_base, "boltz"])
        if r.returncode:
            sys.exit("pip install boltz failed:\n" + (r.stderr or "")[-3000:])

# Whatever route we took, torchvision must import and register its native ops, or
# pytorch-lightning dies on load. Repair by pairing rule (torchvision minor = torch minor + 15
# for torch 2.x) rather than failing outright.
_chk = sh([sys.executable, "-c", "import torchvision.ops, torchvision; print('tv ok')"])
if _chk.returncode:
    _tv2 = sh([sys.executable, "-c", "import torch;print(torch.__version__.split('+')[0])"])
    _tt = (_tv2.stdout or "").strip()
    try:
        _maj, _min = _tt.split(".")[:2]
        _want = f"0.{int(_min) + 15}.0" if _maj == "2" else None
    except Exception:
        _want = None
    if _want:
        print(f"  torchvision broken against torch {_tt}; installing torchvision=={_want}",
              flush=True)
        sh([*_base, f"torchvision=={_want}"])
        _chk = sh([sys.executable, "-c", "import torchvision.ops; print('tv ok')"])
    if _chk.returncode:
        sys.exit(f"torchvision unusable against torch {_tt}:\n{(_chk.stderr or '')[-2000:]}")
print("  torch/torchvision pair OK", flush=True)

# boltz's resolver may have pulled either package forward again. Re-pin unconditionally on
# Pascal: wheels are cached so this is cheap, and it guarantees the pairing survives whatever
# dependency resolution did.
if PASCAL:
    _r = sh([sys.executable, "-m", "pip", "install", "-q", *PIN, *IDX], timeout=3600)
    if _r.returncode:
        sys.exit(f"re-pin after boltz failed rc={_r.returncode}:\n{(_r.stderr or '')[-2000:]}")
    _tv = sh([sys.executable, "-c",
              "import torch,torchvision;print(torch.__version__,torchvision.__version__);"
              "import torchvision.ops as o;print('nms ok')"])
    if _tv.returncode:
        sys.exit(f"torch/torchvision still mismatched after re-pin:\n{(_tv.stderr or '')[-2500:]}")
    print(f"  re-pinned: {_tv.stdout.strip().replace(chr(10), '  ')}", flush=True)

# Verify AFTER installing boltz: if its dependency resolution pulled torch forward again, the
# card is unusable and we must find out now, not 14 minutes deep. get_arch_list() is the
# definitive check -- it reports what this build actually compiled kernels for.
_v = sh([sys.executable, "-c",
         "import torch;print(torch.__version__);print(torch.cuda.is_available());"
         "print(','.join(torch.cuda.get_arch_list()));"
         "print('%d%d'%torch.cuda.get_device_capability(0))"])
if _v.returncode:
    sys.exit(f"torch import failed after install:\n{(_v.stderr or '')[-2000:]}")
_ver, _avail, _arch, _cap = (_v.stdout.strip().split("\n") + ["", "", "", ""])[:4]
print(f"  torch {_ver}  cuda_available={_avail}  card sm_{_cap}", flush=True)
print(f"  kernels: {_arch}", flush=True)
if _avail != "True":
    sys.exit("no CUDA device visible after install")
# Do NOT gate on string membership of sm_<cap> in get_arch_list(). That check is a proxy and
# it FALSELY REJECTED an L4: torch 2.9.1+cu129 lists sm_86 but not sm_89, yet CUDA guarantees
# binary compatibility forward across minor revisions within a major generation, so sm_86
# cubins run on sm_89. It killed a healthy GCE instance in 10 minutes.
#
# The definitive test is to actually execute a kernel. This still catches the P100 case, where
# the matmul raises "no kernel image is available for execution on the device".
_fn = sh([sys.executable, "-c",
          "import torch;x=torch.randn(256,256,device='cuda');"
          "y=(x@x).sum().item();torch.cuda.synchronize();print('CUDA_OK',y is not None)"])
if _fn.returncode or "CUDA_OK" not in (_fn.stdout or ""):
    sys.exit(f"card sm_{_cap} cannot execute a CUDA kernel under torch {_ver} "
             f"(arch list: {_arch}). Prediction would fail. Aborting.\n"
             f"{(_fn.stderr or '')[-2000:]}")
print(f"  CUDA functional check passed on sm_{_cap}", flush=True)

# boltz ships a console script, NOT a runnable package: `python -m boltz` fails with
# "'boltz' is a package and cannot be directly executed". Resolve the entry point and prove it
# runs BEFORE building any input, so this class of error costs seconds rather than an install.
BOLTZ = shutil.which("boltz")
if not BOLTZ:
    cand = [p for p in glob.glob(sys.prefix + "/bin/boltz") + glob.glob("/usr/local/bin/boltz")]
    BOLTZ = cand[0] if cand else None
if not BOLTZ:
    sys.exit("boltz installed but no 'boltz' executable on PATH")
chk = sh([BOLTZ, "--help"])
if chk.returncode:
    sys.exit(f"{BOLTZ} --help failed rc={chk.returncode}:\n{(chk.stderr or '')[-2000:]}")
print(f"  boltz installed -> {BOLTZ}", flush=True)

# Validate every optional flag against `predict --help` before spending an install on a typo.
# Unsupported options are dropped with a note rather than failing the run; --use_msa_server is
# load-bearing (no MSA, no prediction) so its absence is fatal and must be caught here.
_ph = sh([BOLTZ, "predict", "--help"])
PHELP = (_ph.stdout or "") + (_ph.stderr or "")
OPTS = []
for _flag, _val in (("--diffusion_samples", "1"), ("--output_format", "pdb")):
    if _flag in PHELP:
        OPTS += [_flag, _val]
    else:
        print(f"  note: {_flag} not supported by this boltz build - omitting", flush=True)
if "--use_msa_server" not in PHELP:
    sys.exit("boltz predict has no --use_msa_server; cannot build an MSA. predict --help:\n"
             + PHELP[:3000])
print(f"  predict opts: {OPTS or '(none)'}", flush=True)

# Prove the cuequivariance import works BEFORE paying for MSA generation. This is the exact
# failure that wasted 518 s: the module is imported only once prediction starts, so an absent
# package looks like a healthy run right up until the GPU work begins.
if not PASCAL:
    _ce = sh([sys.executable, "-c", "import cuequivariance_torch as c; print('cueq', c.__name__)"])
    if _ce.returncode:
        sys.exit("cuequivariance_torch missing on an sm_70+ card; boltz will dispatch to it at "
                 "prediction time and fail after MSA. Install boltz[cuda].\n"
                 + (_ce.stderr or "")[-1500:])
    print(f"  {_ce.stdout.strip()}", flush=True)

def find(pattern: str) -> str:
    """Resolve a bundle file recursively.

    The launcher locates the bundle with a RECURSIVE glob for manifest.json but hands the worker
    the top-level mount, and Kaggle nests dataset files a directory below it. Joining the
    filename directly onto BUNDLE therefore missed every file and the first run died on
    FileNotFoundError for mpro_sequence.txt.
    """
    hits = sorted(glob.glob(str(BUNDLE / "**" / pattern), recursive=True))
    if not hits:
        listing = [str(p) for p in sorted(BUNDLE.rglob("*"))[:30]]
        sys.exit(f"{pattern!r} not found under {BUNDLE}. Contents: {listing}")
    return hits[0]

seq = open(find("*sequence*.txt")).read().strip()
smiles = json.load(open(find("*smiles*.json")))
labels = {}
for m in sorted(glob.glob(str(BUNDLE / "**" / "manifest*.json"), recursive=True)):
    for c in json.load(open(m))["compounds"]:
        labels[c["name"]] = c["label"]
names = sorted(n for n in smiles if n in labels)
# Kaggle caps a session at ~12 h, so a long panel must be split across sessions and the
# shards recombined. Deterministic slice by sorted name: shard k of N takes every Nth entry,
# which keeps each shard's label mix close to the panel's rather than clustering actives the
# way a contiguous block does (compound IDs are not randomly ordered by label).
_NS = int(os.environ.get("BOLTZ_NSHARDS", "1"))
_SH = int(os.environ.get("BOLTZ_SHARD", "0"))
if _NS > 1:
    _all = names
    names = _all[_SH::_NS]
    _act = sum(labels[n] for n in names) / max(len(names), 1)
    _tot = sum(labels[n] for n in _all) / max(len(_all), 1)
    print(f"  SHARD {_SH}/{_NS}: {len(names)} of {len(_all)} compounds, "
          f"active fraction {_act:.1%} vs panel {_tot:.1%}", flush=True)
print(f"  protein {len(seq)} aa (dimer A+B)   ligands {len(names)}   labelled {len(labels)}",
      flush=True)

# Arm 2: pocket conditioning. Residues within 8.0 A of the Vina box centre [9.05, 8.9, -1.51]
# on chain A -- derived from the receptor, not guessed. Verified to contain Mpro's catalytic
# dyad Cys145 (4.41 A) and His41 (5.02 A); chain B equivalents are 41-49 A away. Set
# BOLTZ_POCKET=0 to reproduce arm 1's blind configuration.
# Pocket comes from the bundle (pocket.json, written by the bundle builder from the receptor
# at the Vina box centre), not hardcoded -- the worker now serves Mpro, Factor Xa and PD-L1.
# Falls back to the Mpro set only if no pocket.json is present, for backward compatibility.
USE_POCKET = os.environ.get("BOLTZ_POCKET", "1") != "0"
_pk = sorted(glob.glob(str(BUNDLE / "**" / "pocket.json"), recursive=True))
if _pk:
    _pj = json.load(open(_pk[0]))
    POCKET = [tuple(c) for c in _pj["pocket"]]
    print(f"  pocket from bundle: {len(POCKET)} residues {POCKET[:4]}...", flush=True)
else:
    POCKET = [("A", r) for r in (164,145,41,189,165,143,142,49,144,166,141,26)]
    print("  pocket: no pocket.json in bundle, using Mpro default", flush=True)

CHAIN_IDS = ", ".join(_pj["chains"]) if _pk else "A, B"

def write_yaml(d: Path, name: str, msa: str | None):
    d.mkdir(parents=True, exist_ok=True)
    msa_line = f"      msa: {msa}\n" if msa else ""
    pocket = ""
    if USE_POCKET:
        contacts = "".join(f"        - [{ch}, {r}]\n" for ch, r in POCKET)
        pocket = ("constraints:\n  - pocket:\n      binder: L\n      contacts:\n"
                  f"{contacts}      max_distance: 5.0\n")
    (d / f"{name}.yaml").write_text(
        "version: 1\nsequences:\n"
        f"  - protein:\n      id: [{CHAIN_IDS}]\n"
        f"      sequence: {seq}\n{msa_line}"
        f"  - ligand:\n      id: L\n      smiles: '{smiles[name]}'\n"
        f"{pocket}"
        "properties:\n  - affinity:\n      binder: L\n")

def run_batch(tag, batch, msa, use_server):
    d = OUT / f"in_{tag}"; o = OUT / f"out_{tag}"
    shutil.rmtree(d, ignore_errors=True); shutil.rmtree(o, ignore_errors=True)
    for n in batch: write_yaml(d, n, msa)
    cmd = [BOLTZ, "predict", str(d), "--out_dir", str(o), *OPTS]
    if use_server: cmd.append("--use_msa_server")
    t0 = time.time(); r = sh(cmd); dt = time.time() - t0
    return dt, o, r

def harvest(o: Path):
    got, err = {}, {}
    for f in glob.glob(str(o / "**" / "affinity_*.json"), recursive=True):
        nm = os.path.basename(f)[len("affinity_"):-len(".json")]
        try:
            j = json.load(open(f))
        except Exception as e:
            err[nm] = f"unreadable affinity json: {e}"; continue
        pb, pv = j.get("affinity_probability_binary"), j.get("affinity_pred_value")
        if pb is None and pv is None:
            err[nm] = f"no affinity fields in {list(j)[:6]}"
        else:
            got[nm] = {"affinity_probability_binary": pb, "affinity_pred_value": pv}
    return got, err

print("\n== STAGE 1 probe", flush=True)
dtA, oA, rA = run_batch("a", names[:PROBE_A], None, True)
gotA, errA = harvest(oA)
if not gotA:
    sys.exit(f"probe batch A produced no affinity output in {dtA:.0f}s\n"
             f"STDERR:\n{(rA.stderr or '')[-4000:]}\nSTDOUT:\n{(rA.stdout or '')[-2000:]}")
def pick_msa():
    """Find the MSA boltz actually produced, and prove it is usable before reusing it.

    v6 globbed for *.a3m and took the first hit, which was a ColabFold scratch file
    (msa/<name>_unpaired_tmp_env/bfd.mgnify30.metaeuk30.smag30.a3m) -- a raw database
    intermediate containing NUL bytes, not the processed alignment. Every ligand in batch B
    then died on KeyError: '\\x00' inside boltz's a3m parser. Boltz writes the real MSA as
    msa/<name>_0.csv, so prefer that, and never accept anything under a tmp_env directory.
    """
    cands = [p for p in sorted(glob.glob(str(OUT / "**" / "msa" / "*.csv"), recursive=True))
             if "tmp_env" not in p]
    if not cands:
        cands = [p for p in sorted(glob.glob(str(OUT / "**" / "*.a3m"), recursive=True))
                 if "tmp_env" not in p]
    for p in cands:
        with open(p, "rb") as fh:
            head = fh.read(65536)
        if b"\x00" in head:
            print(f"  rejecting {p}: contains NUL bytes", flush=True)
            continue
        if not head.strip():
            print(f"  rejecting {p}: empty", flush=True)
            continue
        stable = OUT / ("msa_shared" + os.path.splitext(p)[1])
        shutil.copy(p, stable)
        print(f"  MSA reuse -> {stable} (from {p}, {os.path.getsize(p)} bytes)", flush=True)
        return str(stable)
    print(f"  no reusable MSA among {len(cands)} candidates; falling back to the server per batch",
          flush=True)
    return None

msa = pick_msa()
print(f"  batch A {len(gotA)}/{PROBE_A} in {dtA:.0f}s   msa reuse: {msa or 'NONE (server)'}",
      flush=True)

dtB, oB, rB = run_batch("b", names[PROBE_A:PROBE_A + PROBE_B], msa, msa is None)
gotB, errB = harvest(oB)
if not gotB and msa:
    # The reused MSA was rejected by boltz (wrong format, or an artefact we should not have
    # picked). Fall back to the server so the probe still measures a WORKING path -- that path
    # pays MSA generation per ligand, so the projection will almost certainly exceed the ceiling
    # and abort, which is the honest answer rather than a silent stall.
    print(f"  reused MSA produced nothing; retrying batch B via the server\n"
          f"  stderr: {(rB.stderr or '')[-800:]}", flush=True)
    msa = None
    dtB, oB, rB = run_batch("b2", names[PROBE_A:PROBE_A + PROBE_B], None, True)
    gotB, errB = harvest(oB)
if not gotB:
    sys.exit(f"probe batch B produced no affinity output in {dtB:.0f}s\n"
             f"STDERR:\n{(rB.stderr or '')[-4000:]}")
print(f"  MSA path in use for stage 2: {'reused file' if msa else 'SERVER per ligand (costly)'}",
      flush=True)
marginal = dtB / len(gotB)
projected = marginal * len(names) / 3600.0
print(f"  batch B {len(gotB)}/{PROBE_B} in {dtB:.0f}s", flush=True)
print(f"\n  MARGINAL {marginal:.1f} s/ligand   PROJECTED {projected:.2f} GPU-hours "
      f"for {len(names)}   CEILING {CEILING_HOURS}", flush=True)

report = {"marginal_s_per_ligand": round(marginal, 2), "projected_gpu_hours": round(projected, 3),
          "ceiling_gpu_hours": CEILING_HOURS, "n_ligands": len(names),
          "probe_a_s": round(dtA, 1), "probe_b_s": round(dtB, 1),
          "paper_estimate_s_per_ligand": 20}
if projected > CEILING_HOURS:
    report["decision"] = "ABANDONED - projected cost exceeds pre-registered ceiling"
    json.dump(report, open(OUT / f"boltz2_probe{_TAG}.json", "w"), indent=2)
    print("\n  ABORT per pre-registration. Not spending the quota.", flush=True)
    sys.exit(0)

report["decision"] = "PROCEED"
json.dump(report, open(OUT / f"boltz2_probe{_TAG}.json", "w"), indent=2)

print(f"\n== STAGE 2 full panel ({len(names)} ligands)", flush=True)
scores = dict(gotA); scores.update(gotB)
fails = dict(errA); fails.update(errB)
rest = names[PROBE_A + PROBE_B:]
CH = 25
t0 = time.time()
for i in range(0, len(rest), CH):
    chunk = rest[i:i + CH]
    dt, o, r = run_batch(f"c{i}", chunk, msa, msa is None)
    g, e = harvest(o)
    scores.update(g); fails.update(e)
    for n in chunk:
        if n not in g and n not in e:
            fails[n] = "no affinity json produced; stderr tail: " + (r.stderr or "")[-300:]
    shutil.rmtree(o, ignore_errors=True); shutil.rmtree(OUT / f"in_c{i}", ignore_errors=True)
    done = PROBE_A + PROBE_B + i + len(chunk)
    print(f"  {done}/{len(names)} ok={len(scores)} fail={len(fails)} "
          f"{(time.time()-t0)/60:.1f} min", flush=True)
    json.dump(scores, open(OUT / f"boltz2_scores{_TAG}.json", "w"))
    json.dump(fails, open(OUT / f"boltz2_failures{_TAG}.json", "w"))

json.dump(scores, open(OUT / f"boltz2_scores{_TAG}.json", "w"), indent=1)
json.dump(fails, open(OUT / f"boltz2_failures{_TAG}.json", "w"), indent=1)
print(f"\nDONE scored {len(scores)}  failed {len(fails)}", flush=True)
if fails:
    fa = [labels[n] for n in fails if n in labels]
    sa = [labels[n] for n in scores if n in labels]
    if fa and sa:
        d = sum(fa)/len(fa) - sum(sa)/len(sa)
        print(f"  DROPOUT CONTROL: active fraction failed {sum(fa)/len(fa):.1%} vs "
              f"scored {sum(sa)/len(sa):.1%}  diff {d:+.1%}  "
              f"{'VOID - class-selected' if abs(d) > 0.10 else 'within 10pp, OK'}", flush=True)
