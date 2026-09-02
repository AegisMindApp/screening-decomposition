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
OUT = Path("/kaggle/working")
CEILING_HOURS = 6.0          # pre-registered abort threshold
PROBE_A, PROBE_B = 4, 8

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
if any(g in _smi for g in ("P100", "P40", "K80")):
    print("  Pascal-or-older; installing torch 2.4.1+cu121 (still ships sm_60)", flush=True)
    _r = sh([sys.executable, "-m", "pip", "install", "-q", "torch==2.4.1",
             "--index-url", "https://download.pytorch.org/whl/cu121"], timeout=3600)
    if _r.returncode:
        sys.exit(f"pascal torch bootstrap failed rc={_r.returncode}:\n{(_r.stderr or '')[-2000:]}")
    print("  torch 2.4.1+cu121 installed", flush=True)

print("== install", flush=True)
r = sh(f"{sys.executable} -m pip install -q boltz")
if r.returncode:
    sys.exit("pip install boltz failed:\n" + (r.stderr or "")[-3000:])

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
if f"sm_{_cap}" not in _arch.split(","):
    sys.exit(f"torch {_ver} has no sm_{_cap} kernels (has {_arch}). The card cannot run this "
             f"build and the run would fail once prediction starts. Aborting before spending "
             f"the quota.")

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

seq = open(find("mpro_sequence.txt")).read().strip()
smiles = json.load(open(find("mpro_smiles.json")))
labels = {}
for m in sorted(glob.glob(str(BUNDLE / "**" / "manifest*.json"), recursive=True)):
    for c in json.load(open(m))["compounds"]:
        labels[c["name"]] = c["label"]
names = sorted(n for n in smiles if n in labels)
print(f"  protein {len(seq)} aa (dimer A+B)   ligands {len(names)}   labelled {len(labels)}",
      flush=True)

def write_yaml(d: Path, name: str, msa: str | None):
    d.mkdir(parents=True, exist_ok=True)
    msa_line = f"      msa: {msa}\n" if msa else ""
    (d / f"{name}.yaml").write_text(
        "version: 1\nsequences:\n"
        "  - protein:\n      id: [A, B]\n"
        f"      sequence: {seq}\n{msa_line}"
        f"  - ligand:\n      id: L\n      smiles: '{smiles[name]}'\n"
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
a3m = sorted(glob.glob(str(OUT / "**" / "*.a3m"), recursive=True))
msa = a3m[0] if a3m else None
print(f"  batch A {len(gotA)}/{PROBE_A} in {dtA:.0f}s   msa reuse: {msa or 'NONE (server)'}",
      flush=True)

dtB, oB, rB = run_batch("b", names[PROBE_A:PROBE_A + PROBE_B], msa, msa is None)
gotB, errB = harvest(oB)
if not gotB:
    sys.exit(f"probe batch B produced no affinity output in {dtB:.0f}s\n"
             f"STDERR:\n{(rB.stderr or '')[-4000:]}")
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
    json.dump(report, open(OUT / "boltz2_probe.json", "w"), indent=2)
    print("\n  ABORT per pre-registration. Not spending the quota.", flush=True)
    sys.exit(0)

report["decision"] = "PROCEED"
json.dump(report, open(OUT / "boltz2_probe.json", "w"), indent=2)

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
    json.dump(scores, open(OUT / "boltz2_scores.json", "w"))
    json.dump(fails, open(OUT / "boltz2_failures.json", "w"))

json.dump(scores, open(OUT / "boltz2_scores.json", "w"), indent=1)
json.dump(fails, open(OUT / "boltz2_failures.json", "w"), indent=1)
print(f"\nDONE scored {len(scores)}  failed {len(fails)}", flush=True)
if fails:
    fa = [labels[n] for n in fails if n in labels]
    sa = [labels[n] for n in scores if n in labels]
    if fa and sa:
        d = sum(fa)/len(fa) - sum(sa)/len(sa)
        print(f"  DROPOUT CONTROL: active fraction failed {sum(fa)/len(fa):.1%} vs "
              f"scored {sum(sa)/len(sa):.1%}  diff {d:+.1%}  "
              f"{'VOID - class-selected' if abs(d) > 0.10 else 'within 10pp, OK'}", flush=True)
