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

print("== install", flush=True)
r = sh(f"{sys.executable} -m pip install -q boltz")
if r.returncode:
    sys.exit("pip install boltz failed:\n" + (r.stderr or "")[-3000:])
print("  boltz installed", flush=True)

seq = (BUNDLE / "mpro_sequence.txt").read_text().strip()
smiles = json.load(open(BUNDLE / "mpro_smiles.json"))
labels = {}
for m in glob.glob(str(BUNDLE / "manifest*.json")):
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
    cmd = [sys.executable, "-m", "boltz", "predict", str(d), "--out_dir", str(o),
           "--diffusion_samples", "1", "--output_format", "pdb"]
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
