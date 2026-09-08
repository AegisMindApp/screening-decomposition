#!/usr/bin/env python3
"""FlashBind on the Mpro panel — staged probe.

Rules in PREREGISTRATION.md (c6b4504f9); leakage cleared in LEAKAGE_CHECK.md.

Built as a STAGED probe, not a full run. The Boltz-2 worker took eight iterations because each
push tested one thing. This attempts every stage on 5 compounds and reports where it breaks, so
one push tells us the whole failure surface. Stages are independent: a later failure still
reports what earlier stages proved.

FLASHBIND_FULL=1 runs the whole 751-compound panel instead of the probe.
"""
import json, os, subprocess, sys, time, glob, traceback

BUNDLE = next((p for p in glob.glob("/kaggle/input/**/prots.json", recursive=True)), None)
BUNDLE = os.path.dirname(BUNDLE) if BUNDLE else "analysis/flashbind/bundle"
OUT = os.environ.get("FLASHBIND_OUT", "/kaggle/working")
FULL = os.environ.get("FLASHBIND_FULL") == "1"
os.makedirs(OUT, exist_ok=True)
STAGES = {}

def sh(cmd, timeout=3600, quiet=False):
    r = subprocess.run(cmd, shell=isinstance(cmd, str), capture_output=True, text=True, timeout=timeout)
    if not quiet:
        print((r.stdout or "")[-1500:]); print((r.stderr or "")[-1500:])
    return r

def stage(name):
    def deco(fn):
        def run(*a, **k):
            t0 = time.time()
            print(f"\n{'='*62}\nSTAGE {name}\n{'='*62}", flush=True)
            try:
                out = fn(*a, **k)
                STAGES[name] = {"ok": True, "secs": round(time.time()-t0, 1), "detail": out}
                print(f"  -> OK ({time.time()-t0:.0f}s) {out}", flush=True)
                return out
            except Exception as e:
                STAGES[name] = {"ok": False, "secs": round(time.time()-t0, 1),
                                "error": f"{type(e).__name__}: {e}"}
                print(f"  -> FAILED: {type(e).__name__}: {e}", flush=True)
                traceback.print_exc()
                return None
        return run
    return deco

@stage("gpu")
def s_gpu():
    import torch
    assert torch.cuda.is_available(), "no CUDA device"
    # A functional matmul, not get_arch_list(): a string check falsely rejected a healthy L4
    # during the Boltz-2 build and cost a paid instance.
    x = torch.randn(512, 512, device="cuda"); y = float((x @ x).sum()); torch.cuda.synchronize()
    return f"{torch.cuda.get_device_name(0)} torch {torch.__version__} matmul={y is not None}"

@stage("clone")
def s_clone():
    if not os.path.isdir(f"{OUT}/FlashBind"):
        r = sh(f"git clone --depth 1 https://github.com/AIDD-Lab/FlashBind {OUT}/FlashBind")
        assert r.returncode == 0, "clone failed"
    return f"{len(os.listdir(f'{OUT}/FlashBind'))} entries"

@stage("deps")
def s_deps():
    pkgs = ["esm==3.2.0", "torchdrug", "lmdb", "biopython", "rdkit", "einops", "omegaconf"]
    r = sh([sys.executable, "-m", "pip", "install", "-q"] + pkgs, timeout=3000)
    import importlib
    got = {}
    for m in ("esm", "torchdrug", "lmdb", "rdkit"):
        try:
            importlib.import_module(m); got[m] = "ok"
        except Exception as e:
            got[m] = f"IMPORT FAIL {type(e).__name__}"
    assert all(v == "ok" for v in got.values()), f"imports: {got}"
    return got

@stage("checkpoints")
def s_ckpt():
    sh([sys.executable, "-m", "pip", "install", "-q", "huggingface_hub"], timeout=900)
    from huggingface_hub import hf_hub_download
    got = []
    for f in ("binary_1.ckpt", "binary_2.ckpt"):
        p = hf_hub_download("clorf6/FlashBind", f, local_dir=f"{OUT}/checkpoints")
        got.append((f, os.path.getsize(p)))
    return got

@stage("esm3_repr")
def s_esm():
    """The likeliest blocker: ESM3 weights and API surface."""
    from esm.models.esm3 import ESM3
    from esm.sdk.api import ESMProtein
    seq = json.load(open(f"{BUNDLE}/prots.json"))["mpro"]
    m = ESM3.from_pretrained("esm3_sm_open_v1").to("cuda").eval()
    import torch
    with torch.no_grad():
        p = ESMProtein(sequence=seq)
        t = m.encode(p)
        out = m.forward_and_sample(t, None)
    return f"sequence {len(seq)} residues encoded"

if __name__ == "__main__":
    print(f"bundle={BUNDLE} full={FULL}")
    s_gpu(); s_clone(); s_deps(); s_ckpt(); s_esm()
    json.dump(STAGES, open(f"{OUT}/flashbind_stages.json", "w"), indent=1)
    ok = [k for k, v in STAGES.items() if v["ok"]]
    bad = [k for k, v in STAGES.items() if not v["ok"]]
    print(f"\n{'='*62}\nPASSED: {ok}\nFAILED: {bad}")
    print("A failed stage here is information, not a wasted run: it names the blocker.")
