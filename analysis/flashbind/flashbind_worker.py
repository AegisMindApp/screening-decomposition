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

# Detect the card with nvidia-smi, not torch: Kaggle hands out a P100 (sm_60) and the push API
# cannot request otherwise. Modern torch ships no kernels below sm_70, which is exactly how the
# v1 probe died. torch 2.4.1+cu121 is the last pair that still carries Pascal kernels, and
# torchvision must move with it or its C++ ops fail to register.
_smi = subprocess.run("nvidia-smi --query-gpu=name --format=csv,noheader", shell=True,
                      capture_output=True, text=True).stdout
print(f"== accelerator: {_smi.strip() or '(none)'}", flush=True)
if any(g in _smi for g in ("P100", "P40", "K80")):
    print("  Pascal-or-older; pinning torch==2.4.1 torchvision==0.19.1 (cu121, ships sm_60)",
          flush=True)
    _r = subprocess.run([sys.executable, "-m", "pip", "install", "-q", "torch==2.4.1",
                         "torchvision==0.19.1", "--index-url",
                         "https://download.pytorch.org/whl/cu121"],
                        capture_output=True, text=True, timeout=3600)
    if _r.returncode:
        sys.exit(f"pascal bootstrap failed rc={_r.returncode}\n{(_r.stderr or '')[-1500:]}")
    print("  pinned", flush=True)

BUNDLE = next((p for p in glob.glob("/kaggle/input/**/prots.json", recursive=True)), None)
BUNDLE = os.path.dirname(BUNDLE) if BUNDLE else "analysis/flashbind/bundle"
OUT = os.environ.get("FLASHBIND_OUT", "/kaggle/working")
FULL = os.environ.get("FLASHBIND_FULL") == "1"
os.makedirs(OUT, exist_ok=True)
FB = f"{OUT}/FlashBind"
NWANT = 0
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
    # NOT torchdrug. FlashBind's repr/torchdrug.py is a self-contained rdkit reimplementation
    # emitting 56-dim node features; real torchdrug's default emits 67 and the checkpoint was
    # trained on 56. env.yaml lists torchdrug because the training environment had it.
    pkgs = ["esm==3.2.0", "lmdb", "biopython", "rdkit", "einops", "omegaconf",
            "pandas", "joblib", "tqdm", "fair-esm", "torch-geometric==2.6.1",
            "torchmetrics", "spyrmsd", "accelerate", "mlcrate",
            "pytorch-lightning", "hydra-core"]
    r = sh([sys.executable, "-m", "pip", "install", "-q"] + pkgs, timeout=3000)
    # pyg's compiled ops must match the torch we pinned, not the torch Kaggle shipped.
    # pyg publishes against 2.4.0; 2.4.1 is ABI-compatible with it.
    # Derive the ABI tag from the interpreter that is actually running. v3 hardcoded cp311
    # against a cp312 image; every URL still returned 200, because data.pyg.org publishes both.
    # A reachable wheel is not an installable one, and pip skipped all four in silence.
    cp = f"cp{sys.version_info.major}{sys.version_info.minor}"
    W = "https://data.pyg.org/whl/torch-2.4.0+cu121"
    wheels = [f"{W}/{n}+pt24cu121-{cp}-{cp}-linux_x86_64.whl" for n in
              ("torch_scatter-2.1.2", "torch_sparse-0.6.18", "torch_cluster-1.6.3",
               "torch_spline_conv-1.2.2")]
    rp = sh([sys.executable, "-m", "pip", "install", "--no-cache-dir"] + wheels, timeout=3000)
    assert rp.returncode == 0, f"pyg wheels ({cp}) failed rc={rp.returncode}"
    import importlib
    got = {}
    for m in ("esm", "lmdb", "rdkit", "pandas", "joblib", "tqdm"):
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
    # FABind+'s checkpoint is git-lfs in the repo, so a --depth 1 clone yields a ~130-byte
    # pointer that loads as a file and fails as weights. Pull it from HF and assert on size.
    for f in ("fabind_plus_best_ckpt.bin",):
        p = hf_hub_download("KyGao/FABind_plus_model", f, local_dir=f"{OUT}/fabind_ckpt")
        sz = os.path.getsize(p)
        assert sz > 50_000_000, f"{f} is {sz} bytes - that is an LFS pointer, not weights"
        got.append((f, sz))
    return got

@stage("pyg_ops")
def s_pyg():
    """torch-scatter that imports but has no sm_60 kernels fails at the first op, deep inside
    FABind+ where the traceback is useless. Run the op here, on cuda, where it is readable."""
    import torch
    from torch_scatter import scatter_add
    from torch_cluster import radius_graph
    x = torch.ones(64, 8, device="cuda")
    idx = torch.arange(64, device="cuda") % 4
    out = scatter_add(x, idx, dim=0, dim_size=4)
    assert out.shape == (4, 8) and float(out.sum()) == 512.0, f"scatter_add wrong: {out.sum()}"
    e = radius_graph(torch.randn(50, 3, device="cuda"), r=1.5)
    return f"scatter_add ok, radius_graph {e.shape[1]} edges"

@stage("esm3_repr")
def s_esm():
    """FlashBind's own call, from src/affinity/data/repr/esm3.py - not my guess at the SDK.
    v2 died on forward_and_sample(); the shipped path is logits(LogitsConfig(...))."""
    import torch
    from esm.models.esm3 import ESM3
    from esm.sdk.api import LogitsConfig, ESMProtein
    seqs = json.load(open(f"{BUNDLE}/prots.json"))
    m = ESM3.from_pretrained("esm3_sm_open_v1").to("cuda").eval()
    reprs = {}
    for pid, seq in seqs.items():
        with torch.no_grad():
            t = m.encode(ESMProtein(sequence=seq)).to("cuda")
            out = m.logits(t, LogitsConfig(return_embeddings=True))
        emb = out.embeddings.squeeze(0)[1:-1]           # strip BOS/EOS, as they do
        assert emb.shape[0] == len(seq), f"{pid}: {emb.shape[0]} rows vs {len(seq)} residues"
        # affinity_binary.yaml declares protein_repr_dim 1536. The mirror of the 56-dim
        # ligand assertion: a wrong encoder here is silent, and the checkpoint would still run.
        assert emb.shape[1] == 1536, f"{pid}: {emb.shape[1]} dims, expected 1536"
        reprs[pid] = emb.cpu()
    torch.save(reprs, f"{OUT}/esm3.pt")
    del m; torch.cuda.empty_cache()
    return f"{len(reprs)} proteins, [{emb.shape[0]}, 1536]"

@stage("ligand_features")
def s_ligand():
    """Run FlashBind's own featuriser. Assert 56 dims: 67 would mean torchdrug's default
    leaked in, and the checkpoint expects 56."""
    import importlib.util, torch
    src = f"{OUT}/FlashBind/src/affinity/data/repr/torchdrug.py"
    assert os.path.exists(src), f"missing {src}"
    spec = importlib.util.spec_from_file_location("fb_repr", src)
    mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    smi = json.load(open(f"{OUT}/smiles_work.json"))
    feats, bad = {}, []
    for name, s in smi.items():
        try:
            f = mod.extract_feature(mod.read_smiles(s))
            assert f.shape[1] == 56, f"{name}: {f.shape[1]} dims, expected 56"
            feats[name] = f
        except Exception as e:
            bad.append((name, str(e)[:60]))
    assert feats, f"no compound featurised; first errors {bad[:3]}"
    torch.save(feats, f"{OUT}/ligand_features.pt")
    return f"{len(feats)} featurised, {len(bad)} failed, dims=56"

# The docking box the six Vina/gnina scores shared (bundle_CHEMBL4523582_7VU6 manifests).
# FABind+ picks its own pocket; PREREGISTRATION.md says report the distance, do not gate on it.
BOX_CENTRE = (9.05, 8.90, -1.51)
BOX_HALF = 11.0

@stage("ligand_repr_lmdb")
def s_ligrepr():
    """Their CLI, writing the lmdb predict.py expects. The per-compound assertion already ran
    in ligand_features; this is the container."""
    n = sh([sys.executable, f"{FB}/src/affinity/data/repr/torchdrug.py",
            "--input_json", f"{OUT}/smiles_work.json",
            "--output_lmdb", f"{OUT}/data/repr/torchdrug.lmdb"], timeout=1800)
    assert n.returncode == 0, "torchdrug.py failed"
    import lmdb
    with lmdb.open(f"{OUT}/data/repr/torchdrug.lmdb", readonly=True, lock=False).begin() as t:
        k = t.stat()["entries"]
    assert k >= NWANT * 0.75, f"only {k} of {NWANT} ligand reprs written"
    return f"{k} entries"

@stage("swap_esm")
def s_swap():
    """`esm==3.2.0` (EvolutionaryScale, ESM3) and `fair-esm` (Meta, ESM2) both install a top-level
    `esm` package, so only one can exist at a time. The affinity head needs ESM3; FABind+ calls
    esm.pretrained.esm2_t33_650M_UR50D() from fair-esm. They cannot coexist, and v3 quietly kept
    whichever pip wrote last.

    Runs only after esm3_repr has written esm3.pt to disk, so the swap costs nothing."""
    assert os.path.exists(f"{OUT}/esm3.pt"), "refusing to remove ESM3 before its output exists"
    sh([sys.executable, "-m", "pip", "uninstall", "-y", "-q", "esm", "fair-esm"], timeout=900)
    r = sh([sys.executable, "-m", "pip", "install", "-q", "fair-esm"], timeout=900)
    assert r.returncode == 0, "fair-esm install failed"
    # Checked in a subprocess, not by reloading: this process already imported ESM3's `esm`,
    # and importlib.reload would re-execute a spec pointing at files pip has just deleted.
    # A subprocess is also how FABind+ actually runs, so it tests the real condition.
    c = sh([sys.executable, "-c",
            "import esm; print(hasattr(esm.pretrained, 'esm2_t33_650M_UR50D'))"], timeout=300)
    assert "True" in (c.stdout or ""), \
        f"fair-esm installed but esm2_t33_650M_UR50D unreachable: {(c.stderr or '')[-300:]}"
    return "esm3 -> fair-esm; esm2_t33_650M_UR50D reachable in a fresh interpreter"

@stage("torchdrug")
def s_torchdrug():
    """The real torchdrug, for FABind+ only. Distinct from the affinity module, which ships a
    self-contained reimplementation of the same 56-dim featuriser and needs no torchdrug at all.

    torchdrug 0.2.1 declares python<3.11 and the image is 3.12. The cap is metadata; whether
    the code runs on 3.12 is a separate question and the only one that matters, so install past
    the cap and then make the call FABind+ makes. --no-deps because its pins would drag in a
    torch without Pascal kernels and undo the pin the whole run rests on."""
    sh([sys.executable, "-m", "pip", "install", "-q", "decorator", "easydict", "ninja",
        "jinja2", "pyyaml", "networkx", "matplotlib"], timeout=1200)
    r = sh([sys.executable, "-m", "pip", "install", "-q", "--no-deps",
            "--ignore-requires-python", "torchdrug==0.2.1"], timeout=1200)
    assert r.returncode == 0, f"torchdrug install rc={r.returncode}"
    # Exactly the call in FABind_plus/fabind/utils/inference_mol_utils.py:81. If this works the
    # blocker is gone; if it does not, it names the incompatibility rather than surfacing it
    # 30 minutes later inside docking.
    probe = ("from torchdrug import data as td;"
             "m=td.Molecule.from_smiles('O=C(CCl)c1ccccc1', node_feature='property_prediction');"
             "print('DIMS', tuple(m.node_feature.shape), 'EDGES', tuple(m.edge_list.shape))")
    c = sh([sys.executable, "-c", probe], timeout=900)
    assert "DIMS" in (c.stdout or ""), \
        f"torchdrug installed but unusable on py3.12: {(c.stderr or '')[-600:]}"
    return [l for l in c.stdout.splitlines() if l.startswith("DIMS")][0]

@stage("fabind_prep")
def s_fbprep():
    """Protein (ESM2-650M) and ligand (ETKDG confs) preprocessing. Conformer generation is the
    likeliest source of the coverage failures the 90% floor is there to catch."""
    import csv, shutil
    smi = json.load(open(f"{OUT}/smiles_work.json"))
    os.makedirs(f"{FB}/FABind_plus/work/pdb", exist_ok=True)
    os.makedirs(f"{FB}/FABind_plus/work/repr", exist_ok=True)
    shutil.copy(f"{BUNDLE}/mpro_receptor.pdb", f"{FB}/FABind_plus/work/pdb/mpro.pdb")
    with open(f"{FB}/FABind_plus/work/smiles.csv", "w", newline="") as f:
        w = csv.writer(f); w.writerow(["smiles", "ligand_id"])
        for lid, s in smi.items(): w.writerow([s, lid])
    with open(f"{FB}/FABind_plus/work/data.csv", "w", newline="") as f:
        w = csv.writer(f); w.writerow(["smiles", "prot_id", "ligand_id"])
        for lid, s in smi.items(): w.writerow([s, "mpro", lid])
    fb = f"{FB}/FABind_plus/fabind"
    r1 = sh(f"cd {fb} && {sys.executable} inference_preprocess_protein.py "
            f"--pdb_file_dir ../work/pdb --save_pt_dir ../work/repr", timeout=3600)
    assert os.path.exists(f"{FB}/FABind_plus/work/repr/processed_protein.pt"), "no protein pt"
    r2 = sh(f"cd {fb} && {sys.executable} inference_preprocess_mol_confs.py "
            f"--index_csv ../work/smiles.csv --save_mols_dir ../work/repr "
            f"--num_threads 0 --resume", timeout=3600)
    mols = glob.glob(f"{FB}/FABind_plus/work/repr/*mol*")
    # v4 reported this stage OK with an empty artefact list, and the real failure (torchdrug
    # missing) only surfaced one stage later. An empty output is not a passing stage.
    assert mols, f"no molecule artefacts written (rc={r2.returncode}); conformer stage produced nothing"
    return f"protein pt ok; mol artefacts {[os.path.basename(m) for m in mols][:4]}"

@stage("fabind_dock")
def s_fbdock():
    fb = f"{FB}/FABind_plus/fabind"
    r = sh(f"cd {fb} && {sys.executable} inference_regression_fabind.py "
           f"--ckpt {OUT}/fabind_ckpt/fabind_plus_best_ckpt.bin --batch_size 4 "
           f"--post-optim --write-mol-to-file "
           f"--sdf-output-path-post-optim ../work/out "
           f"--index-csv ../work/data.csv --preprocess-dir ../work/repr --instance-id 0",
           timeout=7200)
    outs = glob.glob(f"{FB}/FABind_plus/work/out/*")
    assert any("ligand_sdf" in o for o in outs), f"no ligand_sdf lmdb; got {outs[:6]}"
    assert any("pocket_indices" in o for o in outs), f"no pocket_indices lmdb; got {outs[:6]}"
    return f"{[os.path.basename(o) for o in outs]}"

@stage("pocket_agreement")
def s_pocket():
    """Pre-registered as a reported quantity. If FlashBind clears the docking ceiling, the two
    methods that cleared it are also the two that were not box-constrained - this is what
    separates 'property of the scoring approach' from 'property of the box'."""
    import lmdb, pickle, statistics
    from Bio.PDB import PDBParser
    ch = PDBParser(QUIET=True).get_structure("r", f"{BUNDLE}/mpro_receptor.pdb")[0]
    ca = [r["CA"].get_coord() for r in ch.get_residues() if "CA" in r]
    db = next(p for p in glob.glob(f"{FB}/FABind_plus/work/out/pocket_indices*"))
    dists, inside = [], 0
    with lmdb.open(db, readonly=True, lock=False).begin() as t:
        for k, v in t.cursor():
            idx = pickle.loads(v)
            pts = [ca[i] for i in idx if i < len(ca)]
            if not pts: continue
            c = [sum(p[j] for p in pts)/len(pts) for j in range(3)]
            d = sum((c[j]-BOX_CENTRE[j])**2 for j in range(3)) ** 0.5
            dists.append(d)
            inside += all(abs(c[j]-BOX_CENTRE[j]) <= BOX_HALF for j in range(3))
    assert dists, "pocket_indices lmdb read but no centroid computed - empty is not a result"
    res = {"n": len(dists), "median_dist_A": round(statistics.median(dists), 2),
           "min": round(min(dists), 2), "max": round(max(dists), 2),
           "frac_in_box": round(inside/len(dists), 3), "box_centre": BOX_CENTRE}
    json.dump(res, open(f"{OUT}/pocket_agreement.json", "w"), indent=2)
    return res

@stage("affinity_predict")
def s_predict():
    """--devices 1, no torchrun: their predict_binary.sh is --nproc_per_node=4 and would hang
    at rendezvous on one P100 rather than fail."""
    import shutil
    work = f"{FB}/FABind_plus/work/out"
    sdf = next(p for p in glob.glob(f"{work}/ligand_sdf*"))
    pkt = next(p for p in glob.glob(f"{work}/pocket_indices*"))
    os.makedirs(f"{OUT}/data/pdb", exist_ok=True)
    shutil.copy(f"{BUNDLE}/mpro_receptor.pdb", f"{OUT}/data/pdb/mpro.pdb")
    r = sh([sys.executable, f"{FB}/scripts/predict.py",
            "--data", f"{OUT}/id_work.json", "--structure", f"{OUT}/data/pdb",
            "--structure_type", "pdb", "--ligand", sdf, "--ligand_type", "sdf",
            "--pocket_indices", pkt, "--protein_repr", f"{OUT}/esm3.pt",
            "--ligand_repr", f"{OUT}/data/repr/torchdrug.lmdb",
            "--distance_threshold", "20.0", "--out_dir", f"{OUT}/binary", "--devices", "1",
            "--affinity_checkpoint", f"{OUT}/checkpoints/binary_1.ckpt",
            f"{OUT}/checkpoints/binary_2.ckpt"], timeout=10800)
    f = f"{OUT}/binary/affinity_predictions_ensemble.json"
    assert os.path.exists(f), f"no ensemble predictions (rc={r.returncode})"
    preds = json.load(open(f))
    ok = {k: v for k, v in preds.items() if v.get("status") == "success"}
    # An empty or thin prediction dict is what a protein/ligand key-format mismatch looks like,
    # and it reads exactly like a model that declined to score. Same guard as residual_band.py.
    assert len(ok) >= NWANT * 0.5, f"only {len(ok)} of {NWANT} scored - check key formats"
    json.dump(preds, open(f"{OUT}/flashbind_scores.json", "w"), indent=2)
    return f"{len(ok)}/{NWANT} scored, coverage {len(ok)/NWANT:.3f}"

if __name__ == "__main__":
    print(f"bundle={BUNDLE} full={FULL}")
    # One working set for every stage, written once, so the FABind+ csv, the ligand lmdb and
    # the id list cannot drift apart. A key mismatch between them scores nothing and looks
    # like a model that declined to predict.
    _s = json.load(open(f"{BUNDLE}/{'smiles.json' if FULL else 'smiles_probe.json'}"))
    _i = json.load(open(f"{BUNDLE}/{'id.json' if FULL else 'id_probe.json'}"))
    assert {i.split("_", 1)[1] for i in _i} == set(_s), "id list and smiles disagree"
    json.dump(_s, open(f"{OUT}/smiles_work.json", "w"))
    json.dump(_i, open(f"{OUT}/id_work.json", "w"))
    NWANT = len(_s)
    print(f"working set: {NWANT} compounds")

    s_gpu(); s_clone(); s_deps(); s_ckpt(); s_pyg()
    s_esm(); s_ligand(); s_ligrepr(); s_swap(); s_torchdrug()
    s_fbprep(); s_fbdock(); s_pocket(); s_predict()

    json.dump(STAGES, open(f"{OUT}/flashbind_stages.json", "w"), indent=1)
    ok = [k for k, v in STAGES.items() if v["ok"]]
    bad = [k for k, v in STAGES.items() if not v["ok"]]
    print(f"\n{'='*62}\nPASSED: {ok}\nFAILED: {bad}")
    print("A failed stage here is information, not a wasted run: it names the blocker.")
