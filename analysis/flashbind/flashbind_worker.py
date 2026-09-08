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
    # torch-geometric 2.4.0 is FABind+'s own README pin. FlashBind's env.yaml says 2.6.1, but
    # that is the affinity environment; on 2.6.1 FABind+'s post_optim_mol hits
    # KeyError: 'complex' in pyg's separate(), because newer Batch.get_example requires every
    # node type in _slice_dict. If predict.py then needs 2.6.1, the two halves get swapped the
    # way the two esm packages already are.
    pkgs = ["esm==3.2.0", "lmdb", "biopython", "rdkit", "einops", "omegaconf",
            "pandas", "joblib", "tqdm", "fair-esm", "torch-geometric==2.4.0",
            "torchmetrics==1.4.3", "spyrmsd", "accelerate", "mlcrate",
            "pytorch-lightning", "hydra-core", "timeout-decorator", "gemmi",
            "xxhash", "scikit-learn", "scipy", "requests"]
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
    # Checked against the imports actually present in src/affinity and scripts/, so a missing
    # one is named here rather than 12 minutes into the run. v8 lost a session to
    # timeout_decorator, which is in their env.yaml and was not in my list.
    for m in ("esm", "lmdb", "rdkit", "pandas", "joblib", "tqdm", "timeout_decorator",
              "gemmi", "xxhash", "sklearn", "scipy", "einops", "omegaconf",
              "pytorch_lightning", "torchmetrics", "torch_geometric"):
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
    """FlashBind's own call, from src/affinity/data/repr/esm3.py - logits(LogitsConfig(...)),
    not the forward_and_sample() I first guessed at.

    Mpro is a homodimer and the receptor PDB carries both chains, so the representation must
    cover all 597 residues FABind+ indexes into. Their prot_id -> sequence format cannot express
    a dimer, so the encoder is run once per chain and the outputs stacked in structure order:
    same call, same model, no fabricated peptide bond at the A/B junction. v13 ran with chain A
    alone and inference indexed past the end of a [298, 1536] tensor.
    """
    import torch
    from esm.models.esm3 import ESM3
    from esm.sdk.api import LogitsConfig, ESMProtein
    from Bio.PDB import PDBParser
    seqs = json.load(open(f"{BUNDLE}/prots.json"))
    m = ESM3.from_pretrained("esm3_sm_open_v1").to("cuda").eval()
    reprs, shapes = {}, {}
    for pid, chains in seqs.items():
        if isinstance(chains, str):
            chains = [chains]
        parts = []
        for ci, seq in enumerate(chains):
            with torch.no_grad():
                t = m.encode(ESMProtein(sequence=seq)).to("cuda")
                out = m.logits(t, LogitsConfig(return_embeddings=True))
            emb = out.embeddings.squeeze(0)[1:-1]        # strip BOS/EOS, as they do
            assert emb.shape[0] == len(seq), \
                f"{pid} chain {ci}: {emb.shape[0]} rows vs {len(seq)} residues"
            # affinity_binary.yaml declares protein_repr_dim 1536; a wrong encoder is otherwise
            # silent and the checkpoint would still run.
            assert emb.shape[1] == 1536, f"{pid} chain {ci}: {emb.shape[1]} dims, expected 1536"
            parts.append(emb.cpu())
        r = torch.cat(parts, dim=0)
        # The control for the v13 failure: the representation must have exactly one row per
        # residue FABind+ can index. Anything else is a receptor-identity mismatch, and it
        # surfaces here rather than as an out-of-bounds error deep inside inference on
        # whichever compounds happen to select a pocket in the second chain.
        st = PDBParser(QUIET=True).get_structure("r", f"{BUNDLE}/{pid}_receptor.pdb")[0]
        n_ca = sum(1 for c in st for res in c if res.id[0] == " " and "CA" in res)
        assert r.shape[0] == n_ca, \
            f"{pid}: representation has {r.shape[0]} rows, receptor has {n_ca} CA residues"
        reprs[pid] = r
        shapes[pid] = (len(chains), tuple(r.shape), n_ca)
    torch.save(reprs, f"{OUT}/esm3.pt")
    del m; torch.cuda.empty_cache()
    return f"{len(reprs)} proteins; " + "; ".join(
        f"{k}: {v[0]} chains -> {v[1]}, matches {v[2]} CA residues" for k, v in shapes.items())

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
    # torchdrug/data/rdkit/draw.py imports rdkit.Chem.Draw.mplCanvas, which current rdkit no
    # longer ships - it is a plotting module we never touch, but it sits on the import path of
    # data/__init__.py so nothing in torchdrug loads without it. torchdrug is unmaintained
    # (v0.2.1, Jul 2023; that import unchanged since 2021), so the fix is on the rdkit side.
    # 2024.3.6 is the newest cp312 build that still carries mplCanvas - a supported version pin,
    # not a patch. Safe here because ligand_features and ligand_repr_lmdb have already run.
    rr = sh([sys.executable, "-m", "pip", "install", "-q", "rdkit==2024.3.6"], timeout=1800)
    assert rr.returncode == 0, f"rdkit pin rc={rr.returncode}"
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

@stage("pin_check")
def s_pin():
    """Every pip after the Pascal bootstrap can silently move torch. If it does, the P100 loses
    its kernels and the failure surfaces as something unrelated - v6 died inside torchvision, a
    library nothing in this pipeline uses, imported transitively by torchmetrics."""
    import torch, importlib
    assert torch.__version__.startswith("2.4.1"), f"torch drifted to {torch.__version__}"
    x = torch.randn(256, 256, device="cuda"); float((x @ x).sum())
    try:
        tv = importlib.import_module("torchvision")
        tvv = tv.__version__
    except Exception as e:
        tvv = f"UNIMPORTABLE {type(e).__name__}"
    return f"torch {torch.__version__} matmul ok; torchvision {tvv}"

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
    # v6 passed this stage on the existence of two lmdbs that FABind+ created and then left
    # empty when it crashed one line later. The container is not the contents.
    import lmdb
    sdfs = sorted(p for p in outs if "ligand_sdf" in p)
    pkts = sorted(p for p in outs if "pocket_indices" in p)
    # FABind+ names its output per instance-id. We pass only --instance-id 0, so there must be
    # exactly one pair; more would mean the panel sharded, and every glob downstream takes the
    # first match - a fraction of the panel scored, reported at coverage 1.000.
    assert len(sdfs) == 1 and len(pkts) == 1, f"expected one lmdb pair, got {sdfs} {pkts}"
    counts = {}
    for o in outs:
        with lmdb.open(o, readonly=True, lock=False).begin() as t:
            counts[os.path.basename(o)] = t.stat()["entries"]
    for k, v in counts.items():
        assert v > 0, f"{k} is empty (rc={r.returncode}) - docking wrote nothing"
    cov = min(counts.values()) / NWANT
    return f"{counts}, coverage {cov:.3f} of {NWANT}"

@stage("pocket_agreement")
def s_pocket():
    """Pre-registered as a reported quantity. If FlashBind clears the docking ceiling, the two
    methods that cleared it are also the two that were not box-constrained - this is what
    separates 'property of the scoring approach' from 'property of the box'."""
    import lmdb, pickle, statistics
    from Bio.PDB import PDBParser
    ch = PDBParser(QUIET=True).get_structure("r", f"{BUNDLE}/mpro_receptor.pdb")[0]
    ca = [r["CA"].get_coord() for r in ch.get_residues() if "CA" in r]
    dbs = glob.glob(f"{FB}/FABind_plus/work/out/pocket_indices*")
    assert len(dbs) == 1, f"expected one pocket_indices lmdb, got {dbs}"
    db = dbs[0]
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
    dists = [float(d) for d in dists]   # Bio.PDB coords are numpy float32; json refuses them
    res = {"n": len(dists), "median_dist_A": round(statistics.median(dists), 2),
           "min": round(min(dists), 2), "max": round(max(dists), 2),
           "frac_in_box": round(inside/len(dists), 3), "box_centre": list(BOX_CENTRE),
           "all_dists_A": [round(d, 2) for d in sorted(dists)]}
    json.dump(res, open(f"{OUT}/pocket_agreement.json", "w"), indent=2)
    return res

@stage("swap_pyg")
def s_swappyg():
    """The two halves need different torch-geometric versions, and each fails loudly on the
    other's. FABind+'s post_optim_mol needs <=2.4 (Batch.get_example, KeyError: 'complex' on
    2.6); FlashBind's affinity EGNN calls inspector.collect_param_data, which is 2.6 API and
    absent from 2.4. Same structure as the esm collision, same handling: run each half against
    the version it was written for, and swap once the first half's output is on disk."""
    for p in ("ligand_sdf_0.lmdb", "pocket_indices_0.lmdb"):
        d = f"{FB}/FABind_plus/work/out/{p}"
        assert os.path.isdir(d), f"refusing to swap pyg before {p} exists"
    r = sh([sys.executable, "-m", "pip", "install", "-q", "torch-geometric==2.6.1"], timeout=1200)
    assert r.returncode == 0, "pyg 2.6.1 install failed"
    c = sh([sys.executable, "-c",
            "import torch_geometric as g; from torch_geometric.inspector import Inspector;"
            "print(g.__version__, hasattr(Inspector('x'), 'collect_param_data'))"], timeout=300)
    assert "True" in (c.stdout or ""), \
        f"pyg swapped but collect_param_data absent: {(c.stdout or '')} {(c.stderr or '')[-300:]}"
    return f"pyg 2.4.0 -> {c.stdout.strip()}"

@stage("affinity_predict")
def s_predict():
    """One subprocess per checkpoint, then average.

    Their ensemble path runs both checkpoints in a single process and reopens the same
    ligand lmdb for the second, which lmdb refuses ("already open in this process"). Model 1
    scored 48/48 before it hit that. Running each checkpoint in its own process is exactly
    what docs/predict.md defines the ensemble to be - "the average of all model predictions" -
    so this changes nothing about the result, only about the process boundary.

    --devices 1, no torchrun: their predict_binary.sh is --nproc_per_node=4 and would hang on
    one P100 rather than fail."""
    import shutil, statistics
    work = f"{FB}/FABind_plus/work/out"
    sdfs = glob.glob(f"{work}/ligand_sdf*"); pkts = glob.glob(f"{work}/pocket_indices*")
    assert len(sdfs) == 1 and len(pkts) == 1, f"expected one lmdb pair, got {sdfs} {pkts}"
    sdf, pkt = sdfs[0], pkts[0]
    os.makedirs(f"{OUT}/data/pdb", exist_ok=True)
    shutil.copy(f"{BUNDLE}/mpro_receptor.pdb", f"{OUT}/data/pdb/mpro.pdb")
    # Score only what was actually docked. One compound failed conformer generation, so its id
    # is in id_work.json with no pose in the lmdb; their dataset raises on the missing key and
    # the batch then dies in collation, taking the whole run with it. Deriving the id list from
    # the lmdb keys makes the dropped set explicit instead of fatal.
    import lmdb as _lmdb
    with _lmdb.open(sdf, readonly=True, lock=False).begin() as t:
        have = {k.decode() for k, _ in t.cursor()}
    ids = json.load(open(f"{OUT}/id_work.json"))
    keep = [i for i in ids if i in have]
    dropped = [i for i in ids if i not in have]
    assert len(keep) >= 0.90 * len(ids), \
        f"only {len(keep)}/{len(ids)} docked - below the pre-registered 90% coverage floor"
    lab = json.load(open(f"{BUNDLE}/labels.json"))
    if dropped:
        da = [lab.get(i.split("_", 1)[1]) for i in dropped]
        print(f"  dropped {len(dropped)} undocked: {dropped[:5]} labels={da[:5]}", flush=True)
    json.dump(keep, open(f"{OUT}/id_predict.json", "w"))
    env_pp = os.environ.get("PYTHONPATH", "")
    os.environ["PYTHONPATH"] = f"{FB}/src" + (f":{env_pp}" if env_pp else "")

    per_model = []
    for i, ck in enumerate(("binary_1.ckpt", "binary_2.ckpt")):
        od = f"{OUT}/binary_{i}"
        r = sh([sys.executable, f"{FB}/scripts/predict.py",
                "--data", f"{OUT}/id_predict.json", "--structure", f"{OUT}/data/pdb",
                "--structure_type", "pdb", "--ligand", sdf, "--ligand_type", "sdf",
                "--pocket_indices", pkt, "--protein_repr", f"{OUT}/esm3.pt",
                "--ligand_repr", f"{OUT}/data/repr/torchdrug.lmdb",
                "--distance_threshold", "20.0", "--out_dir", od, "--devices", "1",
                "--affinity_checkpoint", f"{OUT}/checkpoints/{ck}"], timeout=10800)
        # Their writer names it affinity_predictions_{i}.json only in ensemble mode; with a
        # single checkpoint it is affinity_predictions.json, and "predictions_*" misses it.
        hits = [h for h in glob.glob(f"{od}/**/affinity_predictions*.json", recursive=True)
                if "ensemble" not in os.path.basename(h)]
        assert hits, f"{ck}: no predictions written (rc={r.returncode})"
        d = json.load(open(hits[0]))
        ok = {k: v["binary"] for k, v in d.items()
              if v.get("status") == "success" and "binary" in v}
        assert len(ok) >= len(keep) * 0.9, f"{ck}: only {len(ok)} of {len(keep)} scored"
        per_model.append(ok)
        print(f"  {ck}: {len(ok)}/{NWANT} scored", flush=True)

    shared = set(per_model[0]) & set(per_model[1])
    assert len(shared) >= NWANT * 0.5, f"only {len(shared)} keys shared across checkpoints"
    ens = {k: sum(m[k] for m in per_model) / len(per_model) for k in shared}
    json.dump({"ensemble": ens,
               "per_model": [{k: v for k, v in m.items()} for m in per_model]},
              open(f"{OUT}/flashbind_scores.json", "w"), indent=1)
    # Emitted into the log as well as to disk. /kaggle/working holds the FlashBind clone, so
    # the output-file listing runs to thousands of entries and paginates; the log is the
    # reliable channel back. Markers so the parse is exact rather than a regex over prose.
    print("SCORES_BEGIN", flush=True)
    print(json.dumps({"ensemble": ens, "per_model": per_model}), flush=True)
    print("SCORES_END", flush=True)
    vals = sorted(ens.values())
    return (f"{len(ens)}/{NWANT} panel scored, coverage {len(ens)/NWANT:.3f} "
            f"({len(dropped)} undocked); "
            f"binary min {vals[0]:.4f} median {statistics.median(vals):.4f} max {vals[-1]:.4f}")

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
    s_esm(); s_ligand(); s_ligrepr(); s_swap(); s_torchdrug(); s_pin()
    s_fbprep(); s_fbdock(); s_pocket(); s_swappyg(); s_predict()

    json.dump(STAGES, open(f"{OUT}/flashbind_stages.json", "w"), indent=1)
    ok = [k for k, v in STAGES.items() if v["ok"]]
    bad = [k for k, v in STAGES.items() if not v["ok"]]
    print(f"\n{'='*62}\nPASSED: {ok}\nFAILED: {bad}")
    print("A failed stage here is information, not a wasted run: it names the blocker.")
