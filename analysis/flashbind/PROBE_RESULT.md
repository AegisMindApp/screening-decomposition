# FlashBind probe: two blockers, one of them structural

**8 September 2026.** Staged probe on Kaggle GPU, kernel `oceansparx/flashbind-probe` v1.
Rules in `PREREGISTRATION.md` (`c6b4504f9`). Leakage already cleared in `LEAKAGE_CHECK.md`.

The probe was built to attempt every stage on 5 compounds and report where each breaks, rather
than test one thing per push. The Boltz-2 worker took eight pushes because each one tested a
single hypothesis. One push here named the whole failure surface.

| stage | result | |
|---|---|---|
| clone | **OK** | 13 entries, 37 s |
| checkpoints | **OK** | `binary_1.ckpt`, `binary_2.ckpt`, 43 MB each, ungated |
| gpu | **FAIL** | `no kernel image is available for execution on the device` |
| deps | **FAIL** | `No matching distribution found for torchdrug` |
| esm3_repr | **FAIL** | cascade from deps |

## Blocker 1 — Kaggle gave a P100, torch has no sm_60 kernels

> Tesla P100-PCIE-16GB with CUDA capability **sm_60** is not compatible with the current PyTorch
> installation. The current install supports sm_70 sm_75 sm_80 sm_86 sm_90 sm_100 sm_120.

Known and already solved for Boltz-2: install `torch==2.4.1 torchvision==0.19.1` from the cu121
index, the last release shipping Pascal kernels, and pin both together so pip cannot orphan
torchvision. The push API cannot request a card, so the worker must detect and adapt.
**Tractable — the code exists in `analysis/boltz2/boltz2_worker.py`.**

## Blocker 2 — torchdrug cannot be installed on Kaggle's Python at all

    ERROR: Ignored the following versions that require a different python version:
      ... 0.2.1 Requires-Python >=3.7,<3.11
    ERROR: No matching distribution found for torchdrug

torchdrug's newest release caps at **Python < 3.11**. Kaggle runs 3.11+. There is no version
that installs. This is not a pin to loosen; the package does not support the interpreter.

FlashBind needs torchdrug only to build `ligand_repr` (`torchdrug.lmdb`) — a **precomputed**
feature file consumed by `predict.py`. So the work can be split: generate ligand features off
Kaggle, ship them in the bundle, and let the GPU kernel do ESM3, docking and inference.

**But the obvious host does not work either.** Local Python 3.10 satisfies the interpreter
constraint and still fails: `Failed to build 'torch-scatter' when getting requirements to build
wheel` — torch-scatter needs a compiler toolchain matched to an installed torch, and this box
does not have one.

## Where this leaves the arm

Not blocked in principle: checkpoints are public and MIT, ESM3 is ungated, the panel bundle is
built and verified (751 ligands, 257 active, matching the panel), and leakage is ruled out for
the target. Two concrete engineering problems remain, and the second is the real one — a
dependency that installs on neither the GPU host nor this machine.

Three routes, none attempted yet:

1. **Prebuilt torch-scatter wheels.** `data.pyg.org` publishes wheels per torch/CUDA version. If
   one matches a torch that also ships sm_60, both blockers fall to the same pin.
2. **Bypass torchdrug.** Inspect `src/affinity/data/repr/torchdrug.py`. If it uses a narrow slice
   of the API — a standard molecule featuriser — reimplement that slice with RDKit and assert
   the output shape against the shipped `mf-pcba/repr/torchdrug.lmdb`.
3. **Ask the authors.** The repo is active (pushed 31 Aug 2026) and the paper is TMLR 2026.

Route 2 is the most likely to work and the easiest to get quietly wrong: a reimplemented
featuriser that produces plausible tensors of the right shape but different semantics would
give a number that looks fine and means nothing. It needs an equality check against their
shipped features before any score is read.

## What the probe cost

One Kaggle GPU session, a few minutes of an 11-hour weekly quota. Two blockers named, two stages
proven, nothing wasted on a full run that would have died at the first import.

---

# Follow-up: blocker 2 does not exist, and using torchdrug would have been wrong

**8 September 2026.** Route 1 (prebuilt wheels) was attempted and succeeded, and in doing so
showed the blocker was misdiagnosed.

## The wheels work

`data.pyg.org` publishes `torch_scatter-2.1.2+pt24cpu-cp310`. With `torch==2.4.1+cpu` plus
`torch-scatter`, `torch-cluster`, `torch-sparse`, `torch-spline-conv` from that index and
`numpy<2`, **torchdrug 0.2.1 installs and featurises** on local Python 3.10 — 67-dim node
features, 18-dim edge features.

## But FlashBind never uses torchdrug

`src/affinity/data/repr/torchdrug.py` is a **self-contained reimplementation**, headed *"Code
from source code of torchdrug"*. Across 170 lines the string "torchdrug" appears three times: two
comments and an output filename. Its imports are rdkit, lmdb, pandas, joblib, tqdm, torch —
every one installable on Kaggle's Python 3.11.

`env.yaml` lists torchdrug because the training environment had it. The inference path does not.

## And the two featurisers are not the same

| | node features |
|---|---|
| their `extract_feature` | **56** |
| torchdrug `node_feature="default"` | **67** |

Their 56 decomposes exactly as their vocabularies dictate: symbol 17+1, degree 7+1, num_Hs 7+1,
total valence 8+1, formal charge 11+1, aromatic 1. torchdrug's default adds chirality and
hybridisation terms they dropped. They also compute **no edge features**, where torchdrug returns
18.

**The checkpoint was trained on the 56-dim variant.** Had the env.yaml dependency been taken at
face value — install torchdrug, call the default featuriser — the model would have been fed
67-dim tensors: either a crash, or a silent mismatch producing scores that look like scores.

That is the failure this file warned about before either route was tried: *"a reimplemented
featuriser that produces plausible tensors of the right shape but different semantics would give
a number that looks fine and means nothing."* The warning was aimed at the wrong direction — the
risk was not in reimplementing their code, it was in using the library their config names.

## Revised plan

- **Blocker 2: dissolved.** Use their `repr/torchdrug.py` unchanged. No torchdrug, no PyG wheels,
  no split between hosts.
- **Blocker 1: unchanged.** The P100 still needs the `torch==2.4.1`/`torchvision==0.19.1` cu121
  pin from `analysis/boltz2/boltz2_worker.py`.
- The local py3.10 venv is retained only as the reference that established the 56-vs-67
  difference. It is not part of the pipeline.

---

## v3 — the affinity input side is solved; two environment blockers remain

Kernel `oceansparx/flashbind-probe` v3, Tesla P100, 48-compound working set.

| stage | result |
|---|---|
| gpu | OK — torch 2.4.1+cu121, `matmul=True` on the P100 |
| clone | OK |
| deps | OK |
| checkpoints | OK — binary_1/2 43 MB each, `fabind_plus_best_ckpt.bin` **180,353,461 B** |
| pyg_ops | **FAILED** — `No module named 'torch_scatter'` |
| esm3_repr | **OK** — 1 protein, **[298, 1536]** |
| ligand_features | OK — 48 featurised, 0 failed, dims=56 |
| ligand_repr_lmdb | OK — 48 entries |
| fabind_prep | **FAILED** — `module 'esm' has no attribute 'pretrained'` |
| fabind_dock, pocket_agreement, affinity_predict | not reached |

**What is now settled.** The affinity model's entire input surface works. ESM3 produces
exactly what `affinity_binary.yaml` declares — 298 rows for a 298-residue sequence, 1536
dims — using FlashBind's own `logits(LogitsConfig(return_embeddings=True))` call rather than
my guess at the SDK. Ligands featurise 48/48 at 56 dims through their reimplementation, and
the lmdb container holds all 48. The FABind+ checkpoint is real weights, not the LFS pointer a
`--depth 1` clone would have handed us.

**Both remaining blockers are mine, and neither is about FlashBind.**

*cp311 wheels on a cp312 image.* I hardcoded the ABI tag and verified only that the URLs
returned 200. data.pyg.org publishes cp311 **and** cp312, so the 200 never discriminated
between them — pip skipped all four wheels without erroring and `torch_scatter` was simply
absent. A reachable wheel is not an installable wheel. Same shape as every other silent
failure in this project: the check returned a clean answer to a question I had not asked.

*Two packages, one namespace.* `esm==3.2.0` (EvolutionaryScale's ESM3) and `fair-esm`
(Meta's ESM2) both install a top-level `esm`. Installing both leaves whichever pip wrote last,
and ESM3 won — so FABind+'s `esm.pretrained.esm2_t33_650M_UR50D()` found no such attribute.
The pipeline genuinely needs both models: ESM3 for the affinity head, ESM2 for the pose
generator. They cannot coexist in one interpreter, so v4 swaps packages *after* `esm3.pt` is
on disk, which costs nothing.

**Cost so far:** three sessions, ~5 GPU-minutes each, against an 11-hour weekly quota. The
staged design keeps paying — v3 cleared four unknowns and named two blockers in 310 seconds.

---

## v12 — the probe passes end to end

Kernel `oceansparx/flashbind-probe` v12, Tesla P100, 48-compound stratified probe.
**All 16 stages pass; `FAILED: []`.**

| stage | result |
|---|---|
| gpu / pin_check | torch 2.4.1+cu121, live matmul, torchvision 0.19.1+cu121 |
| pyg_ops | `scatter_add` correct, `radius_graph` 492 edges |
| esm3_repr | 1 protein, **[298, 1536]** |
| ligand_features / ligand_repr_lmdb | 48 featurised, 0 failed, **dims=56**; 48 lmdb entries |
| swap_esm | ESM3 → fair-esm, `esm2_t33_650M_UR50D` reachable |
| torchdrug | `DIMS (10, 56) EDGES (20, 3)` |
| fabind_prep / fabind_dock | 48 poses + 48 pocket index sets, **coverage 1.000** |
| pocket_agreement | median **7.89 Å**, all 48 in box |
| swap_pyg | 2.4.0 → 2.6.1, `collect_param_data` present |
| affinity_predict | both checkpoints **48/48**; binary min 0.0800, median 0.4469, max 0.9203 |

**Coverage is 1.000, not merely above the floor.** The pre-registered bar was ≥90% and no
compound was lost at any stage — no ETKDG failure, no docking failure, no scoring failure. The
dropout-bias control that voided the DiffDock arm has nothing to fire on here.

**The pocket control reproduced four times** across v9–v12: median 7.90, 7.88, 7.88, 7.89 Å,
`frac_in_box` 1.000 every run. The spread is FABind+'s own sampling and is far inside anything
that would change the reading.

### What the twelve iterations actually were

Nine of the twelve failures were environment; three were my own glue. None were FlashBind
producing a wrong answer. The environment problems clustered into one pattern worth recording:
**this pipeline is two codebases with incompatible dependency requirements**, and each half
fails loudly on the other's versions.

| conflict | FABind+ needs | affinity head needs |
|---|---|---|
| `esm` namespace | fair-esm (`esm.pretrained.esm2_t33_650M_UR50D`) | esm 3.2.0 (`esm.models.esm3`) |
| torch-geometric | ≤2.4 (`KeyError: 'complex'` on 2.6) | 2.6 (`inspector.collect_param_data`) |
| torchdrug | the real package, for `Molecule.from_smiles` | its own 56-dim reimplementation |

Both collisions are handled the same way: run each half against the version it was written
for, and swap once the first half's output is on disk. Neither half is modified.

The three glue bugs were all the same bug — **an empty or missing artefact passing as a
result.** An lmdb that existed but held nothing; a stage that returned OK on an empty file
list; a glob that missed the file because splitting the ensemble had renamed it. Every one was
caught by an assertion rather than by inspection, and the one that was not asserted (v6's
`fabind_dock`) produced a clean `rc=0` and two empty databases that the next stage's guard had
to catch.

**Cost:** twelve sessions, 5–10 GPU-minutes each, against an 11-hour weekly quota.

### What is not yet established

Nothing about the arm's question. The probe proves the plumbing. The pre-registered positive
control — descriptors reproducing **0.7654 ± 0.005** on the panel — is only readable on the
full 751, and until it passes, any AUROC from 48 compounds has no reference frame.
