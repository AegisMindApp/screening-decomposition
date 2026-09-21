"""Validate the gate itself, because a gate that fails everything is a broken gate.

Each target is docked into a box centred on its OWN co-crystal ligand -- the most favourable
geometry available -- so a failure is a property of the protocol, not of box placement.

The suite deliberately spans three outcomes:
  PASS            a protocol/target the gate must not reject (3PTB, the textbook redock)
  FAIL            the known-bad case the gate was built for (3N4C), plus any it finds itself
  NOT_APPLICABLE  covalent native ligands, where a non-covalent redock cannot be a control
"""
import json, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).parent))
from redock_control import redock_control, fetch_pdb, find_ligand, prepare_receptor

VINA = str(Path(__file__).parent.parent / "retrospective_benchmark/kaggle/bundle_CHEMBL612545/shard0/vina")
R = Path(__file__).parent.parent / "amr_glass/docking/results"
TARGETS = [
    ("3PTB", None,                        (22,)*3, "trypsin/benzamidine — textbook easy redock"),
    ("1STP", None,                        (22,)*3, "streptavidin/biotin — classic tight binder"),
    ("3N4C", R/"3N4C_ctss_receptor.pdbqt", (24,)*3, "cathepsin S — the case this gate was built for"),
    ("1AU0", R/"1AU0_catk_receptor.pdbqt", (24,)*3, "cathepsin K"),
    ("2YJ2", R/"2YJ2_ctsl_receptor.pdbqt", (24,)*3, "cathepsin L"),
    ("3AI8", R/"3AI8_ctsb_receptor.pdbqt", (24,)*3, "cathepsin B"),
]
rows = []
for pdb, rec, size, note in TARGETS:
    wd = Path(f"/tmp/gate_{pdb.lower()}"); wd.mkdir(exist_ok=True)
    p = fetch_pdb(pdb, wd)
    rec = str(rec) if rec and Path(rec).exists() else str(prepare_receptor(p, wd/f"{pdb}_receptor.pdbqt"))
    f = find_ligand(p)
    if not f:
        print(f"{pdb:5s} {note}: no co-crystal ligand — gate not applicable"); continue
    resn, ch, atoms, _ = f
    cen = np.mean([a[1] for a in atoms], axis=0)
    v = redock_control(rec, pdb, cen, size, exhaustiveness=16, vina_bin=VINA, workdir=str(wd))
    rows.append({"pdb": pdb, "note": note, "verdict": v})
    cov = v.get("covalent", {}).get("is_covalent")
    print(f"{pdb:5s} {resn:4s} ({len(atoms):2d} heavy){' COVALENT' if cov else '         '} "
          f"-> {v['verdict']:14s} best {str(v.get('best_rmsd')):>5s} A @ rank {v.get('best_rank')}, "
          f"top {str(v.get('top_pose_rmsd')):>5s} A   [{note}]")
    for r in v["reasons"]:
        print(f"        - {r}")
    sys.stdout.flush()

json.dump(rows, open(Path(__file__).parent/"gate_validation.json", "w"), indent=1)
from collections import Counter
c = Counter(r["verdict"]["verdict"] for r in rows)
print(f"\n{dict(c)}")
print("Gate is usable only if PASS >= 1 (it can pass) AND FAIL >= 1 (it can fail):",
      "YES" if c.get("PASS",0) and c.get("FAIL",0) else "NO")
