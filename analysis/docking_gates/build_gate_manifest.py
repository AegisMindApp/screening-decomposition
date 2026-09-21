"""Measure the redocking control for every receptor the pipeline can dock against.

Run offline; commit the resulting gate_manifest.json. Receptors are taken from
docking_prep.RECEPTOR_REGISTRY so the manifest cannot drift from what the pipeline uses.

The box used here is the receptor's own co-crystal ligand centroid -- the most favourable
geometry available. A target that fails under those conditions cannot be rescued by box tuning,
which is what makes a FAIL here decisive rather than provisional.
"""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).parent))

from redock_control import (redock_control, fetch_pdb, find_ligand,      # noqa: E402
                            covalent_to_protein)
from gate_manifest import MANIFEST, md5, pdb_id_from_path, load          # noqa: E402

VINA = str(ROOT / "analysis/retrospective_benchmark/kaggle/bundle_CHEMBL612545/shard0/vina")


def receptors():
    from analysis.dispatch.docking_prep import RECEPTOR_REGISTRY
    seen = {}
    for key, entry in RECEPTOR_REGISTRY.items():
        p = Path(entry["path"])
        seen.setdefault(p.name, {"path": p, "keywords": []})["keywords"].append(key)
    return seen


def main():
    # MERGE, DO NOT OVERWRITE. This loop only measures receptors in RECEPTOR_REGISTRY, but the
    # manifest also holds entries from other sources -- the retrospective-benchmark receptor set
    # contributed CHEMBL4523582_7VU6_receptor.pdbqt, which was the manifest's ONLY PASS. Writing
    # `out` wholesale deleted it, stripping the deployed gate of its sole positive control and
    # leaving a manifest that could only refuse. That failure was self-concealing: validate_gate.py
    # builds its own targets and never reads the manifest, so it would still report
    # "PASS >= 1 AND FAIL >= 1: YES" while nothing deployed could pass.
    #
    # Preserving an unmeasured entry cannot go stale-open, because gate_status() re-checks
    # receptor_md5 and file existence on EVERY lookup: a preserved verdict is honoured only while
    # the receptor is byte-identical to when its gate ran, and reads UNKNOWN otherwise. So the
    # merge inherits the fail-closed property rather than working around it.
    #
    # The mirror risk is silent preservation -- an entry quietly surviving is as invisible as one
    # quietly deleted -- so every preserved entry is printed, and any entry whose verdict CHANGED
    # under re-measurement is called out loudly rather than being folded into the summary.
    prior = load()
    out = {}
    for name, info in sorted(receptors().items()):
        p = info["path"]
        if not p.exists():
            print(f"{name}: receptor file missing, skipped"); continue
        pdb_id = pdb_id_from_path(p)
        if not pdb_id:
            out[name] = {"verdict": "UNKNOWN", "reason": "no PDB id inferable from filename",
                         "receptor_md5": md5(p), "keywords": info["keywords"]}
            print(f"{name}: no PDB id in filename"); continue

        wd = Path(f"/tmp/gatebuild_{pdb_id.lower()}"); wd.mkdir(exist_ok=True)
        try:
            pdb = fetch_pdb(pdb_id, wd)
            f = find_ligand(pdb)
        except Exception as e:                                           # noqa: BLE001
            out[name] = {"verdict": "UNKNOWN", "reason": f"structure fetch failed: {e}",
                         "receptor_md5": md5(p), "keywords": info["keywords"]}
            print(f"{name}: fetch failed {e}"); continue

        if not f:
            # No co-crystal ligand means no redocking control is possible. That is NOT a pass --
            # a blind box on an apo structure is precisely the MSH3 failure mode.
            out[name] = {"verdict": "NOT_APPLICABLE", "pdb_id": pdb_id,
                         "reason": "no co-crystal ligand in the structure, so the protocol cannot "
                                   "be validated against a known answer; a blind box on an apo "
                                   "structure is unvalidatable by construction",
                         "receptor_md5": md5(p), "keywords": info["keywords"]}
            print(f"{name} ({pdb_id}): NOT_APPLICABLE — apo / no ligand"); continue

        resn, ch, atoms, _ = f
        cen = np.mean([a[1] for a in atoms], axis=0)
        v = redock_control(str(p), pdb_id, cen, (24, 24, 24), exhaustiveness=16,
                           vina_bin=VINA, workdir=str(wd), seed=42)
        out[name] = {"verdict": v["verdict"], "pdb_id": pdb_id, "ligand": resn,
                     "top_pose_rmsd": v.get("top_pose_rmsd"), "best_rmsd": v.get("best_rmsd"),
                     "best_rank": v.get("best_rank"), "failure_mode": v.get("failure_mode"),
                     "covalent": v.get("covalent", {}).get("is_covalent"),
                     "rmsd_method": v.get("rmsd_method", "CCD template (exact)"),
                     "reasons": v["reasons"], "receptor_md5": md5(p),
                     "keywords": info["keywords"]}
        print(f"{name} ({pdb_id}/{resn}): {v['verdict']} — top pose {v.get('top_pose_rmsd')} A")
        sys.stdout.flush()

    changed = [(n, prior[n].get("verdict"), out[n].get("verdict")) for n in out
               if n in prior and prior[n].get("verdict") != out[n].get("verdict")]
    preserved = {n: e for n, e in prior.items() if n not in out}

    merged = {**prior, **out}
    MANIFEST.write_text(json.dumps(merged, indent=1) + "\n")

    from collections import Counter
    print(f"\nre-measured {len(out)}: {dict(Counter(e['verdict'] for e in out.values()))}")
    if preserved:
        # Printed, not summarised. These verdicts were not re-measured this run, so the reader has
        # to be able to see exactly what is being carried on trust -- including whether the file
        # backing it still exists and still hashes the same.
        print(f"preserved {len(preserved)} entry(s) not in RECEPTOR_REGISTRY:")
        for n, e in sorted(preserved.items()):
            # Scoped to analysis/ deliberately: every receptor set lives there, and an unbounded
            # rglob from ROOT crawls node_modules/ and .next/ for no benefit.
            hits = list((ROOT / "analysis").rglob(n))
            if not hits:
                state = "receptor file NOT FOUND — reads UNKNOWN at lookup"
            elif e.get("receptor_md5") and e["receptor_md5"] != md5(hits[0]):
                state = "receptor CHANGED since gated — reads UNKNOWN at lookup"
            else:
                state = "receptor unchanged — verdict still live"
            print(f"    {n}: {e.get('verdict')} ({state})")
    if changed:
        print("\n*** VERDICT CHANGED under re-measurement — investigate before committing ***")
        for n, was, now in changed:
            print(f"    {n}: {was} -> {now}")
    print(f"\nwrote {MANIFEST}: {dict(Counter(e['verdict'] for e in merged.values()))}")


if __name__ == "__main__":
    main()
