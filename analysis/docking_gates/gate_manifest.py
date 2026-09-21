"""Precomputed redocking-gate verdicts, and the lookup the docking pipeline calls.

WHY A MANIFEST RATHER THAN A LIVE CHECK

The gate docks a ligand and takes minutes; it also needs network access to RCSB. Neither is
acceptable inside a validator that runs per-discovery on Render. So verdicts are computed offline
by `build_manifest()` and committed, and the runtime lookup is a dict read.

The cost is that the manifest goes stale if a receptor file changes. `receptor_md5` is recorded
and checked on lookup: a receptor edited after its gate ran reads as UNKNOWN, not as its old
verdict. A stale PASS is the one failure mode that would reintroduce exactly the problem the gate
exists to prevent.

FAIL-CLOSED

`gate_status()` returns UNKNOWN for anything not in the manifest, and the caller must treat
UNKNOWN as "not permitted to report a numerical binding result". That is deliberate. The lesson
from CTSS is that an unvalidated docking assay produces confident numbers that mean nothing, so
the default for an unmeasured target is refusal, not permission.
"""
from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path
from typing import Any, Dict, Optional

HERE = Path(__file__).parent
MANIFEST = HERE / "gate_manifest.json"
_PDB_IN_NAME = re.compile(r"(?:^|[_/])([0-9][A-Za-z0-9]{3})(?:[_.]|$)")


def md5(path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def pdb_id_from_path(path) -> Optional[str]:
    """Infer the source PDB ID from a receptor filename (3RXX_receptor.pdbqt -> 3RXX).

    Inferred rather than required in the registry because RECEPTOR_REGISTRY is documented as
    freely editable, and a gate that silently skips entries missing a hand-added field would be
    fail-open.
    """
    m = _PDB_IN_NAME.search(Path(path).name)
    return m.group(1).upper() if m else None


def load() -> Dict[str, Any]:
    if not MANIFEST.exists():
        return {}
    try:
        return json.loads(MANIFEST.read_text())
    except Exception:                                                    # noqa: BLE001
        return {}


def gate_status(receptor_path) -> Dict[str, Any]:
    """Verdict for a receptor: PASS / FAIL / NOT_APPLICABLE / UNKNOWN.

    UNKNOWN covers three cases the caller must treat identically: never gated, receptor file
    changed since it was gated, and manifest missing.
    """
    p = Path(receptor_path)
    entry = load().get(p.name)
    if not entry:
        return {"verdict": "UNKNOWN",
                "reason": f"no redocking control on record for {p.name}. Run "
                          f"analysis/docking_gates/build_gate_manifest.py to measure it."}
    if not p.exists():
        return {"verdict": "UNKNOWN", "reason": f"{p.name} is missing"}
    if entry.get("receptor_md5") and entry["receptor_md5"] != md5(p):
        return {"verdict": "UNKNOWN",
                "reason": f"{p.name} has changed since its redocking control ran — the recorded "
                          f"{entry.get('verdict')} no longer applies. Re-run the manifest builder."}
    return entry


def permits_numerical_result(receptor_path) -> tuple[bool, Dict[str, Any]]:
    """Is a numerical binding claim admissible for this receptor? Only on a PASS."""
    st = gate_status(receptor_path)
    return st.get("verdict") == "PASS", st
