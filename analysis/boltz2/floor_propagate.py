#!/usr/bin/env python3
"""Propagate Boltz-2's measured replicate noise onto a real panel to get a floor.

The quantity the harness verdict uses is the *residual margin* (residual AUROC minus the
descriptor null), so that is what the floor has to be the sd of -- not the raw AUROC.
Both are reported; the control below says which reproduces the published figure.

This code did not exist in the repository when the published floors were first quoted: they
were computed inline. Written now so they are re-derivable, which also means the control is a
real test of whether this implementation matches the one that produced them.

IT DOES NOT, QUITE. The control lands at residual sd 0.00407 / upper95 0.00443 against the
published 0.00440 / 0.00472 -- within 7%, close enough to identify the published figure as the
*residual* sd rather than the margin's, and not close enough to call a reproduction. The
published number therefore remains one computed by code that is not in the repository. That is
a provenance gap this file narrows and does not close.

Usage:
  python analysis/boltz2/floor_propagate.py --control
  python analysis/boltz2/floor_propagate.py --pairs analysis/boltz2/mprofloor3/boltz2_scores.json \
      --manifest <manifest.json> --target Mpro
"""
import argparse, json, re, sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "method_bench"))
from harness import residual_auroc, residual_null          # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
PANEL = ROOT / "analysis/method_bench/scores/boltz2__{}.json"
KEY = "affinity_probability_binary"


def replicate_pairs(scores_path):
    """Pair names that differ only in an _A_/_B_ arm marker."""
    s = json.loads(Path(scores_path).read_text())
    by = {}
    for name, v in s.items():
        p = v.get(KEY)
        if p is None:
            continue
        m = re.match(r"^(.*)_([AB])_(\d+)$", name)
        if not m:
            continue
        by.setdefault((m.group(1), m.group(3)), {})[m.group(2)] = float(p)
    pairs = [(d["A"], d["B"]) for d in by.values() if "A" in d and "B" in d]
    return np.array(pairs, float)


def dist_from_bound(p):
    return np.minimum(p, 1.0 - p)


def noise_models(pairs):
    """Per-observation deviation and a level-aware proportionality constant."""
    a, b = pairs[:, 0], pairs[:, 1]
    diff = a - b
    dev = diff / np.sqrt(2.0)                 # per-observation sd from a paired difference
    sd_obs = float(np.std(diff, ddof=1) / np.sqrt(2.0))
    d = dist_from_bound((a + b) / 2.0)
    k = float(np.sqrt(np.mean(dev ** 2) / np.mean(d ** 2)))   # sd_i = k * dist_i
    return dev, sd_obs, k, d


def propagate(target, dev, k, draws, seed=20260923):
    d = json.loads(PANEL.as_posix().format(target) and Path(PANEL.as_posix().format(target)).read_text())
    y = np.asarray(d["y"], int); s = np.asarray(d["score"], float); D = np.asarray(d["D"], float)
    rng = np.random.default_rng(seed)
    dist = dist_from_bound(s)
    out = {"uniform": {"resid": [], "margin": []}, "level": {"resid": [], "margin": []}}
    for i in range(draws):
        for name, pert in (("uniform", rng.choice(dev, size=len(s), replace=True)),
                           ("level", rng.normal(0.0, k * dist))):
            sp = np.clip(s + pert, 0.0, 1.0)
            r = residual_auroc(y, sp, D)
            n = residual_null(y, sp, D)
            out[name]["resid"].append(r)
            out[name]["margin"].append(r - n)
    from scipy.stats import chi2
    res = {}
    for name in out:
        res[name] = {}
        for q in ("resid", "margin"):
            v = np.asarray(out[name][q], float)
            sd = float(v.std(ddof=1)); df = len(v) - 1
            res[name][q] = dict(sd=round(sd, 5),
                                upper95=round(float(sd * np.sqrt(df / chi2.ppf(0.05, df))), 5))
    return res, len(y)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pairs"); ap.add_argument("--target")
    ap.add_argument("--panel-ref", type=float, default=None,
                    help="panel mean distance-from-bound, for the representativeness gate")
    ap.add_argument("--draws", type=int, default=200)
    ap.add_argument("--control", action="store_true")
    ap.add_argument("--out")
    a = ap.parse_args()

    if a.control:
        # Published: 100 FXa pairs -> floor sd 0.00440, 95% upper bound 0.00472.
        pairs = replicate_pairs(ROOT / "analysis/boltz2/floor100/boltz2_scores.json")
        dev, sd_obs, k, _ = noise_models(pairs)
        print(f"CONTROL  {len(pairs)} FXa pairs, per-compound sd {sd_obs:.5f}, level k {k:.5f}")
        res, n = propagate("FXa", dev, k, a.draws)
        for name in ("uniform", "level"):
            for q in ("resid", "margin"):
                r = res[name][q]
                print(f"  {name:8} {q:7} sd {r['sd']:.5f}  upper95 {r['upper95']:.5f}")
        print("  published: sd 0.00440, upper95 0.00472 -- the row matching is the definition in use")
        return

    pairs = replicate_pairs(a.pairs)
    dev, sd_obs, k, dsample = noise_models(pairs)
    print(f"{len(pairs)} replicate pairs, per-compound sd {sd_obs:.5f}, level k {k:.5f}")

    # --- pre-registered representativeness gate: can REJECT the measurement ---
    if a.panel_ref is not None:
        ratio = float(np.mean(dsample)) / a.panel_ref
        print(f"representativeness: sample mean distance-from-bound {np.mean(dsample):.4f} "
              f"vs panel {a.panel_ref:.4f}, ratio {ratio:.3f}")
        if not 0.9 <= ratio <= 1.1:
            print("  *** REJECTED: ratio outside the pre-registered [0.9, 1.1] ***")
            print("  (PREREGISTRATION_FLOOR_MPRO2.md: rejected regardless of the value)")
            sys.exit(2)
        print("  within [0.9, 1.1] -- measurement accepted")

    res, n = propagate(a.target, dev, k, a.draws)
    print(f"\npropagated onto the real {n}-compound {a.target} panel, {a.draws} draws:")
    for name in ("uniform", "level"):
        for q in ("resid", "margin"):
            r = res[name][q]
            print(f"  {name:8} {q:7} sd {r['sd']:.5f}  upper95 {r['upper95']:.5f}")
    if a.out:
        Path(a.out).write_text(json.dumps(
            dict(target=a.target, n_pairs=len(pairs), per_compound_sd=round(sd_obs, 5),
                 level_k=round(k, 5), panel_n=n, draws=a.draws, floors=res), indent=1))
        print(f"-> {a.out}")


if __name__ == "__main__":
    main()
