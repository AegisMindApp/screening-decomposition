#!/usr/bin/env python3
"""Method-agnostic transferability harness: can a scoring method be pointed at a NEW target?

Built 19 Sep 2026. The point is not Boltz-2. The point is that every future method — a new
co-folding model, a rescorer, an enhancement to one we already have — gets evaluated the same
way, so a method that survives can be applied to objects and compound.

ADDING A METHOD costs one file. Drop `scores/<method>__<target>.json` holding
{compound_id: score}, list it in `registry.json`, and every metric below is computed for it.

WHAT IT MEASURES, and why each exists
-------------------------------------
1. **In-target AUROC.** The number everyone quotes. Necessary, not sufficient.
2. **Residual AUROC** — descriptors regressed out. Answers "is this more than molecular
   properties?" Molecular weight alone explains ~half of a docking score, so a method can look
   good while ranking by size.
3. **Cold transfer, reported TWO ways.** This is the harness's reason to exist:
     * `cold_absolute`  — fit on another target, apply here. What you would actually get.
     * `cold_increment` — the same, minus the descriptor-only cold baseline.
   A method can have a large increment and useless absolute performance. Measured on
   Boltz-2 (`analysis/boltz2/TRANSFER.json`): Mpro->FXa increment +0.088, clearing its floor,
   while the absolute cold AUROC is **0.489** — chance. The stored verdict reads "TRANSFERS",
   which is true of the increment and false of the model. Reporting only one of these is how a
   method gets deployed on a target it cannot rank.
4. **Every comparison is charged against a reproducibility floor**, passed in per method. A
   difference smaller than the method's own run-to-run noise is not a difference.

VERDICTS
--------
* `DEPLOYABLE`    — cold_absolute clears `usable_auroc` on EVERY held-out target. Can be
                    pointed at a new object.
* `CONTRIBUTES`   — cold_increment clears the floor but cold_absolute does not. Useful as a
                    component where the target has its own labels; NOT deployable cold.
* `IN_TARGET_ONLY`— in-target only; nothing survives transfer.
* `NOT_BEYOND_PROPERTIES` — residual AUROC does not clear chance + floor. It is ranking by size.
* `INSUFFICIENT`  — fewer than `min_targets` targets. Two targets cannot show transferability.

    python analysis/method_bench/harness.py --registry analysis/method_bench/registry.json
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy import stats

HERE = Path(__file__).resolve().parent
MIN_TARGETS = 3          # two targets give one transfer direction each way and no replication
USABLE_AUROC = 0.60      # below this a ranking is not worth acting on


def auroc(y, s):
    y = np.asarray(y, int)
    s = np.asarray(s, float)
    pos, neg = s[y == 1], s[y == 0]
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    r = stats.rankdata(np.concatenate([pos, neg]))
    return float((r[:len(pos)].sum() - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg)))


def residual_auroc(y, s, D, seed=0):
    """AUROC of the method score after regressing out descriptors, out-of-fold.

    Out-of-fold matters: fitting the descriptor model on the same rows it is removed from
    absorbs method signal into the descriptor fit and understates the residual.
    """
    from sklearn.ensemble import GradientBoostingRegressor
    from sklearn.model_selection import KFold
    s = np.asarray(s, float)
    D = np.asarray(D, float)
    resid = np.zeros_like(s)
    for tr, te in KFold(5, shuffle=True, random_state=seed).split(D):
        m = GradientBoostingRegressor(random_state=seed).fit(D[tr], s[tr])
        resid[te] = s[te] - m.predict(D[te])
    return auroc(y, resid)


def residual_null(y, s, D, seed=0):
    """The residual AUROC a method would show if it contained NOTHING but descriptors.

    The null is NOT 0.5. Residualisation is imperfect out-of-fold -- a gradient-boosted fit
    cannot fully remove descriptor signal -- so a pure-property score retains apparent residual
    signal. Measured on a planted pure-descriptor method: 0.548, which cleared a 0.5-based
    threshold and was wrongly called 'beyond properties' until this null was added.

    Construction: take the out-of-fold DESCRIPTOR PREDICTION of the method's own score. That is
    by definition a pure function of descriptors with the same scale and shape as the method,
    and its residual is what imperfect residualisation leaves behind.
    """
    from sklearn.ensemble import GradientBoostingRegressor
    from sklearn.model_selection import KFold
    s = np.asarray(s, float)
    D = np.asarray(D, float)
    pred = np.zeros_like(s)
    for tr, te in KFold(5, shuffle=True, random_state=seed + 1).split(D):
        m = GradientBoostingRegressor(random_state=seed).fit(D[tr], s[tr])
        pred[te] = m.predict(D[te])
    return residual_auroc(y, pred, D, seed=seed)


def cold_transfer(donor, recipient, seed=0):
    """Fit on the donor target in full, apply cold to the recipient. Returns (with, without).

    `with` uses descriptors + the method score; `without` uses descriptors alone. The pair is
    what separates 'the method contributes' from 'the transferred model is usable'.
    """
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler

    def fit_apply(cols):
        Xd = np.column_stack([donor[c] for c in cols])
        Xr = np.column_stack([recipient[c] for c in cols])
        sc = StandardScaler().fit(Xd)
        m = LogisticRegression(max_iter=2000).fit(sc.transform(Xd), donor["y"])
        return auroc(recipient["y"], m.decision_function(sc.transform(Xr)))

    return fit_apply(["D", "score"]), fit_apply(["D"])


MAX_DESCRIPTOR_BASELINE = 0.80   # above this a target cannot discriminate method from properties


def target_admissible(y, D, seed=0):
    """Can this target distinguish a method from molecular properties at all?

    If seven free descriptors already separate the panel, nothing measured on it can show a
    method is 'beyond properties' -- the headroom does not exist. PD-L1 is the case that forced
    this: 79% active, inactives 127 Da heavier than actives, descriptor baseline 0.8740, Vina
    at 0.5661 with a CI including chance. Recorded as compromised in
    analysis/boltz2/PREREGISTRATION_TARGETS.md BEFORE it ran, and then nearly used anyway as
    the third target for a transferability verdict.

    Returns (admissible, descriptor_baseline).
    """
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import cross_val_predict
    from sklearn.preprocessing import StandardScaler
    from sklearn.pipeline import make_pipeline
    y = np.asarray(y, int)
    D = np.asarray(D, float)
    pipe = make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000))
    pred = cross_val_predict(pipe, D, y, cv=5, method="decision_function")
    base = auroc(y, pred)
    return bool(base <= MAX_DESCRIPTOR_BASELINE), base


def evaluate(method, targets, floor, usable=USABLE_AUROC, min_targets=MIN_TARGETS):
    targets = dict(targets)
    out = {"method": method, "floor": floor, "n_targets_supplied": len(targets),
           "in_target": {}, "residual": {}, "cold": {}}
    # Drop inadmissible targets BEFORE anything is computed on them, and say which.
    out["inadmissible"] = {}
    for t in list(targets):
        ok, base = target_admissible(targets[t]["y"], targets[t]["D"])
        out.setdefault("descriptor_baseline", {})[t] = base
        if not ok:
            out["inadmissible"][t] = (
                f"descriptor baseline {base:.4f} > {MAX_DESCRIPTOR_BASELINE}; properties "
                f"already separate this panel, so nothing measured here can show a method "
                f"exceeds them")
            del targets[t]

    for t, d in targets.items():
        out["in_target"][t] = auroc(d["y"], d["score"])
        out["residual"][t] = residual_auroc(d["y"], d["score"], d["D"])
        out.setdefault("residual_null", {})[t] = residual_null(d["y"], d["score"], d["D"])

    names = sorted(targets)
    for donor in names:
        for rec in names:
            if donor == rec:
                continue
            w, wo = cold_transfer(targets[donor], targets[rec])
            out["cold"][f"{donor}->{rec}"] = {
                "cold_absolute": w, "descriptors_only": wo, "cold_increment": w - wo,
                "absolute_usable": bool(w >= usable),
                "increment_clears_floor": bool((w - wo) > floor),
            }

    res = np.array(list(out["residual"].values()), float)
    # Clear the MEASURED null for each target, not 0.5.
    nulls = np.array([out["residual_null"][t] for t in out["residual"]], float)
    margins = res - nulls
    cold = out["cold"].values()
    abs_ok = all(c["absolute_usable"] for c in cold) if cold else False
    inc_ok = all(c["increment_clears_floor"] for c in cold) if cold else False
    beyond = bool(np.nanmin(margins) > floor) if len(margins) else False

    out["n_targets"] = len(targets)
    if len(targets) < min_targets:
        v = (f"INSUFFICIENT — {len(targets)} targets; {min_targets} needed before "
             f"transferability means anything")
    elif not beyond:
        v = (f"NOT_BEYOND_PROPERTIES — worst residual margin {np.nanmin(margins):+.4f} over "
             f"the measured pure-descriptor null does not clear the floor ({floor}); "
             f"it is ranking by molecular properties")
    elif abs_ok:
        v = ("DEPLOYABLE — cold absolute AUROC clears "
             f"{usable} on every held-out target; can be pointed at a new object")
    elif inc_ok:
        v = ("CONTRIBUTES, NOT DEPLOYABLE — the increment clears the floor everywhere but the "
             "cold absolute does not. Useful only where the target has its own labels.")
    else:
        v = "IN_TARGET_ONLY — nothing survives transfer"
    out["verdict"] = v
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--registry", default=str(HERE / "registry.json"))
    ap.add_argument("--out", default=str(HERE / "method_bench.json"))
    a = ap.parse_args()
    reg = json.loads(Path(a.registry).read_text())

    results = {}
    for method, spec in reg["methods"].items():
        targets = {}
        for t, path in spec["scores"].items():
            p = Path(path)
            if not p.exists():
                continue
            d = json.loads(p.read_text())
            targets[t] = {"y": d["y"], "score": d["score"], "D": d["D"]}
        results[method] = evaluate(method, targets, spec.get("floor", 0.0))
        r = results[method]
        print(f"\n=== {method}  ({r['n_targets']} targets, floor {r['floor']})")
        for t in sorted(r["in_target"]):
            print(f"    {t:10} in-target {r['in_target'][t]:.4f}   "
                  f"residual {r['residual'][t]:.4f}  null {r['residual_null'][t]:.4f}  "
                  f"margin {r['residual'][t]-r['residual_null'][t]:+.4f}")
        for k, c in r["cold"].items():
            print(f"    {k:18} cold_absolute {c['cold_absolute']:.4f}  "
                  f"increment {c['cold_increment']:+.4f}  "
                  f"usable={c['absolute_usable']}")
        print(f"    VERDICT: {r['verdict']}")
    Path(a.out).write_text(json.dumps(results, indent=1))
    print(f"\n-> {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
