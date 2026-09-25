#!/usr/bin/env python3
"""Two objections the manuscript does not currently answer, both from committed data.

Rules fixed in PREREGISTRATION_METRIC_SCOPE.md before this ran.

  1. Is the floor metric-specific?  Screening papers report EF@1%, not AUROC.
  2. Is the floor level-specific?   Ours was measured where Vina is BELOW CHANCE (0.4175).
"""
import json, math, sys
from pathlib import Path
import numpy as np
from scipy.stats import chi2

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from paired_exhaustiveness import load           # noqa: E402


def auroc(y, s):
    y = np.asarray(y, int); s = np.asarray(s, float)
    # Vina affinity is lower-is-better, so rank on the negated score.
    o = np.argsort(-(-s), kind="mergesort")
    r = np.empty(len(s), float); r[np.argsort(-s, kind="mergesort")] = np.arange(len(s))
    from scipy.stats import rankdata
    rk = rankdata(-s)
    na, nn = y.sum(), len(y) - y.sum()
    return (rk[y == 1].sum() - na * (na + 1) / 2) / (na * nn)


def ef_at(y, s, frac):
    """Enrichment factor at the top `frac` of the ranking (lower score = better)."""
    y = np.asarray(y, int); s = np.asarray(s, float)
    n = len(y); k = max(1, int(round(frac * n)))
    top = np.argsort(s, kind="mergesort")[:k]
    hit = y[top].mean()
    base = y.mean()
    return hit / base if base > 0 else float("nan")


def bedroc(y, s, alpha=20.0):
    y = np.asarray(y, int); s = np.asarray(s, float)
    n = len(y); na = int(y.sum()); ra = na / n
    order = np.argsort(s, kind="mergesort")
    ranks = np.where(y[order] == 1)[0] + 1
    rie_num = np.exp(-alpha * ranks / n).sum()
    rie_den = na / n * (1 - math.exp(-alpha)) / (math.exp(alpha / n) - 1)
    rie = rie_num / rie_den
    fac = ra * math.sinh(alpha / 2) / (math.cosh(alpha / 2) - math.cosh(alpha / 2 - alpha * ra))
    return rie * fac + 1 / (1 - math.exp(alpha * (1 - ra)))


def bound(v):
    v = np.asarray(v, float); sd = v.std(ddof=1); df = len(v) - 1
    return sd, sd * math.sqrt(df / chi2.ppf(0.05, df))


def part1():
    print("=" * 78)
    print("1. IS THE FLOOR METRIC-SPECIFIC?  Same seeds, four metrics.")
    print("=" * 78)
    rows = []
    for exh, d in ((4, HERE / "results"), (32, HERE / "results")):
        sc, lab = load(exh, str(d))
        seeds = sorted(sc)
        common = sorted(set.intersection(*[set(sc[s]) for s in seeds]))
        y = [lab[c] for c in common]
        print(f"\nexhaustiveness {exh}: {len(seeds)} seeds, {len(common)} compounds, "
              f"{sum(y)} active ({np.mean(y):.1%})")
        print(f"  {'metric':10} {'mean':>9} {'sd':>9} {'95% bound':>10} {'CV = sd/mean':>13}")
        for name, fn in (("AUROC", lambda yy, ss: auroc(yy, ss)),
                         ("EF@1%", lambda yy, ss: ef_at(yy, ss, 0.01)),
                         ("EF@5%", lambda yy, ss: ef_at(yy, ss, 0.05)),
                         ("BEDROC20", lambda yy, ss: bedroc(yy, ss))):
            v = np.array([fn(y, [sc[s][c] for c in common]) for s in seeds])
            sd, up = bound(v)
            cv = sd / abs(v.mean()) if v.mean() else float("nan")
            print(f"  {name:10} {v.mean():9.4f} {sd:9.5f} {up:10.5f} {cv:13.4f}")
            rows.append(dict(exh=exh, metric=name, mean=round(float(v.mean()), 5),
                             sd=round(float(sd), 6), bound95=round(float(up), 6),
                             cv=round(float(cv), 5), n=len(common), seeds=len(seeds)))
    a = next(r for r in rows if r["exh"] == 4 and r["metric"] == "AUROC")
    e = next(r for r in rows if r["exh"] == 4 and r["metric"] == "EF@1%")
    ratio = e["cv"] / a["cv"]
    print(f"\n  PRE-REGISTERED: CV(EF@1%) / CV(AUROC) > 3 at exh=4?")
    print(f"  measured {e['cv']:.4f} / {a['cv']:.4f} = {ratio:.2f}x  -> "
          f"{'CONFIRMED' if ratio > 3 else 'NOT CONFIRMED'}")
    return rows, ratio


def _panel_at(y, pool, delta, rng):
    """Assign the observed score multiset to compounds so the panel attains an AUROC set by delta.

    latent = delta*label + N(0,1); rank on latent; hand the sorted scores out in that order.
    Labels, panel size, class balance and the score multiset are all untouched -- only which
    compound holds which score changes. delta<0 gives below-chance panels, which is where this
    benchmark actually sits.
    """
    latent = delta * y + rng.standard_normal(len(y))
    order = np.argsort(-latent, kind="mergesort")     # best first
    s = np.empty(len(y))
    s[order] = pool                                    # pool ascending == best first (lower=better)
    return s


def _delta_for(y, pool, target, rng):
    lo, hi = -6.0, 6.0
    for _ in range(40):
        mid = (lo + hi) / 2
        a = np.mean([auroc(y, _panel_at(y, pool, mid, rng)) for _ in range(6)])
        if a < target:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def part2(n_draws=300, seed=20260924):
    print("\n" + "=" * 78)
    print("2. IS THE FLOOR LEVEL-SPECIFIC?  Everything fixed but the AUROC level.")
    print("=" * 78)
    sc, lab = load(4, str(HERE / "results"))
    seeds = sorted(sc)
    common = sorted(set.intersection(*[set(sc[s]) for s in seeds]))
    y = np.array([lab[c] for c in common], int)
    M = np.array([[sc[s][c] for c in common] for s in seeds], float)   # seeds x compounds
    dev = M - M.mean(axis=0)            # each compound's OWN observed seed deviations
    pool = np.sort(M[0])
    rng = np.random.default_rng(seed)
    real = float(np.mean([auroc(y, M[i]) for i in range(len(seeds))]))
    real_sd = float(np.std([auroc(y, M[i]) for i in range(len(seeds))], ddof=1))

    def induced(s):
        """Perturb by RESAMPLING each compound's own measured deviations -- heteroscedastic,
        not an assumed Gaussian."""
        v = []
        for _ in range(n_draws):
            k = rng.integers(0, dev.shape[0], dev.shape[1])
            v.append(auroc(y, s + dev[k, np.arange(dev.shape[1])]))
        return bound(v)

    print(f"n={len(y)} ({int(y.sum())} active), score multiset and labels preserved, "
          f"noise resampled from the measured per-compound deviations, {n_draws} draws\n")

    # --- POSITIVE CONTROL: at the real panel, does the simulation reproduce the real floor?
    sd_c, _ = induced(M.mean(axis=0))
    print(f"  POSITIVE CONTROL at the real panel (AUROC {real:.4f})")
    print(f"    simulated sd {sd_c:.5f}   measured across seeds {real_sd:.5f}   "
          f"ratio {sd_c/real_sd:.2f}")
    ok = 0.5 <= sd_c / real_sd <= 2.0
    print(f"    {'within 2x -- the simulation tracks reality, trend below is usable' if ok else '*** OUTSIDE 2x -- trend below is NOT usable ***'}\n")

    print(f"  {'target':>7} {'realised':>9} {'induced sd':>11} {'95% bound':>10} {'vs peak':>8}")
    out = []
    for target in (0.30, 0.42, 0.50, 0.60, 0.70, 0.80, 0.90):
        d = _delta_for(y, pool, target, rng)
        s = _panel_at(y, pool, d, rng)
        realised = auroc(y, s)
        sd, up = induced(s)
        out.append(dict(target=target, realised=round(float(realised), 4),
                        sd=round(float(sd), 6), bound95=round(float(up), 6)))
    peak = max(r["sd"] for r in out)
    for r in out:
        print(f"  {r['target']:7.2f} {r['realised']:9.4f} {r['sd']:11.5f} "
              f"{r['bound95']:10.5f} {r['sd']/peak:8.2f}")
    pk = max(out, key=lambda r: r["sd"])
    print(f"\n  PRE-REGISTERED: sd largest near 0.5, falling away from it?")
    print(f"  peak at realised AUROC {pk['realised']:.3f}, |A-0.5| = {abs(pk['realised']-0.5):.3f}")
    hi = min(out, key=lambda r: abs(r["realised"] - 0.80))
    print(f"  sd at {hi['realised']:.3f} is {hi['sd']:.5f} = {hi['sd']/peak:.2f}x the peak")
    print(f"  -> {'CONFIRMED' if abs(pk['realised']-0.5) < 0.15 and hi['sd'] < peak else 'NOT CONFIRMED'}")
    out.append(dict(control_simulated_sd=round(sd_c, 6), control_measured_sd=round(real_sd, 6),
                    control_ratio=round(sd_c / real_sd, 3), real_auroc=round(real, 4)))
    return out


if __name__ == "__main__":
    rows, ratio = part1()
    lv = part2()
    (HERE / "METRIC_SCOPE.json").write_text(json.dumps(
        dict(metrics=rows, ef_cv_ratio=round(ratio, 3), level_dependence=lv,
             preregistration="analysis/seed_floor/PREREGISTRATION_METRIC_SCOPE.md"), indent=1))
    print("\n-> analysis/seed_floor/METRIC_SCOPE.json")
