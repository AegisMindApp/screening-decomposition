# Pre-registration: replicating the interventions on a second target

Written 24 September 2026, **before any Factor Xa intervention run**. The manuscript measures
five interventions on one target. The first question an examiner asks is whether they replicate,
and the manuscript currently has no answer. This fixes what will be run, what will be read from
it, and — equally important — **which of the five cannot honestly be replicated and why**.

## Which of the five can be replicated, and which cannot

| intervention | Mpro result | replicable on Factor Xa? | why |
|---|---|---|---|
| receptor preparation repair | +0.045 | **yes, CPU** | both receptors already exist. `bundle_CHEMBL244`'s `receptor_md5` is `3c58477e…`, which is exactly the `receptor_md5_before_protonation` recorded in `bundle_CHEMBL244_PROT`. Same compounds, same box, one variable |
| eight-fold search effort | +0.008 | **yes, CPU** | exhaustiveness 4 is already run at ten seeds; exhaustiveness 32 is a re-run of the same bundle |
| pose ensemble vs top pose | +0.055 | **no** | checked before promising it: the bundle docks `num_modes: 5` but the stored output keeps **only the top-mode affinity** (`{"CHEMBL100732": -10.78, …}`). Replicating it needs a fresh run that retains modes, which is a larger job than the two below |
| scoring function (gnina CNN) | +0.140 | **not now** | needs gnina on Factor Xa poses — GPU, and a pose set that is not yet stored |
| binding site to a co-folding model | −0.013 | **not now** | needs a blind Boltz-2 Factor Xa arm, ~9 GPU-hours, and the Mpro blind arm is still running |

**Three of five will remain single-target and the paper must say so**, rather than describing
itself as replicated. A partial replication reported as a replication is the failure this
document exists to prevent.

## Arm 1 — receptor preparation repair

Dock the compounds common to `bundle_CHEMBL244` (unprotonated, 894) and
`bundle_CHEMBL244_PROT` (protonated, 886) against each receptor, exhaustiveness 4, three seeds
per arm, everything else identical. ΔAUROC = protonated − unprotonated, paired within seed.

This is the *cleanest* of the five to replicate: the two bundles differ in one documented step
(`obabel -xr -p 7.4`) and the md5s prove it.

## Arm 2 — eight-fold search effort

Exhaustiveness 32 on `bundle_CHEMBL244_PROT`, logically sharded (`names[shard::n_shards]`, so no
new bundle), paired within seed against the existing exhaustiveness-4 seeds. Six seeds.

Six rather than ten is a declared compromise: ten would be ~120 kernel-runs. Six is enough to
detect an effect the size Mpro's was, and **the Mpro experience — the verdict moved between six
and ten seeds — is a reason to treat any six-seed Factor Xa verdict as provisional**, which is
stated here rather than discovered later.

## Reading rules, fixed now

Each arm is charged the floor for **its own target, protocol and panel**. Factor Xa's
exhaustiveness-4 floor is measured (sd 0.00427 on 886, 95% bound 0.00703); its exhaustiveness-32
floor will be measured from the same run and must not be borrowed from Mpro's.

- **REPLICATES** — the Factor Xa effect has the same sign as Mpro's and its 95% CI overlaps
  Mpro's CI. The intervention is a property of the method, not of the target.
- **DOES NOT REPLICATE** — same sign but non-overlapping intervals, or opposite sign with the
  Factor Xa interval clearing its floor. The manuscript's single-target results stand as measured
  and must be described as target-specific.
- **INCONCLUSIVE** — the Factor Xa interval straddles its floor. The replication is
  underpowered and says nothing either way; it does **not** count as support.

A replication that returns INCONCLUSIVE is reported as prominently as one that succeeds. The
whole point of adding a second target is to find out, and an underpowered null that gets quietly
dropped would make the exercise worse than not running it.

## What a successful replication would and would not establish

It would show the effect is not a peculiarity of one protein. It would **not** establish that the
intervention is useful: Vina is below chance on Mpro (0.4175) and at 0.6712 on Factor Xa, so an
intervention can replicate and still leave a screen not worth running. The manuscript's SURVIVES
verdict is about resolvability throughout, and adding a target does not change that.

Nor does it make the *floor* results more general than they already are — those are separately
measured on both targets (§3.2) and do not depend on this.

---

## Amendment, 24 September 2026 — scope widened to four of five, before any of the new runs

Two things changed the day this was written, and both are recorded here **before the runs they
authorise**.

**1. Pose ensemble and gnina come into scope.** The table above rules both out. That was wrong
about the cause in the first case and pessimistic in the second: `kaggle_shard_worker.py` already
passes `--out {name}_out.pdbqt` and Vina writes **all five modes** into it. The Mpro pose set
(`exh4_poses/`, 753 files) is exactly that output, kept. Nothing about the protocol needs
changing — the seed-floor path simply discards the directory.

So **one Factor Xa docking run at exhaustiveness 4 with the output retained** supplies both:

- **pose ensemble** — parse the five modes per compound, compare the ensemble score against the
  top mode, exactly as on Mpro;
- **gnina CNN rescoring** — `--score_only` over those same poses, so the pose noise is shared
  between arms, which is the condition that made the Mpro comparison a paired one.

The replication therefore covers **four of five interventions**. Only the Boltz-2 binding-site
row remains single-target, and it is the row whose verdict is INCONCLUSIVE under every floor
tried, so it is the least informative one to add.

**2. Arm 2 moves to GCP, and arm 1 stays on Kaggle.** Kaggle's GPU quota is exhausted until
Saturday and its five CPU slots would take ~2.5 days for arm 2's 72 shards. The same
`run_seeds.py` runs on a GCP instance against a copy of the same bundle. This is a change of
*where*, not of *what*: the runner asserts the receptor and vina md5s against the manifest on
every host, so a cross-platform difference cannot pass silently. A prior cross-platform check on
this pipeline gave r = 0.9842 with mean |Δ| 0.101 kcal/mol between local and Kaggle docking
(`pose_ensemble/PROTONATED_RESULT.md`), which is the magnitude to expect and is recorded now so
it is not discovered as a surprise.

### Reading rules for the two new arms, fixed now

Both are paired comparisons on the same compounds, so both take **Factor Xa's own
exhaustiveness-4 floor** (sd 0.00427 on 886, 95% bound 0.00703), rescaled to whatever panel
survives pose parsing. Neither takes Mpro's floor and neither takes the exhaustiveness-32 floor.

- **pose ensemble** — Mpro gave +0.055 [+0.035, +0.075]. REPLICATES if the Factor Xa effect is
  positive and its interval overlaps that one.
- **gnina** — Mpro gave +0.140 [+0.089, +0.191], on a panel where Vina sits below chance.
  **Factor Xa is at 0.6712, above chance**, so there is less room for a large gain and a smaller
  effect here would not be a failure to replicate in the interesting sense. That asymmetry is
  stated now because it is exactly the kind of thing that is easy to rationalise afterwards: a
  smaller gnina effect on Factor Xa is **expected**, and the pre-registered test is therefore
  the **sign and the clearing of the floor**, not the overlap of intervals.

### What would make this amendment dishonest

Widening scope after seeing results. Nothing from any Factor Xa intervention run exists yet: arm
1 has not started (Kaggle's CPU slots are still occupied by Mpro seeds 11-15) and arm 2 has not
been pushed. The only Factor Xa data in hand is the exhaustiveness-4 seed floor, which is an
input to all of this and was published before any of it was contemplated.
