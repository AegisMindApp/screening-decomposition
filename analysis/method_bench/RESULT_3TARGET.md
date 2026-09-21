# Boltz-2 on three admissible targets: NOT_BEYOND_PROPERTIES

21 Sep 2026. The ALDH1 panel exists so that a transferability verdict becomes *possible*; the
harness has refused one since September on the grounds that two targets give one transfer
direction each way and no replication. With ALDH1 scored, the rule fixed in `harness.py` before
these numbers existed returns:

> **NOT_BEYOND_PROPERTIES** — worst residual margin +0.0042 over the measured pure-descriptor
> null does not clear the floor (0.0201); it is ranking by molecular properties.

## The numbers

| target | in-target | residual | measured null | margin | clears 0.0201? |
|---|---|---|---|---|---|
| Mpro | 0.7910 | 0.6406 | 0.5835 | **+0.0571** | yes |
| ALDH1 | 0.5629 | 0.5818 | 0.5218 | **+0.0600** | yes |
| Factor Xa | 0.7227 | 0.5474 | 0.5431 | **+0.0042** | **no** |

Two of three clear comfortably. Factor Xa does not, and the rule requires every admissible
target to clear — one target where a method reduces to properties is enough to stop it being
pointed at a new object, because nothing tells you in advance which kind of target you have.

## The verdict's weakest link, stated plainly

**Boltz-2's floor of 0.0201 is inherited, not measured by us.** It comes from
`analysis/boltz2/TRANSFER.json` and has never been re-derived by seed replication the way
Vina's was (`analysis/seed_floor/RESULT.md`, 10 replicates). The verdict turns on whether
+0.0042 clears it, so a materially smaller floor would flip Factor Xa and with it the
conclusion.

That is exactly the substitution this project's own resubmission gate blocks on, and it is
flagged here rather than quietly relied upon. **The honest status is: NOT_BEYOND_PROPERTIES
conditional on a floor we did not measure.** Measuring Boltz-2's floor by seed replication is
now the single highest-value outstanding job on this line — it is the difference between a
verdict and a verdict-shaped object.

Note the asymmetry in what that would change. A larger true floor only strengthens the negative.
A smaller one weakens it. So the current conclusion is the *conservative* reading, not a
convenient one.

## What transfers, which is not what was expected

| direction | cold absolute | increment | usable (>=0.60)? |
|---|---|---|---|
| ALDH1 -> Mpro | **0.7271** | +0.1832 | yes |
| FXa -> Mpro | **0.6800** | +0.3247 | yes |
| ALDH1 -> FXa | **0.6653** | +0.1347 | yes |
| FXa -> ALDH1 | 0.5732 | +0.0537 | no |
| Mpro -> ALDH1 | 0.5055 | +0.0101 | no |
| Mpro -> FXa | 0.4970 | +0.0924 | no |

Three of six directions are usable. The surprise is **ALDH1**: the target where Boltz-2 is
weakest in-target (0.5629) is its *best donor*, transferring usefully to both proteases. And
ALDH1 is the hardest *recipient* — nothing transfers into it above 0.58.

That asymmetry is informative rather than noise. ALDH1 was chosen precisely because its
descriptor baseline is lowest (0.5698), i.e. it has the least property signal for a model to
lean on. A model fitted there learns something less property-shaped, which travels; models
fitted on the proteases do not travel into it.

## What this closes

The generation constraint and the "point a validated method at objects" plan rested on Boltz-2
surviving this harness. It does not. Boltz-2 should **not** be pointed at new objects on the
strength of in-target AUROC, and `method_ladder.json` now records that as a measured verdict on
three admissible targets rather than an `INSUFFICIENT` awaiting data.

Vina remains `INSUFFICIENT` at two admissible targets; PD-L1 stays excluded (descriptor
baseline 0.874).

## Provenance

750/750 compounds scored, **zero failures**, 300 actives (40.0%), no overlap between the four
shards — each checked rather than assumed. Descriptor baseline on the *delivered* panel is
0.5698, identical to the figure `THIRD_TARGET.md` pre-registered on the designed panel, which
is the re-check that document committed to. Run cost ~$20 of L4 time across spot and on-demand.
