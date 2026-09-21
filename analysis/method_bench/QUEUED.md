# Queued: method transferability, after the gnina determinism run

Queued 19 Sep 2026. Runs on GPU when `analysis/gnina_determinism/` finishes.

## Why a harness rather than a Boltz-2 test

A one-off Boltz-2 result answers one question once. The harness answers it for **every method we
ever build or enhance**, at the cost of one score file each — which is where compounding comes
from. Boltz-2 is simply the first method through it.

## The design decision, earned from an existing discrepancy

`analysis/boltz2/TRANSFER.json` records verdict **`TRANSFERS`**, and the manuscript reports a
cold-transfer AUROC of **0.489** on Factor Xa. Both are correct: the stored verdict measures the
*increment* Boltz-2 adds over descriptors (+0.088, clears its floor), while 0.489 is the
*absolute* performance of the transferred model — chance.

For pointing a method at a new object, only the absolute matters. A verdict of `TRANSFERS` reads
as "the model transfers" while measuring "the contribution transfers". The harness therefore
reports **both**, and `DEPLOYABLE` requires the absolute to clear 0.60 on every held-out target.

## What must happen before it can run

1. **Three targets minimum.** Mpro, Factor Xa and PD-L1 all have Boltz-2 arms; two targets give
   one transfer direction each way and no replication, so `MIN_TARGETS = 3`.
2. **Extract per-compound scores** into `scores/<method>__<target>.json` as
   `{"y": [...], "score": [...], "D": [[...]]}`. The existing result JSONs hold summaries, not
   per-compound vectors, for some targets — that extraction is the real work.
3. **Each method needs a measured floor.** Vina's is 0.00817 from 10 replicates. Boltz-2's is
   0.0201 from `TRANSFER.json`. gnina's is being measured now. A method with an unmeasured floor
   must not be entered as 0.

## Self-test first, always

`selftest.py` plants four methods with known answers. It already caught a real defect: a
**pure-descriptor** method scored residual 0.548 and passed a 0.5-based "beyond properties"
check, because out-of-fold residualisation cannot fully remove descriptor signal. The null is
now **measured** per target rather than assumed to be 0.5, and the planted case is correctly
flagged `NOT_BEYOND_PROPERTIES`.

Run `selftest.py` before trusting any result from `harness.py`.
