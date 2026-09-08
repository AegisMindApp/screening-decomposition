# The protein representation covers one chain; FABind+ indexes two

**Found 8 September 2026, in the v13 full-panel run. No FlashBind score has been read.**

## What happened

The 751-compound run cleared every stage until inference, then:

```
Error processing mpro_CHEMBL3348861: index 298 is out of bounds for axis 0 with size 298
Error processing mpro_CHEMBL151:     index 298 is out of bounds for axis 0 with size 298
KeyError: 'edge_index'
```

`bundle/mpro_receptor.pdb` is the Mpro **homodimer** — chain A 298 residues, chain B 299,
**597 total**. `bundle/prots.json` holds chain A's sequence only, so the ESM3 representation is
**[298, 1536]**. FABind+ parses the whole file and its pocket indices run 0–596. Any pocket
touching chain B indexes past the end of the representation.

This is a receptor-identity error on my side, the same class as the MSH3 chain selection
recorded in the pipeline pre-flight notes: both sides of the experiment must be the same
object, and here they were not.

## The probe passed by luck, not by correctness

All 48 probe compounds happened to select pockets entirely inside chain A, so no index exceeded
297 and the run completed. **A 48/48 success and a `coverage 1.000` were produced by a
misaligned setup.** The stratified probe was built to measure dropout, and it cannot see a
defect that only fires on compounds it does not contain. Probe scores are not to be reported.

## What the other arms used

| arm | receptor |
|---|---|
| six docking scores (Vina, gnina) | dimer — `receptor.pdbqt`, 4 chain segments, ~591 residues |
| Boltz-2 | dimer — `id: [A, B]`, two copies of the same 298-residue sequence |
| FlashBind as run | **dimer structure, monomer representation** — inconsistent |

Both comparators use the biological dimer. Mpro's active site sits at the dimer interface and
the N-finger of one protomer completes the S1 pocket of the other, so this is not a formality.

## What survives

`pocket_agreement` is unaffected. It computed centroids from the same full-model CA ordering
FABind+ indexes into, so its arithmetic was consistent throughout: median 7.88–7.90 Å across
four runs, `frac_in_box` 1.000. The docking-box control stands.

## The options

1. **Per-chain ESM3, concatenated in structure order → [597, 1536].** Same encoder, same call
   as `src/affinity/data/repr/esm3.py`, run once per chain and stacked. Keeps the dimer, matches
   both comparators, and invents no peptide bond. Deviates from FlashBind's `prot_id -> sequence`
   format only because a homodimer cannot be written in it.
2. **Concatenate the two sequences and encode as one 597-residue chain.** Simplest, but ESM3
   then sees a protein that does not exist and the junction embeddings are fiction.
3. **Restrict the receptor to chain A.** Everything aligns at 298, but FlashBind would then be
   scoring a monomer while docking and Boltz-2 scored the dimer — a receptor difference that
   sits directly on the comparison the arm exists to make.

Option 1 is the recommendation. Options 2 and 3 each introduce a defect on the axis being
measured.
