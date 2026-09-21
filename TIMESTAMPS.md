# What the timestamps in this repository do and do not establish

The papers this repository supports claim that reading rules were fixed **before** the data
existed. A reader should be able to check that rather than take it on trust. This file states
exactly how far the evidence goes, because the answer differs by date.

## Two tiers of evidence

**Tier 1 — material public since 5–8 September 2026.** The repository was created
2026-09-05 and last pushed 2026-09-08 before this addition. Everything present at that point
carries a GitHub-side record that it existed by then, independent of anything the author
controls. For pre-registrations whose corresponding experiments ran after 8 September, that is
genuine third-party evidence of precedence.

**Tier 2 — material added 21 September 2026** (`analysis/seed_floor/`,
`analysis/method_bench/`, the Boltz-2 floor work, and `docs/papers/resolution_v3_draft.md`).
These artifacts were written and committed in a private repository between 18 and 21 September,
and became public on 21 September. **The public record therefore establishes only that they
existed by 21 September — not that each pre-registration predated its own experiment.**

For tier 2, precedence rests on the commit ordering in the private history. Commit dates can be
set arbitrarily by whoever makes the commit, so that ordering is author-controlled and is *not*
equivalent to third-party attestation. It is offered as a record, not as proof.

## Why the repository was not re-exported

The original export used `git filter-repo`, which preserves dates but rewrites hashes. Re-running
it would have changed every existing public hash and broken the hash map recorded in the
supplementary material. The 21 September material was therefore added as ordinary commits on top
of the existing history, so all previously published hashes remain valid.

## What would be stronger

Committing future pre-registrations directly to this public repository, before the corresponding
run starts, would put them in tier 1 by construction. That is the intended practice from here.
An external timestamping authority, or an OpenTimestamps proof on the pre-registration commits,
would be stronger still and costs nothing.
