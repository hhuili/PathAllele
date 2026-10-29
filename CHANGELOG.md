# Changelog

## 0.2.0rc1

- Define complete allele units from resolved physical graph context and unique internal sequence anchors before classifying edits. Equal nonempty lengths are MNV; unequal nonempty lengths are replacement. Size and graph orientation remain separate dimensions.
- Save every assessable observed sequence and its joint genotype in compressed catalogue sidecars. Keep the main VCF compact, with sequence-verified component projection and explicit missingness.
- Retain stable local small-variant projections in source alleles crossing unresolved territory. Complete units are never recovered from a target-site list or external callset.
- Bound normalization by source flanks, close touching units, and verify joint physical sequence replay to prevent independent repeat-copy edits from collapsing.
- Validate no overlapping or duplicate called ALT edits, complete catalogue genotypes, final VCF replay, and the dependency/code/input identities used by each run.

## 0.1.0rc2

* Define per-ALT edit `TYPE` separately from the 50-bp `SIZECLASS` boundary, including large equal-length MNVs and unequal replacements.
* Keep hierarchy unassessed at header level; clarify that `PID` denotes source contexts, not nested parent variants.
* Remove repeated constants and derivable metadata; keep only the compressed/indexed final VCF at output root, with raw serialization under `work/`.
* Validate classification for every ALT and test mixed multiallelic ordering, padding and size boundaries. Calling, allele sequence and genotype rules are unchanged.

## 0.1.0rc1

First independent research pre-release, extracted from a native graph-path calling workflow.

* Preserve primary alignment identities and physical query ownership, preventing unsupported cross-record bridges.
* Recover eligible unaligned native sequence with conservative graph-only, uniquely oriented contexts.
* Retain complete sequence representations when primitive homology is unstable; keep unresolved genotypes missing.
* Add a portable CLI, explicit input validation, dependency/version checks, hashed stage resume, output locking, runtime memory guard, source/VCF invariants and synthetic integration tests.
* Add MIT licensing, build instructions, a pinned native-header dependency, installation metadata and continuous integration.

This release is restricted to compact haploid genomes and one sequence-resolved GFA1 reference. Broader biological accuracy and scaling claims are not made.
