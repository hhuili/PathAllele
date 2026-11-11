# Changes

## 0.2.0-dev.20261111

- Remove the standalone manual page; keep usage in CLI help and README.
- Define sequence units from complete source alleles before projection. Local reverse-complement evidence constrains graph and sequence cuts, including partial alignment rescue, so short forward matches cannot fragment supported inversion components.
- Keep source orientation evidence distinct from sequence classification; sequence-only support remains `NATIVE=SEQUENCE`.
- Add a PGGB-derived inversion regression with equivalent fragmented, unshared and cyclic representations, plus controls for remote copies, short reversals, multiple inversions and nonunique sequence evidence.


## 0.2.0 development (2026-11-10)

- Use `pathallele -g graph.gfa -r ref.fa -o out` directly; remove `call`,
  `validate`, `export`, `doctor` and the legacy adapter build command.
- Select all non-reference samples by default; `-s` remains an optional subset/order.
- Remove the duplicate full-run audit and audit-only checkpoints,
  `validation.json`, unused fields and repeated replay/hash computations.
- Strip release symbols. Keep source-sequence consistency and missingness rules.
- Work checkpoint format is now v2; start a fresh run after upgrading.

## 0.1.1 development (2026-11-10)

- Embed the existing vcflib/WFA2 adapter with the `native-aligner` feature;
  preserve process isolation through an internal worker.
- Remove `--aligner`, external-adapter discovery and environment overrides.
  Calling always uses the package-owned engine; scientific rules are unchanged.
- Add Bioconda recipe preparation and installed-package smoke tests.

## 0.1.0 (2026-11-03)

First public Rust release. Source-resolved haploid calling, complete-locus and
structural-group VCFs, explicit missingness, source validation and resumable runs.
