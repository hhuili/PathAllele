# Changes

## 0.2.0 development (2026-11-10)

- Remove the standalone man page; keep usage and output semantics in `--help`
  and README.
- Conda build 1 uses vcflib build 2 development files; remove bundled WFA2
  headers and the HTSlib include-path workaround.

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
