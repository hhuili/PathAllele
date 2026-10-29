# Changelog

## 0.1.0rc1

First independent research pre-release, extracted from a native graph-path calling workflow.

* Preserve primary alignment identities and physical query ownership, preventing unsupported cross-record bridges.
* Recover eligible unaligned native sequence with conservative graph-only, uniquely oriented contexts.
* Retain complete sequence representations when primitive homology is unstable; keep unresolved genotypes missing.
* Add a portable CLI, explicit input validation, dependency/version checks, hashed stage resume, output locking, runtime memory guard, source/VCF invariants and synthetic integration tests.
* Add MIT licensing, build instructions, a pinned native-header dependency, installation metadata and continuous integration.

This release is restricted to compact haploid genomes and one sequence-resolved GFA1 reference. Broader biological accuracy and scaling claims are not made.
