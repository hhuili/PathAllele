# Third-party components

PathAllele's original Rust code and thin C++ API adapter are MIT licensed. The source includes the unmodified WFA2 interface-header include closure and its MIT license under native/wfa2_headers; alignment implementation code and external binaries are supplied by separately installed packages.

The `native-aligner` feature compiles the same adapter into PathAllele and
dynamically links the separately installed vcflib/WFA2 libraries. Their algorithms
and adapter scoring/decomposition settings are unchanged. The external adapter
is retained only as a developer regression reference, not as a calling backend.

* **vcflib 1.0.15**, native allele alignment API: [upstream MIT license](https://github.com/vcflib/vcflib/blob/v1.0.15/LICENSE), copyright Erik Garrison and contributors. PathAllele calls the original vcflib/WFA2 API rather than copying or rewriting it.
* **WFA2**, native library and build-time headers: [upstream MIT license at the pinned header commit](https://github.com/pjotrp/WFA2-lib/blob/f86a35abf89933605d81df6e67e42edc0aaf7f55/LICENSE), copyright Santiago Marco-Sola and contributors. Header download commit: `f86a35abf89933605d81df6e67e42edc0aaf7f55`; archive SHA256: `901070d8cc5ebaafb372a0b1bb95edd164a10e3db4e50d1da15d9f3116e8e73a`.
* **minimap2**, **bcftools**, **samtools** and their transitive dependencies are separately installed dependencies and retain their own licenses. Refer to their distributions for full notices.

Downloaded development headers retain their original LICENSE. If redistributing native binaries or dependency environments, review the licenses of every included component; this source-only package does not grant rights over those external distributions.

Rust dependencies are pinned in `Cargo.lock`; license metadata is available
through `cargo metadata`. No Python runtime or Python FFI is used.
