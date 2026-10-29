# Third-party components

PathAllele's original Python code and thin C++ API adapter are MIT licensed. This source distribution does not bundle external tool binaries or their source trees.

* **vcflib 1.0.15**, native allele alignment API: [upstream MIT license](https://github.com/vcflib/vcflib/blob/v1.0.15/LICENSE), copyright Erik Garrison and contributors. PathAllele calls the original BiWFA API rather than copying or rewriting it.
* **WFA2**, native library and build-time headers: [upstream MIT license at the pinned header commit](https://github.com/pjotrp/WFA2-lib/blob/f86a35abf89933605d81df6e67e42edc0aaf7f55/LICENSE), copyright Santiago Marco-Sola and contributors. Header download commit: `f86a35abf89933605d81df6e67e42edc0aaf7f55`; archive SHA256: `901070d8cc5ebaafb372a0b1bb95edd164a10e3db4e50d1da15d9f3116e8e73a`.
* **minimap2**, **bcftools**, **samtools**, **NumPy**, **pysam** and their transitive dependencies are separately installed dependencies and retain their own licenses. Refer to their distributions for full notices.

Downloaded development headers retain their original LICENSE. If redistributing native binaries or dependency environments, review the licenses of every included component; this source-only package does not grant rights over those external distributions.
