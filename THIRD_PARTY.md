# Third-party components

PathAllele's original Rust code and thin C++ API adapter are MIT licensed.
The `native-aligner` feature compiles the adapter and dynamically links the
separately installed vcflib/WFA2 libraries. The dependency packages also supply
all native headers; PathAllele does not vendor their source or headers.

* **vcflib 1.0.15**, native allele alignment API: [upstream MIT license](https://github.com/vcflib/vcflib/blob/v1.0.15/LICENSE), copyright Erik Garrison and contributors.
* **WFA2**, alignment library: [upstream MIT license](https://github.com/pjotrp/WFA2-lib/blob/f86a35abf89933605d81df6e67e42edc0aaf7f55/LICENSE), copyright Santiago Marco-Sola and contributors. The Linux vcflib package supplies its matching WFA2 library and development headers.
* **minimap2**, **bcftools**, **samtools** and their transitive dependencies are separately installed dependencies and retain their own licenses. Refer to their distributions for full notices.

Rust dependencies are pinned in `Cargo.lock`; Conda builds collect their licenses
with `cargo-bundle-licenses`. No Python runtime or Python FFI is used by PathAllele.

## Native component license notices

MIT License

Copyright (c) 2012-2022 Erik Garrison (vcflib)
Copyright (c) 2017 Santiago Marco-Sola (WFA2)

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

AUTHOR/CONTACT: Santiago Marco-Sola <santiagomsola@gmail.com>
