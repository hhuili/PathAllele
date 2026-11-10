# Conda recipe

Build from a **clean source checkout**, as for TeloXplorer:

```sh
conda build recipe --override-channels -c conda-forge -c hhuili -c bioconda
```

The development package uses `source: path: ..`. Use a clean checkout/staged source tree containing no analysis inputs, results, environments or build caches. For formal Bioconda submission, replace this source stanza with the published release archive URL and its real SHA256.

Rust and C++ compilers are build dependencies. `vcflib` is both the native host dependency and an ABI-pinned runtime dependency. The tested Linux vcflib package supplies `libvcflib` and `libwfa2`, so a second WFA2 runtime package is not added. vcflib 1.0.15 build 2 supplies matching WFA2 headers under `$PREFIX/include/wfa2lib` and fixes its exported development files. This corrected dependency build is distributed through `hhuili` pending upstream discussion in [vcflib issue #456](https://github.com/vcflib/vcflib/issues/456); no upstream recipe PR has been submitted. PathAllele neither vendors nor downloads WFA2 headers. Put `hhuili` before `bioconda` when testing this build, including with strict channel priority.

The recipe currently targets Linux x86-64. It is not a `noarch` package: the installed Rust executable and embedded C++ interface contain platform-specific machine code. A pure Rust executable would also require platform-specific builds. `run_test.sh` checks installed-package operation without an external adapter or obsolete option; compiler packages are not runtime dependencies. macOS/native-platform support requires separate build validation before enabling it. Conda handles runtime dependencies; users need only `conda install` and the regular PathAllele commands.

Install the Linux development package from https://anaconda.org/hhuili/pathallele:

```sh
conda create -n pathallele --override-channels -c conda-forge -c hhuili -c bioconda pathallele=0.2.0.dev20261110
conda activate pathallele
pathallele --help
```

The backend remains Rust plus the embedded vcflib/WFA2 adapter, with an internal worker for process isolation. Direct Rust-to-WFA2 FFI is a possible later refactor, not an implemented feature of this package.

The current 0.2.0 development source removes all public subcommands and makes
`-s` optional. Its CLI is `pathallele -g graph.gfa -r ref.fa -o out`.
