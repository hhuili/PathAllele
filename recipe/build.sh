#!/usr/bin/env bash
set -euo pipefail
export PATHALLELE_NATIVE_PREFIX="$PREFIX"
export CARGO_HOME="$SRC_DIR/cargo-home"
export CARGO_BUILD_JOBS="${CPU_COUNT:-2}"
cargo-bundle-licenses --format yaml --output THIRDPARTY.yml
cargo install --locked --no-track --root "$PREFIX" --path . --features native-aligner
install -Dm644 pathallele.1 "$PREFIX/share/man/man1/pathallele.1"
mkdir -p "$PREFIX/share/pathallele/examples/demo"
cp examples/demo/{graph.gfa,reference.fa,samples.txt,truth.tsv} "$PREFIX/share/pathallele/examples/demo/"
