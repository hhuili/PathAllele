#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
[[ "${1:-integrated}" == integrated ]] || { echo 'Usage: build_native.sh [integrated]' >&2; exit 2; }
: "${CONDA_PREFIX:?Activate the calling environment first}"
PATHALLELE_NATIVE_PREFIX="$CONDA_PREFIX" cargo build --release --locked \
  --features native-aligner --manifest-path "$ROOT/Cargo.toml"
