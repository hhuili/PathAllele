#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
: "${CONDA_PREFIX:?Activate the calling environment first}"
BUILD="$ROOT/tools/build"
mkdir -p "$BUILD/wfa-headers"
ARCHIVE="$BUILD/wfa.tar.gz"
curl -fL --retry 2 'https://codeload.github.com/pjotrp/WFA2-lib/tar.gz/f86a35abf89933605d81df6e67e42edc0aaf7f55' -o "$ARCHIVE"
EXPECTED=901070d8cc5ebaafb372a0b1bb95edd164a10e3db4e50d1da15d9f3116e8e73a
if command -v sha256sum >/dev/null 2>&1; then
  ACTUAL=$(sha256sum "$ARCHIVE")
else
  ACTUAL=$(shasum -a 256 "$ARCHIVE")
fi
[[ "${ACTUAL%% *}" == "$EXPECTED" ]] || { echo 'WFA2 header archive checksum mismatch' >&2; exit 1; }

tar xzf "$ARCHIVE" --strip-components=1 -C "$BUILD/wfa-headers"
pathallele build-aligner --prefix "$CONDA_PREFIX" --wfa-headers "$BUILD/wfa-headers" \
  --output "$ROOT/tools/pathallele-align" --cxx "${CXX:-c++}"
