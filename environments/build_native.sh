#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
: "${CONDA_PREFIX:?Activate the calling environment first}"
BUILD="$ROOT/tools/build"
mkdir -p "$BUILD/wfa-headers"
ARCHIVE="$BUILD/wfa.tar.gz"
curl -fL --retry 2 'https://codeload.github.com/pjotrp/WFA2-lib/tar.gz/f86a35abf89933605d81df6e67e42edc0aaf7f55' -o "$ARCHIVE"
python - "$ARCHIVE" <<'PY'
import hashlib, sys
from pathlib import Path
expected='901070d8cc5ebaafb372a0b1bb95edd164a10e3db4e50d1da15d9f3116e8e73a'
actual=hashlib.sha256(Path(sys.argv[1]).read_bytes()).hexdigest()
if actual!=expected:raise SystemExit('WFA2 header archive checksum mismatch')
PY
tar xzf "$ARCHIVE" --strip-components=1 -C "$BUILD/wfa-headers"
pathallele build-aligner --prefix "$CONDA_PREFIX" --wfa-headers "$BUILD/wfa-headers" \
  --output "$ROOT/tools/pathallele-align" --cxx "${CXX:-c++}"
