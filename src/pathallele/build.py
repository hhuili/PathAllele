"""Build only our thin adapter; dependencies remain separately installed."""

from pathlib import Path
import json
import shutil
import subprocess
import tempfile
from .cli import digest


def build_aligner(args):
    prefix = args.prefix.resolve()
    headers = args.wfa_headers.resolve()
    output = args.output.resolve()
    if output.exists():
        raise ValueError("Refusing to overwrite an existing native adapter")
    for file in (prefix / "include/vcflib/Variant.h", headers / "wavefront/wfa.hpp"):
        if not file.is_file():
            raise ValueError("Missing development header: " + str(file))
    compiler = shutil.which(args.cxx)
    if not compiler:
        raise ValueError("C++ compiler not found: " + args.cxx)
    output.parent.mkdir(parents=True, exist_ok=True)
    source = Path(__file__).parent / "native/align_alleles.cpp"
    with tempfile.TemporaryDirectory(prefix="pathallele-build-") as temp:
        include = Path(temp) / "htslib"
        include.mkdir()
        (include / "htslib").symlink_to(prefix / "include/htslib", target_is_directory=True)
        command = [
            compiler,
            "-std=c++17",
            "-O2",
            "-I" + temp,
            "-I" + str(prefix / "include/vcflib"),
            "-I" + str(prefix / "include"),
            "-I" + str(headers),
            str(source),
            "-L" + str(prefix / "lib"),
            "-Wl,-rpath," + str(prefix / "lib"),
            "-lvcflib",
            "-lwfa2",
            "-o",
            str(output),
        ]
        try:
            subprocess.run(command, check=True)
        except BaseException:
            output.unlink(missing_ok=True)
            raise
    (output.parent / (output.name + ".build.json")).write_text(
        json.dumps(
            dict(command=command, source_sha256=digest(source), binary_sha256=digest(output)),
            indent=2,
        )
        + "\n"
    )
    print(output)
