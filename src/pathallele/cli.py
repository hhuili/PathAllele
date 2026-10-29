"""Portable command line, provenance and resumable stage execution."""

import argparse
import contextlib
import fcntl
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
from . import __version__


def digest(path):
    value = hashlib.sha256()
    with open(path, "rb") as stream:
        for chunk in iter(lambda: stream.read(2**20), b""):
            value.update(chunk)
    return value.hexdigest()


def save_json(path, value):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2) + "\n")
    temporary.replace(path)


def identities(paths, root=None):
    return {str(p.relative_to(root) if root else p): digest(p) for p in sorted(paths)}


def runtime(aligner):
    result = {}
    for name in ("minimap2", "bcftools", "samtools"):
        found = shutil.which(name)
        if not found:
            raise ValueError(
                f"{name} is absent from PATH; activate the documented calling environment"
            )
        path = Path(found).resolve()
        run = subprocess.run([str(path), "--version"], text=True, capture_output=True, check=True)
        result[name] = dict(
            path=str(path), sha256=digest(path), version=(run.stdout + run.stderr).splitlines()[0]
        )
    path = Path(aligner).resolve()
    if not path.is_file() or not os.access(path, os.X_OK):
        raise ValueError(
            "An executable native adapter is required; see pathallele build-aligner --help"
        )
    result["align_alleles"] = dict(path=str(path), sha256=digest(path))
    import numpy, pysam

    result["python"] = sys.version.split()[0]
    result["numpy"] = numpy.__version__
    result["pysam"] = pysam.__version__
    return result


def checked_pipeline(commands, log):
    processes = []
    try:
        previous = None
        for i, command in enumerate(commands):
            process = subprocess.Popen(
                command,
                stdin=previous,
                stdout=subprocess.PIPE if i + 1 < len(commands) else log,
                stderr=log,
            )
            if previous is not None:
                previous.close()
            processes.append(process)
            previous = process.stdout
        codes = [process.wait() for process in processes]
        if any(codes):
            raise RuntimeError("VCF serialization pipeline failed: " + str(codes))
    finally:
        for process in processes:
            if process.poll() is None:
                process.terminate()
                process.wait()


def call(args):
    from .graph_io import Graph
    from .sources import build_sources
    from .validation import validate

    if (
        min(args.threads, args.jobs, args.minimum_anchor, args.mapping_k, args.maximum_graph_gap)
        <= 0
    ):
        raise ValueError("Thread counts and anchor/gap lengths must be positive")
    if not 0 <= args.minimum_mapq <= 60 or args.minimum_anchor < 19 or args.mapping_k < 19:
        raise ValueError("MAPQ must be 0..60; supported sequence anchors are at least 19 bp")
    if sys.flags.optimize:
        raise ValueError("Run without Python -O; scientific invariant checks must remain enabled")
    inputs = {key: Path(getattr(args, key)).resolve() for key in ("graph", "reference", "samples")}
    if args.path_map:
        inputs["path_map"] = Path(args.path_map).resolve()
    tools = runtime(args.aligner)
    graph = Graph(
        inputs["graph"],
        inputs["reference"],
        inputs["samples"],
        args.reference_path,
        inputs.get("path_map"),
    )
    estimate = len(graph.reference) * len(graph.samples) * 24 / 2**30
    if estimate > args.max_matrix_gb:
        raise ValueError(
            f"Estimated dense projection arrays need {estimate:.2f} GiB, exceeding --max-matrix-gb; this version targets compact haploid genomes"
        )
    package = Path(__file__).parent
    settings = {
        key: getattr(args, key)
        for key in (
            "threads",
            "jobs",
            "minimum_anchor",
            "mapping_k",
            "maximum_graph_gap",
            "minimum_mapq",
            "reference_path",
            "max_matrix_gb",
        )
    }
    signature = dict(
        version=__version__,
        inputs={key: dict(path=str(p), sha256=digest(p)) for key, p in inputs.items()},
        settings=settings,
        tools=tools,
        code=identities(package.glob("*.py"), package),
    )
    output = args.output.resolve()
    manifest_path = output / "run.json"
    existing = None
    if output.exists() and any(output.iterdir()):
        if not args.resume or not manifest_path.is_file():
            raise ValueError(
                "Output is nonempty; use a fresh directory or --resume for the same recorded run"
            )
        existing = json.loads(manifest_path.read_text())
        if existing["signature"] != signature:
            raise ValueError("Resume refused: input, code, settings or tool identity changed")
    output.mkdir(parents=True, exist_ok=True)
    with (output / ".lock").open("a+") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise ValueError("Another process owns this output directory") from exc
        manifest = existing or dict(
            signature=signature,
            stages={},
            started=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        )
        if manifest.get("status") == "complete":
            verify_manifest(output, manifest)
            print("Verified complete run: " + str(output))
            return
        manifest["status"] = "running"
        save_json(manifest_path, manifest)
        work = output / "work"
        work.mkdir(exist_ok=True)
        try:
            with (
                (output / "run.log").open("a", buffering=1) as log,
                contextlib.redirect_stdout(log),
                contextlib.redirect_stderr(log),
            ):
                print(
                    f"PathAllele {__version__}; {len(graph.samples)} samples; estimated projection arrays {estimate:.2f} GiB",
                    flush=True,
                )

                def stage_done(name):
                    stage = manifest["stages"].get(name)
                    if not stage:
                        return False
                    for rel, sha in stage.items():
                        p = output / rel
                        if not p.is_file() or digest(p) != sha:
                            raise ValueError("Completed stage output changed: " + rel)
                    return True

                def finish_stage(name, files):
                    manifest["stages"][name] = identities(files, output)
                    save_json(manifest_path, manifest)

                if not stage_done("sources"):
                    (work / "reference.fa").write_text(
                        ">" + graph.chrom + "\n" + graph.reference + "\n"
                    )
                    subprocess.run(
                        ["samtools", "faidx", str(work / "reference.fa")], stderr=log, check=True
                    )
                    build_sources(
                        graph,
                        work,
                        output,
                        args.threads,
                        args.minimum_anchor,
                        args.mapping_k,
                        args.maximum_graph_gap,
                        args.minimum_mapq,
                    )
                    finish_stage(
                        "sources",
                        [
                            work / n
                            for n in (
                                "reference.fa",
                                "reference.fa.fai",
                                "parents.pkl",
                                "coverage.pkl",
                                "paths.fa",
                                "context.paf",
                            )
                        ]
                        + [output / n for n in ("Graph.paths.tsv.gz", "summary.context.json")],
                    )
                if not stage_done("projection"):
                    alignment = work / "projection"
                    if alignment.exists():
                        shutil.rmtree(alignment)
                    subprocess.run(
                        [
                            sys.executable,
                            "-m",
                            "pathallele.project_variants",
                            "--blocks",
                            str(work / "parents.pkl"),
                            "--coverage",
                            str(work / "coverage.pkl"),
                            "--reference",
                            str(work / "reference.fa"),
                            "--output",
                            str(output),
                            "--work",
                            str(alignment),
                            "--jobs",
                            str(args.jobs),
                            "--aligner",
                            tools["align_alleles"]["path"],
                        ],
                        stdout=log,
                        stderr=log,
                        check=True,
                    )
                    finish_stage(
                        "projection",
                        [
                            output / n
                            for n in (
                                "Graph.vcf",
                                "Graph.alleles.tsv.gz",
                                "Graph.anchors.tsv.gz",
                                "Graph.alignment_ambiguity.tsv.gz",
                                "summary.projection.json",
                            )
                        ]
                        + [alignment / n for n in ("aligned.pkl", "projection.pkl")],
                    )
                if not stage_done("vcf"):
                    checked_pipeline(
                        [
                            [
                                "bcftools",
                                "norm",
                                "-N",
                                "-f",
                                str(work / "reference.fa"),
                                "-c",
                                "e",
                                str(output / "Graph.vcf"),
                                "-Ou",
                            ],
                            ["bcftools", "sort", "-T", str(work / "sort.XXXXXX"), "-Ou"],
                            ["bcftools", "+fill-tags", "-Ou", "--", "-t", "AC,AN,AF,NS"],
                            [
                                "bcftools",
                                "view",
                                "-c",
                                "1",
                                "-Oz",
                                "-o",
                                str(output / "Graph.vcf.gz"),
                            ],
                        ],
                        log,
                    )
                    subprocess.run(
                        ["bcftools", "index", "-ft", str(output / "Graph.vcf.gz")],
                        stdout=log,
                        stderr=log,
                        check=True,
                    )
                    finish_stage("vcf", [output / "Graph.vcf.gz", output / "Graph.vcf.gz.tbi"])
                validation = validate(output, graph)
                save_json(output / "validation.json", validation)
                manifest.update(
                    status="complete",
                    finished=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                    validation=validation,
                    outputs=identities(
                        [
                            p
                            for p in output.iterdir()
                            if p.is_file() and p.name not in ("run.json", "run.log", ".lock")
                        ],
                        output,
                    ),
                )
                save_json(manifest_path, manifest)
        except BaseException as exc:
            manifest.update(status="failed", error=type(exc).__name__ + ": " + str(exc))
            save_json(manifest_path, manifest)
            raise
    print("Completed and validated: " + str(output))


def verify_manifest(output, manifest):
    for rel, sha in manifest.get("outputs", {}).items():
        if not (output / rel).is_file() or digest(output / rel) != sha:
            raise ValueError("Output identity differs: " + rel)
    for stage in manifest.get("stages", {}).values():
        for rel, sha in stage.items():
            if not (output / rel).is_file() or digest(output / rel) != sha:
                raise ValueError("Stage identity differs: " + rel)


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="PathAllele: source-resolved haploid graph-path variant calling"
    )
    parser.add_argument("--version", action="version", version=__version__)
    commands = parser.add_subparsers(dest="command", required=True)
    run = commands.add_parser(
        "call", help="Call a fresh cohort or resume an identical interrupted run"
    )
    for key in ("graph", "reference", "samples", "output"):
        run.add_argument("--" + key, type=Path, required=True)
    run.add_argument("--path-map", type=Path)
    run.add_argument("--reference-path")
    run.add_argument("--aligner", required=True)
    run.add_argument("--threads", type=int, default=16)
    run.add_argument("--jobs", type=int, default=8)
    run.add_argument("--minimum-anchor", type=int, default=31)
    run.add_argument("--mapping-k", type=int, default=19)
    run.add_argument("--maximum-graph-gap", type=int, default=1000)
    run.add_argument("--minimum-mapq", type=int, default=20)
    run.add_argument("--max-matrix-gb", type=float, default=8)
    run.add_argument("--resume", action="store_true")
    doctor = commands.add_parser("doctor", help="Check executable dependencies and record versions")
    doctor.add_argument("--aligner", required=True)
    check = commands.add_parser(
        "validate", help="Recheck hashes and source/VCF invariants of a trusted PathAllele run"
    )
    check.add_argument("output", type=Path)
    build = commands.add_parser(
        "build-aligner", help="Build the native adapter against installed vcflib and WFA2"
    )
    build.add_argument("--prefix", type=Path, required=True, help="Dependency environment prefix")
    build.add_argument(
        "--wfa-headers",
        type=Path,
        required=True,
        help="WFA2 source tree containing wavefront/wfa.hpp",
    )
    build.add_argument("--output", type=Path, required=True)
    build.add_argument("--cxx", default=os.environ.get("CXX", "c++"))
    args = parser.parse_args(argv)
    try:
        if args.command == "call":
            call(args)
        elif args.command == "doctor":
            print(json.dumps(runtime(args.aligner), indent=2))
        elif args.command == "build-aligner":
            from .build import build_aligner

            build_aligner(args)
        else:
            from .graph_io import Graph
            from .validation import validate

            output = args.output.resolve()
            manifest = json.loads((output / "run.json").read_text())
            if manifest.get("status") != "complete":
                raise ValueError("Run is not complete")
            verify_manifest(output, manifest)
            inputs = manifest["signature"]["inputs"]
            for item in inputs.values():
                if digest(item["path"]) != item["sha256"]:
                    raise ValueError("Original input changed: " + item["path"])
            graph = Graph(
                inputs["graph"]["path"],
                inputs["reference"]["path"],
                inputs["samples"]["path"],
                manifest["signature"]["settings"]["reference_path"],
                inputs.get("path_map", {}).get("path"),
            )
            print(json.dumps(validate(output, graph), indent=2))
    except (ValueError, OSError, RuntimeError, subprocess.CalledProcessError) as exc:
        parser.exit(2, "PathAllele: " + str(exc) + "\n")


if __name__ == "__main__":
    main()
