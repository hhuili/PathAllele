"""Extract auditable, nonoverlapping native source intervals from graph paths."""

import csv
import gzip
import json
import multiprocessing as mp
from contextlib import contextmanager
import pickle
import subprocess
from bisect import bisect_left, bisect_right
from collections import Counter
from .components import plan_path, read_alignments
from .path_context import bounded_anchors, mapped_copy_intervals


def _initialize_planner(graph, context, settings):
    global _GRAPH, _CONTEXT, _SETTINGS
    _GRAPH, _CONTEXT, _SETTINGS = graph, context, settings


def _plan(name):
    return plan_path(_GRAPH, name, _CONTEXT[name], *_SETTINGS)


@contextmanager
def planned_paths(graph, context, settings, threads):
    """Fork read-only graph state; emit results in deterministic input order."""
    names = [name for sample in graph.samples for name in graph.paths[sample]]
    if threads == 1:
        yield iter(plan_path(graph, name, context[name], *settings) for name in names)
    else:
        with mp.get_context("fork").Pool(
            min(threads, len(names)),
            initializer=_initialize_planner,
            initargs=(graph, context, settings),
        ) as pool:
            yield pool.imap(_plan, names, chunksize=1)


def build_sources(
    graph,
    work,
    output,
    threads=16,
    minimum_anchor=31,
    mapping_k=19,
    maximum_graph_gap=1000,
    minimum_mapq=20,
):
    reference = graph.reference
    nodes, rw, rp = graph.nodes, graph.refwalk, graph.ref_offsets
    with (work / "paths.fa").open("w") as stream:
        for sample in graph.samples:
            for name in graph.paths[sample]:
                stream.write(">" + name + "\n" + graph.spell(name) + "\n")
    with (work / "context.paf").open("w") as out, (work / "context.log").open("w") as log:
        subprocess.run(
            [
                "minimap2",
                "-cx",
                "asm20",
                "-k",
                str(mapping_k),
                "-s",
                "40",
                "--end-bonus",
                "100",
                "--no-hash-name",
                "--eqx",
                "--secondary=yes",
                "-N",
                "20",
                "-t",
                str(threads),
                str(work / "reference.fa"),
                str(work / "paths.fa"),
            ],
            stdout=out,
            stderr=log,
            check=True,
        )
    context = read_alignments(work / "context.paf", minimum_mapq)
    path_samples = {name: sample for sample, paths in graph.paths.items() for name in paths}
    covered = [[] for _ in graph.samples]
    chains = [[] for _ in graph.samples]
    uncertain = mapped_copy_intervals(
        work / "context.paf", graph.samples, len(reference), minimum_mapq, path_samples
    )
    windows = {}
    stats = Counter(competing_physical_copy_bp=sum(e - s for iv in uncertain for s, e in iv))
    with (
        planned_paths(
            graph, context, (minimum_anchor, mapping_k, maximum_graph_gap), threads
        ) as planned,
        gzip.open(output / "Graph.paths.tsv.gz", "wt") as stream,
    ):
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "sample",
                "path",
                "query_start0",
                "query_end0",
                "ref_start0",
                "ref_end0",
                "block_id",
                "allele_index",
                "start_step0",
                "end_step0",
                "status",
                "strand",
                "component",
                "context_kind",
            ]
        )
        for si, sample in enumerate(graph.samples):
            for name in graph.paths[sample]:
                blocks, missing, oriented = next(planned)
                length = len(oriented["+"][2])
                stats["paths"] += 1
                stats["native_bp"] += length
                for lo, hi in missing:
                    qp = oriented["+"][1]
                    writer.writerow(
                        [
                            sample,
                            name,
                            lo,
                            hi,
                            ".",
                            ".",
                            ".",
                            -1,
                            max(0, bisect_right(qp, lo) - 1),
                            bisect_left(qp, hi),
                            "UNBOUNDED",
                            "+",
                            ".",
                            "unresolved",
                        ]
                    )
                    stats["unplaced_bp"] += hi - lo
                for block in blocks:
                    strand = block["strand"]
                    runs = block["runs"]
                    walk, qp, query = oriented[strand]
                    stats[block["kind"] + "_components"] += 1

                    def emit(qs, qe, rs, re, pid, ai, status):
                        start = max(0, bisect_right(qp, qs) - 1)
                        end = bisect_left(qp, qe)
                        if strand == "-":
                            qs, qe = length - qe, length - qs
                            start, end = len(walk) - end, len(walk) - start
                        writer.writerow(
                            [
                                sample,
                                name,
                                qs,
                                qe,
                                rs,
                                re,
                                pid,
                                ai,
                                start,
                                end,
                                status,
                                strand,
                                block["identity"],
                                block["kind"],
                            ]
                        )

                    for r, q, n, payload in runs:
                        if reference[r : r + n] != query[q : q + n]:
                            raise ValueError("Source match differs from reference: " + name)
                        covered[si].append((r, r + n))
                        stats["context_matched_bp"] += n
                        emit(q, q + n, r, r + n, ".", 0, "MATCH")
                    first, last = runs[0], runs[-1]
                    chains[si].append((first[0], last[0] + last[2]))
                    if block["record"] is not None:
                        index = block["record"].index
                    else:
                        # Graph-only internal homology is constrained to the same
                        # bounded graph context, never borrowed from another PAF.
                        raw = [(q, q + n, r) for r, q, n, p in runs]
                        index = (raw, [v[0] for v in raw], [v[1] for v in raw])
                    for left, right in zip(runs, runs[1:]):
                        rs, re = left[0] + left[2], right[0]
                        qs, qe = left[1] + left[2], right[1]
                        if rs > re or qs > qe:
                            raise ValueError("Source component is not collinear: " + name)
                        covered[si].append((rs, re))
                        if reference[rs:re] == query[qs:qe]:
                            if qe > qs:
                                emit(qs, qe, rs, re, ".", 0, "MATCH")
                                stats["gap_matched_bp"] += qe - qs
                            continue
                        ref, alt = reference[rs - 1 : re + 1], query[qs - 1 : qe + 1]
                        if not ref or not alt or ref[0] != alt[0] or ref[-1] != alt[-1]:
                            raise ValueError("Source allele lacks two exact flanks: " + name)
                        ra, rb = bisect_left(rp, rs), bisect_right(rp, re) - 1
                        qa, qb = bisect_left(qp, qs), bisect_right(qp, qe) - 1
                        ri, qi = list(rw[ra:rb]), list(walk[qa:qb])
                        rnodes = set(ri)
                        back = sum(
                            len(nodes[abs(n)]) for n in qi if -n in rnodes and n not in rnodes
                        )
                        forward = sum(len(nodes[abs(n)]) for n in qi if n in rnodes)
                        complete = rp[ra] == rs and rp[rb] == re and qp[qa] == qs and qp[qb] == qe
                        kind = (
                            "INV"
                            if complete and ri and qi == [-n for n in ri[::-1]]
                            else "COMPLEX"
                            if back >= minimum_anchor and back > forward
                            else "SEQUENCE"
                        )
                        anchors = (
                            []
                            if kind != "SEQUENCE"
                            else bounded_anchors(
                                ref,
                                alt,
                                ri,
                                qi,
                                nodes,
                                rp[ra] - rs + 1,
                                qp[qa] - qs + 1,
                                index,
                                rs - 1,
                                qs - 1,
                            )
                        )
                        signature = tuple((r, q, n) for r, q, n, node in anchors)
                        key = (rs - 1, re + 1, alt, signature, kind)
                        if key not in windows:
                            ident = f"PC{len(windows):06d}"
                            windows[key] = dict(
                                chrom=graph.chrom,
                                pos=rs,
                                gid=ident,
                                pid=ident,
                                alleles=[ref, alt],
                                carriers=[[], []],
                                anchors=[[], anchors],
                                native_kind=["REFERENCE", kind],
                            )
                        parent = windows[key]
                        parent["carriers"][1].append(si)
                        emit(
                            qs,
                            qe,
                            rs,
                            re,
                            parent["gid"],
                            1,
                            "PARTIAL_SEQUENCE" if set(alt) - set("ACGT") else "ALLELE",
                        )
                        stats["allele_bp"] += qe - qs
                        stats["allele_observations"] += 1
                        stats["native_" + kind] += 1
                        stats[block["kind"] + "_allele_observations"] += 1
            if (si + 1) % 100 == 0:
                print(f"Physical path context: {si + 1}/{len(graph.samples)}", flush=True)
    parents = list(windows.values())
    for parent in parents:
        parent["carriers"] = [sorted(set(v)) for v in parent["carriers"]]
    with (work / "parents.pkl").open("wb") as stream:
        pickle.dump((graph.samples, parents), stream, protocol=5)
    with (work / "coverage.pkl").open("wb") as stream:
        pickle.dump((covered, uncertain, chains), stream, protocol=5)
    if stats["native_bp"] != sum(
        stats[k] for k in ("context_matched_bp", "gap_matched_bp", "allele_bp", "unplaced_bp")
    ):
        raise ValueError("Native path bases do not partition exactly into source states")
    stats["blocks"] = len(parents)
    (output / "summary.context.json").write_text(json.dumps(dict(stats), indent=2) + "\n")
    print(json.dumps(dict(stats), indent=2), flush=True)
