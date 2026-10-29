"""Source sequence, physical partition and emitted genotype invariants.

This is internal consistency validation, not an independent biological truth set.
Only read work pickles created by this trusted run (never third-party pickles).
"""

from collections import Counter, defaultdict
import csv
import gzip
import pickle
import numpy as np
import pysam
from .allele_model import Edit, replay, trim_edit
from .graph_io import reverse_complement


def require(condition, message):
    if not condition:
        raise ValueError("Validation failed: " + message)


def validate(output, graph):
    work = output / "work"
    with (work / "parents.pkl").open("rb") as stream:
        samples, parents = pickle.load(stream)
    with (work / "projection/aligned.pkl").open("rb") as stream:
        aligned = pickle.load(stream)
    with (work / "projection/projection.pkl").open("rb") as stream:
        models, positive, callable_mask = pickle.load(stream)
    parent_by_id = {p["gid"]: p for p in parents}
    samples_index = {sample: i for i, sample in enumerate(samples)}
    observed = defaultdict(set)
    partitions = defaultdict(list)
    stats = Counter()
    reference = graph.reference
    current_name = None
    query = ""
    with gzip.open(output / "Graph.paths.tsv.gz", "rt") as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            name = row["path"]
            lo = int(row["query_start0"])
            hi = int(row["query_end0"])
            if name != current_name:
                query = graph.spell(name)
                current_name = name
            require(name in graph.paths[row["sample"]], "path/sample ownership")
            require(0 <= lo <= hi <= len(query), "physical query bounds")
            if lo < hi:
                partitions[name].append((lo, hi))
            stats["source_" + row["status"]] += 1
            if row["status"] == "UNBOUNDED":
                continue
            rs, re = int(row["ref_start0"]), int(row["ref_end0"])
            sequence = query if row["strand"] == "+" else reverse_complement(query)
            qs, qe = (lo, hi) if row["strand"] == "+" else (len(query) - hi, len(query) - lo)
            if row["status"] == "MATCH":
                require(sequence[qs:qe] == reference[rs:re], "source exact match")
            else:
                parent = parent_by_id[row["block_id"]]
                require(parent["alleles"][0] == reference[rs - 1 : re + 1], "source REF sequence")
                require(parent["alleles"][1] == sequence[qs - 1 : qe + 1], "source ALT sequence")
                observed[parent["gid"]].add(samples_index[row["sample"]])
    for sample in graph.samples:
        for name in graph.paths[sample]:
            intervals = sorted(partitions[name])
            cursor = 0
            for lo, hi in intervals:
                require(lo == cursor, "native bases overlap or disappear: " + name)
                cursor = hi
            require(cursor == len(graph.spell(name)), "unaccounted path end: " + name)
            stats["fully_partitioned_paths"] += 1
    for parent, edits in zip(parents, aligned):
        require(set(parent["carriers"][1]) == observed[parent["gid"]], "source carrier identity")
        offset = parent["pos"] - 1
        for allele, es in zip(parent["alleles"], edits):
            lo = min([offset] + [offset + e.start for e in es])
            hi = max([offset + len(parent["alleles"][0])] + [offset + e.end for e in es])
            expected = (
                reference[lo:offset] + allele + reference[offset + len(parent["alleles"][0]) : hi]
            )
            actual = replay(
                reference[lo:hi],
                [Edit(offset + e.start - lo, offset + e.end - lo, e.alt, e.complex) for e in es],
            )
            require(actual == expected, "normalized source allele replay: " + parent["gid"])
            stats["replayed_source_alleles"] += 1
    model_keys = [{(e.start, e.end, e.alt) for e in model} for model in models]
    emitted = [set() for _ in samples]
    ids = set()
    alleles_seen = set()
    last = -1
    with pysam.VariantFile(output / "Graph.vcf.gz") as vcf:
        require(list(vcf.header.samples) == samples, "VCF cohort order")
        require(list(vcf.header.contigs) == [graph.chrom], "VCF reference identity")
        for record in vcf:
            require(record.start >= last, "VCF sorting")
            last = record.start
            require(record.id not in ids, "duplicate VCF ID")
            ids.add(record.id)
            require(record.ref == reference[record.start : record.stop], "VCF reference sequence")
            allele_edits = {0: None}
            ac = [0] * len(record.alts)
            an = 0
            for ai, alt in enumerate(record.alts, 1):
                key = (record.start, record.ref, alt)
                require(key not in alleles_seen, "duplicate ALT")
                alleles_seen.add(key)
                require(not (set(alt) - set("ACGT")), "unresolved bases in confident catalogue")
                allele_edits[ai] = trim_edit(record.start, record.ref, alt)
            for si, sample in enumerate(samples):
                call = record.samples[sample]
                gt = call["GT"]
                ft = call.get("FT")
                require(len(gt) == 1, "non-haploid GT")
                value = gt[0]
                if value is None:
                    require(ft != "PASS", "missing GT marked PASS")
                    stats["missing_" + str(ft)] += 1
                    continue
                require(ft == "PASS", "resolved GT without PASS")
                an += 1
                require(
                    bool(np.all(callable_mask[si, record.start : record.stop])),
                    "confident GT in unresolved coverage",
                )
                if value:
                    require(1 <= value <= len(record.alts), "invalid ALT index")
                    edit = allele_edits[value]
                    require(
                        (edit.start, edit.end, edit.alt) in model_keys[si],
                        "ALT has no sequence-derived model witness",
                    )
                    emitted[si].add((edit.start, edit.end, edit.alt))
                    ac[value - 1] += 1
                    stats["source_supported_ALT_genotypes"] += 1
                else:
                    stats["reference_genotypes"] += 1
            require(
                list(record.info["AC"]) == ac and record.info["AN"] == an, "AC/AN inconsistency"
            )
            require(all(ac), "unobserved ALT in VCF")
            stats["vcf_records"] += 1
            stats["vcf_ALT_alleles"] += len(record.alts)
    for si in range(len(samples)):
        expected = {(e.start, e.end, e.alt) for e in positive[si]}
        require(emitted[si] == expected, "final VCF differs from positive projection model")
        stats["unemitted_model_edits"] += len(model_keys[si] - emitted[si])
    return dict(
        status="pass",
        samples=len(samples),
        counts=dict(stats),
        scope="Every native path base partitioned; every source allele matches its physical graph occurrence and replays after normalization; every emitted ALT has a source-derived model witness; confident genotypes require resolved coverage; REF, sample order, IDs, sort, AC/AN and projection serialization checked.",
        limitation="Unresolved source models remain missing or compound. A pass does not establish sensitivity, specificity, whole-genome haplotype completeness, or independent biological accuracy.",
    )
