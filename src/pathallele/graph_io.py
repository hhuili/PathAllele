"""Validated GFA1 P paths and explicit haploid cohort identities."""

from array import array
from collections import defaultdict
import csv
import gzip
from pathlib import Path

COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")
IUPAC = frozenset("ACGTRYSWKMBDHVN")


def reverse_complement(sequence):
    return sequence.translate(COMPLEMENT)[::-1]


def text_open(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def read_fasta(path):
    records = {}
    name = None
    with text_open(path) as stream:
        for line in stream:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if not line[1:].strip():
                    raise ValueError("Empty FASTA header")
                name = line[1:].split()[0]
                if name in records:
                    raise ValueError(f"Duplicate FASTA record: {name}")
                records[name] = []
            elif name is None:
                raise ValueError("FASTA sequence precedes its header")
            else:
                if set(line.upper()) - IUPAC:
                    raise ValueError(f"Unsupported FASTA sequence symbols: {name}")
                records[name].append(line.upper())
    if not records or any(not value for value in records.values()):
        raise ValueError("FASTA must contain nonempty records")
    return {key: "".join(value) for key, value in records.items()}


def read_samples(path):
    rows = Path(path).read_text().splitlines()
    samples = [line.strip() for line in rows if line.strip() and not line.startswith("#")]
    if not samples or len(samples) != len(set(samples)):
        raise ValueError("Samples must be nonempty and unique")
    if any(any(c.isspace() for c in s) for s in samples):
        raise ValueError("One sample identifier is required per line")
    return samples


def offsets(nodes, walk):
    result = [0]
    for node in walk:
        result.append(result[-1] + len(nodes[abs(node)]))
    return result


class Graph:
    """Zero-overlap, sequence-resolved GFA1 paths; arbitrary segment identifiers."""

    def __init__(self, graph, reference, samples, reference_path=None, path_samples=None):
        refs = read_fasta(reference)
        if len(refs) != 1:
            raise ValueError("Exactly one reference FASTA record is supported")
        self.chrom, self.reference = next(iter(refs.items()))
        if set(self.reference) - set("ACGT"):
            raise ValueError("The reference must be completely ACGT resolved")
        self.samples = read_samples(samples)
        if self.chrom in self.samples:
            raise ValueError("The reference must not be a genotyped sample")
        explicit = None
        if path_samples:
            explicit = {}
            with open(path_samples) as stream:
                reader = csv.DictReader(stream, delimiter="\t")
                if not {"path", "sample"} <= set(reader.fieldnames or []):
                    raise ValueError("Path map needs path and sample columns")
                for row in reader:
                    if row["path"] in explicit or row["sample"] not in self.samples:
                        raise ValueError("Duplicate path or undeclared sample in path map")
                    explicit[row["path"]] = row["sample"]
        self.nodes, names, edges = {}, {}, set()

        def handle(name):
            if name not in names:
                names[name] = len(names) + 1
            return names[name]

        with text_open(graph) as stream:
            for line in stream:
                fields = line.rstrip("\n").split("\t")
                if fields[0] == "S":
                    key = handle(fields[1])
                    sequence = fields[2].upper()
                    if key in self.nodes or not sequence or set(sequence) - IUPAC:
                        raise ValueError(
                            "Duplicate, missing or invalid node sequence: " + fields[1]
                        )
                    self.nodes[key] = sequence
                elif fields[0] == "L":
                    if fields[5] not in ("0M", "*"):
                        raise ValueError("Overlapping GFA segments are not supported")
                    if fields[2] not in ("+", "-") or fields[4] not in ("+", "-"):
                        raise ValueError("Invalid GFA link orientation")
                    left = handle(fields[1]) * (1 if fields[2] == "+" else -1)
                    right = handle(fields[3]) * (1 if fields[4] == "+" else -1)
                    edges.add((left, right))
                    edges.add((-right, -left))
        self.paths = defaultdict(list)
        self.walks = {}
        reference_candidates = []
        haplotypes = defaultdict(set)
        with text_open(graph) as stream:
            for line in stream:
                if not line.startswith("P\t"):
                    continue
                _, name, tokens, overlaps, *_ = line.rstrip("\n").split("\t")
                sample = explicit.get(name) if explicit is not None else name.split("#")[0]
                is_reference = (
                    name == reference_path
                    if reference_path
                    else (name == self.chrom or name.startswith(self.chrom + "#"))
                )
                if sample not in self.samples and not is_reference:
                    continue
                if name in self.walks:
                    raise ValueError("Duplicate GFA path: " + name)
                if overlaps != "*" and any(x != "0M" for x in overlaps.split(",")):
                    raise ValueError("Nonzero path overlaps are not supported: " + name)
                walk = array("i")
                for token in tokens.split(","):
                    if len(token) < 2 or token[-1] not in "+-" or token[:-1] not in names:
                        raise ValueError("Unknown node or malformed path step: " + token)
                    key = names[token[:-1]]
                    if key not in self.nodes:
                        raise ValueError("Path node has no sequence: " + token)
                    walk.append(key if token[-1] == "+" else -key)
                if any((x, y) not in edges for x, y in zip(walk, walk[1:])):
                    raise ValueError("GFA path has an edge absent from the graph: " + name)
                self.walks[name] = walk
                if is_reference:
                    if explicit is not None and name in explicit:
                        raise ValueError(
                            "The reference path cannot be assigned to a genotyped sample"
                        )
                    reference_candidates.append(name)
                else:
                    self.paths[sample].append(name)
                    if explicit is None and "#" in name:
                        haplotypes[sample].add(name.split("#")[1])
        if len(reference_candidates) != 1:
            raise ValueError("Select one complete GFA reference path with --reference-path")
        self.reference_path = reference_candidates[0]
        if self.spell(self.reference_path) != self.reference:
            raise ValueError("Spelled graph reference differs from the reference FASTA")
        if set(self.samples) - set(self.paths):
            raise ValueError(
                "Samples have no selected graph paths: "
                + ",".join(sorted(set(self.samples) - set(self.paths)))
            )
        if any(len(haps) > 1 for haps in haplotypes.values()):
            raise ValueError("Multiple haplotypes per sample are unsupported; use a haploid cohort")
        if explicit is not None and set(explicit) - set(self.walks):
            raise ValueError("Explicit sample paths are absent from the graph")
        self.refwalk = self.walks[self.reference_path]
        self.ref_offsets = offsets(self.nodes, self.refwalk)
        self.ref_occurrences = {
            (node, self.ref_offsets[i]): i for i, node in enumerate(self.refwalk)
        }
        self.ref_places = defaultdict(list)
        for i, node in enumerate(self.refwalk):
            self.ref_places[abs(node)].append(i)

    def spell(self, name):
        return "".join(
            self.nodes[n] if n > 0 else reverse_complement(self.nodes[-n]) for n in self.walks[name]
        )
