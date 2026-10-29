"""Independent checks of complete joint observations and compact VCF replay."""

from bisect import bisect_right
from collections import Counter, defaultdict
import csv
import gzip
import numpy as np

from .allele_model import Edit, intersects, replay


def load_catalogue(output, samples, reference):
    with np.load(output / "Graph.locus_genotypes.npz", allow_pickle=False) as archive:
        data = {key: archive[key] for key in archive.files}
    assert list(data["samples"]) == samples
    alleles = defaultdict(list)
    bounds = {ident: tuple(map(int, b)) for ident, b in zip(data["loci"], data["bounds"])}
    with gzip.open(output / "Graph.loci.tsv.gz", "rt") as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            ident = row["locus_id"]
            assert bounds[ident] == (int(row["ref_start0"]), int(row["ref_end0"]))
            assert int(row["allele_index"]) == len(alleles[ident])
            alleles[ident].append(row["sequence"])
    assert len(set(data["loci"])) == len(data["loci"]) == len(alleles)
    for ident, seqs in alleles.items():
        s, e = bounds[ident]
        assert seqs[0] == reference[s:e]
        assert len(seqs) == len(set(seqs))
        assert all(not (set(seq) - set("ACGT")) for seq in seqs)
    return data, alleles, {ident: i for i, ident in enumerate(data["loci"])}


def validate_catalogue(data, alleles, models, emitted, callable_mask, reference, chains):
    counts = Counter()
    bins = []
    concrete = np.broadcast_to(
        np.frombuffer(reference.encode(), dtype=np.uint8), callable_mask.shape
    ).copy()
    chain_index = []
    for si, model in enumerate(models):
        index = defaultdict(set)
        for e in model:
            for k in range(e.start // 1000, e.end // 1000 + 1):
                index[k].add(e)
            concrete[si, e.start : e.end] = (
                ord("?") if e.structural else np.frombuffer(e.alt.encode(), dtype=np.uint8)
            )
        bins.append(index)
        intervals = sorted(chains[si])
        ends = []
        end = -1
        for s, e in intervals:
            end = max(end, e)
            ends.append(end)
        chain_index.append(([s for s, e in intervals], ends))
    for li, ident in enumerate(data["loci"]):
        s, e = map(int, data["bounds"][li])
        seqs = alleles[ident]
        gt = data["gt"][li]
        projected = data["projected_gt"][li]
        valid = gt >= 0
        assert np.array_equal(valid, data["source_status"][li] == 0)
        assert np.array_equal(projected >= 0, data["export_status"][li] == 0)
        assert np.all(valid[projected >= 0])
        assert set(range(1, len(seqs))).issubset(set(gt[valid]))
        counts["joint_resolved_genotypes"] += int(valid.sum())
        counts["joint_alleles"] += len(seqs) - 1
        counts["joint_loci"] += 1
        if data["kinds"][li] == "SNP":
            expected = np.array([ord(seq) for seq in seqs], dtype=np.uint8)
            assert np.array_equal(concrete[valid, s], expected[gt[valid]]), ident
            assert np.all(callable_mask[valid, s]), ident
            continue
        for si in np.flatnonzero(valid):
            lo, hi = max(0, s - 1), min(len(reference), e + 1)
            starts, ends = chain_index[si]
            ci = bisect_right(starts, lo) - 1
            assert (
                ci >= 0
                and ends[ci] >= hi
                and np.all(callable_mask[si, s:e] if s < e else callable_mask[si, lo:hi])
            ), (ident, si)
            found = set().union(
                *(bins[si].get(k, set()) for k in range(lo // 1000, hi // 1000 + 1))
            )
            local = [x for x in found if intersects(Edit(s, e, ""), x) and (s != e or x.structural)]
            assert all(s <= x.start <= x.end <= e for x in local), (ident, si)
            observed = replay(
                reference[s:e], [Edit(x.start - s, x.end - s, x.alt, x.complex) for x in local]
            )
            assert observed == seqs[gt[si]], (ident, si, "joint sequence")
            if projected[si] >= 0:
                actual = [x for x in local if (x.start, x.end, x.alt) in emitted[si]]
                represented = replay(
                    reference[s:e], [Edit(x.start - s, x.end - s, x.alt, x.complex) for x in actual]
                )
                assert represented == observed, (ident, si, "VCF component replay")
                counts["joint_complete_VCF_replays"] += 1
                if gt[si] > 0 and projected[si] == 0:
                    counts["reference_backbone_with_emitted_components"] += 1
    return counts
