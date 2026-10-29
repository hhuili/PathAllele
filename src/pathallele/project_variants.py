#!/usr/bin/env python3
"""Align native allele intervals and project one shared haploid variant catalogue."""

import argparse
import copy
from collections import Counter, defaultdict
from bisect import bisect_right
import csv
import gzip
import json
from pathlib import Path
import pickle
import numpy as np
from . import __version__
from .graph_io import reverse_complement
from .allele_model import DNA, Edit, atomic_edits, intersects, replay, vcf_alleles
from .graph_anchors import anchored_gaps, unique_runs
from .sequence_anchors import refine_anchors
from .joint_catalogue import JointCatalogue, factor_state
from .normalize_alleles import (
    normalize,
    directional_disagreements,
    preserve_complex_gaps,
    cohere_complex_models,
    close_touching_units,
)


def decompose_alleles(parents):
    """Complete sequence units bounded before calling, without a gapped CIGAR.

    Graph identity and unique internal homology define units. Shared padding is
    trimmed; a nonempty mixed unit stays intact regardless of its best-scoring
    alignment. A one-base substitution and a pure insertion/deletion remain
    naturally simple. No coordinates, cohort labels or known variants enter.
    """
    edits = [[[] for _ in p["alleles"]] for p in parents]
    for p, ee in zip(parents, edits):
        kinds = p.get("native_kind", ["REFERENCE"] + ["SEQUENCE"] * (len(ee) - 1))
        for ai, allele in enumerate(p["alleles"][1:], 1):
            if kinds[ai] in ("INV", "COMPLEX"):
                ee[ai] = [Edit(1, len(p["alleles"][0]) - 1, allele[1:-1], True)]
            else:
                for start, ref, alt in anchored_gaps(p["alleles"][0], allele, p["anchors"][ai]):
                    ee[ai].extend(atomic_edits(start, ref, alt, keep_complex=bool(ref and alt)))
            if replay(p["alleles"][0], ee[ai]) != allele:
                raise ValueError(f"Sequence units failed to replay: {p['pid']} allele {ai}")
    return edits


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--blocks", type=Path, required=True)
    p.add_argument("--coverage", type=Path, required=True)
    p.add_argument("--reference", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--work", type=Path, required=True)
    p.add_argument("--aligner", required=True)
    p.add_argument("--jobs", type=int, default=8)
    p.add_argument("--sequence-anchor", type=int, default=31)
    args = p.parse_args(argv)
    args.work.mkdir(parents=True, exist_ok=False)
    samples, parents = pickle.load(args.blocks.open("rb"))
    coverage, uncertain, chains = pickle.load(args.coverage.open("rb"))
    reference = "".join(
        (l for l in args.reference.read_text().splitlines() if not l.startswith(">"))
    )
    n = len(samples)
    length = len(reference)
    known = np.zeros((n, length), dtype=np.uint8)
    for si, intervals in enumerate(coverage):
        for s, e in intervals:
            known[si, s:e] = np.minimum(known[si, s:e] + 1, 2)
        for s, e in uncertain[si]:
            known[si, s:e] = 2
    primitive_parents = copy.deepcopy(parents)
    for b in primitive_parents:
        b["anchors"] = [
            unique_runs(b["alleles"][0], a, anchors)
            for a, anchors in zip(b["alleles"], b["anchors"])
        ]
    for b in parents:
        b["anchors"] = [
            refine_anchors(
                b["alleles"][0], a, unique_runs(b["alleles"][0], a, anchors), args.sequence_anchor
            )
            if b.get("native_kind", ["REFERENCE", "SEQUENCE"])[ai] == "SEQUENCE"
            else anchors
            for ai, (a, anchors) in enumerate(zip(b["alleles"], b["anchors"]))
        ]
    forward = decompose_alleles(parents)
    rc = reverse_complement
    reverse_parents = [
        dict(
            b,
            alleles=[rc(s) for s in b["alleles"]],
            anchors=[
                [
                    (len(b["alleles"][0]) - r - n, len(s) - q - n, n, node)
                    for r, q, n, node in anchors
                ][::-1]
                for s, anchors in zip(b["alleles"], b["anchors"])
            ],
        )
        for b in parents
    ]
    reverse_raw = decompose_alleles(reverse_parents)
    reverse_original = [
        [
            [
                Edit(
                    len(b["alleles"][0]) - e.end,
                    len(b["alleles"][0]) - e.start,
                    rc(e.alt),
                    e.complex,
                )
                for e in es
            ]
            for es in ee
        ]
        for b, ee in zip(parents, reverse_raw)
    ]
    edits, forward_failed = normalize(
        parents, forward, reference, args.reference, args.work / "normalized_forward"
    )
    reverse, reverse_failed = normalize(
        parents, reverse_original, reference, args.reference, args.work / "normalized_reverse"
    )
    stable_f, stable_r, compound = preserve_complex_gaps(
        parents,
        forward,
        reverse_original,
        edits,
        reverse,
        forward_failed | reverse_failed,
        reference,
    )
    if compound:
        edits, forward_failed = normalize(
            parents, stable_f, reference, args.reference, args.work / "compound_forward"
        )
        reverse, reverse_failed = normalize(
            parents, stable_r, reference, args.reference, args.work / "compound_reverse"
        )
        unresolved = forward_failed | reverse_failed
        for pi, b in enumerate(parents):
            for ai in range(1, len(b["alleles"])):
                if directional_disagreements(
                    reference, b["pos"] - 1, edits[pi][ai], reverse[pi][ai]
                ):
                    unresolved.add((pi, ai))
        if unresolved:
            # If two gaps interact after normalization, retain their complete
            # parent allele with its actual two-sided source context.
            subset = []
            raw = []
            keys = []
            for pi, ai in sorted(unresolved):
                b = parents[pi]
                subset.append(dict(b, alleles=[b["alleles"][0], b["alleles"][ai]]))
                raw.append([[], [Edit(1, len(b["alleles"][0]) - 1, b["alleles"][ai][1:-1], True)]])
                keys.append((pi, ai))
            full, fail = normalize(
                subset, raw, reference, args.reference, args.work / "compound_parents"
            )
            assert not fail, "A complete allele must normalize without losing its sequence"
            for (pi, ai), ee in zip(keys, full):
                edits[pi][ai] = ee[1]
                reverse[pi][ai] = ee[1]
        forward_failed = reverse_failed = set()
    print(
        f"Preserved complete complex gaps: {sum(c for pi, ai, c in compound)} in {len(compound)} alleles",
        flush=True,
    )
    for b, fw, rv in zip(parents, edits, reverse):
        for ai in range(1, len(fw)):
            fw[ai] = close_touching_units(b["alleles"][0], fw[ai])
            rv[ai] = close_touching_units(b["alleles"][0], rv[ai])
    from .primitive_projection import partial_projections

    partial_samples, partial = partial_projections(
        primitive_parents,
        edits,
        known,
        reference,
        args.reference,
        args.work / "partial_projection",
        args.jobs,
        args.aligner,
    )
    with (args.work / "partial_projection.pkl").open("wb") as f:
        pickle.dump((partial_samples, partial), f, protocol=5)
    print(
        f"Partially resolved source projection: {len(partial)} source alleles; {sum(map(len, partial_samples.values()))} sample/source pairs",
        flush=True,
    )
    ambiguous = np.zeros((n, length), dtype=bool)
    for key, (_, positions) in partial.items():
        for si in partial_samples[key]:
            if positions:
                ambiguous[si, list(positions)] = True
    ambiguous_rows = []
    for pi, (b, fw, rv) in enumerate(zip(parents, edits, reverse)):
        offset = b["pos"] - 1
        for ai, (x, y) in enumerate(zip(fw, rv)):
            positions = directional_disagreements(reference, offset, x, y)
            if (pi, ai) in forward_failed | reverse_failed:
                positions.update(range(offset, offset + len(b["alleles"][0])))
            for si in b["carriers"][ai]:
                if positions:
                    ambiguous[si, list(positions)] = True
            if positions:
                ambiguous_rows.append((b["gid"], ai, len(positions), len(b["carriers"][ai])))
    with gzip.open(args.output / "Graph.alignment_ambiguity.tsv.gz", "wt") as f:
        out = csv.writer(f, delimiter="\t", lineterminator="\n")
        out.writerow(["block_id", "allele_index", "ambiguous_positions", "observed_samples"])
        out.writerows(ambiguous_rows)
    print(
        f"Canonical direction ambiguity: {len(ambiguous_rows)} alleles; {int(ambiguous.sum())} sample/reference positions",
        flush=True,
    )
    with gzip.open(args.output / "Graph.anchors.tsv.gz", "wt") as f:
        out = csv.writer(f, delimiter="\t", lineterminator="\n")
        out.writerow(
            ["block_id", "allele_index", "ref_offset0", "alt_offset0", "length", "oriented_node"]
        )
        for b in parents:
            for ai, anchors in enumerate(b["anchors"]):
                for anchor in anchors:
                    out.writerow([b["gid"], ai, *anchor])
    with (args.work / "aligned.pkl").open("wb") as f:
        pickle.dump(edits, f, protocol=5)
    models = [set() for _ in samples]
    origins = defaultdict(set)
    event_kinds = defaultdict(set)
    with gzip.open(args.output / "Graph.alleles.tsv.gz", "wt") as f:
        out = csv.writer(f, delimiter="\t", lineterminator="\n")
        out.writerow(["block_id", "allele_index", "ref_start0", "ref_end0", "sequence"])
        for pi, (b, ee) in enumerate(zip(parents, edits)):
            offset = b["pos"] - 1
            for ai, (seq, es) in enumerate(zip(b["alleles"], ee)):
                out.writerow([b["gid"], ai, offset, offset + len(b["alleles"][0]), seq])
                absolute = {Edit(e.start + offset, e.end + offset, e.alt, e.complex) for e in es}
                assert all((0 <= e.start <= e.end <= length for e in absolute)), b["gid"]
                fallback = {
                    Edit(e.start + offset, e.end + offset, e.alt, e.complex)
                    for e in partial.get((pi, ai), ([], set()))[0]
                }
                for si in b["carriers"][ai]:
                    models[si].update(
                        fallback if si in partial_samples.get((pi, ai), set()) else absolute
                    )
                for e in absolute | fallback:
                    if e.complex:
                        native = b.get("native_kind", [])[ai]
                        event_kinds[e.start, e.end, e.alt].add(native)
                    origins[e.start, e.end, "EDIT" if e.structural else "SNP"].add(b["gid"])
    complete_alleles = {
        Edit(e.start + b["pos"] - 1, e.end + b["pos"] - 1, e.alt, e.complex)
        for b, ee in zip(parents, edits)
        for es in ee
        for e in es
    }
    promoted = cohere_complex_models(models, reference, known, complete_alleles)
    print(
        f"Unified equivalent complete allele representations: {len(promoted)} sample/locus pairs",
        flush=True,
    )
    with (args.work / "representation_changes.json").open("w") as f:
        json.dump(promoted, f)
    callable_mask = known == 1
    bases = np.broadcast_to(np.frombuffer(reference.encode(), dtype=np.uint8), (n, length)).copy()
    bases[~callable_mask] = ord("?")
    all_bins = []
    for si, mm in enumerate(models):
        all_bb = defaultdict(set)
        for e in sorted(mm):
            for k in range(e.start // 1000, e.end // 1000 + 1):
                all_bb[k].add(e)
            if e.structural:
                bases[si, e.start : e.end] = ord("-")
            else:
                bases[si, e.start : e.end] = np.frombuffer(e.alt.encode(), dtype=np.uint8)
        all_bins.append(all_bb)
        # Individually equivalent normalizations can interact with a nearby
        # allele. Such coupled edits must not yield confident primitive calls.
        previous = None
        for e in sorted(mm):
            if previous and (
                e.start < previous.end or e.start == e.end == previous.start == previous.end
            ):
                ambiguous[
                    si,
                    max(0, min(e.start, previous.start) - 1) : min(
                        length, max(e.end, previous.end) + 1
                    ),
                ] = True
            if previous is None or e.end > previous.end:
                previous = e
    from .validate_models import validate_models

    model_validation = validate_models(parents, models, reference, known, ambiguous)
    (args.work / "model_validation.json").write_text(json.dumps(model_validation, indent=2) + "\n")
    print(json.dumps(model_validation), flush=True)
    bases[~callable_mask] = ord("?")
    ambiguity_prefix = np.pad(np.cumsum(ambiguous, axis=1, dtype=np.int32), ((0, 0), (1, 0)))
    multiple_prefix = np.pad(np.cumsum(known >= 2, axis=1, dtype=np.int32), ((0, 0), (1, 0)))
    prefix = np.pad(np.cumsum(~callable_mask, axis=1, dtype=np.int32), ((0, 0), (1, 0)))
    chain_index = []
    for intervals in chains:
        intervals = sorted(intervals)
        starts = [s for s, e in intervals]
        ends = []
        end = -1
        for s, e in intervals:
            end = max(end, e)
            ends.append(end)
        chain_index.append((starts, ends))

    def observe(si, s, e, kind):
        if kind == "SNP":
            if known[si, s] >= 2:
                return None, "MULTIPLE_OCCURRENCES", []
            if ambiguous[si, s]:
                return None, "AMBIGUOUS_ALIGNMENT", []
            base = chr(int(bases[si, s]))
            if base in DNA:
                return base, "PASS", []
            return None, {"?": "NO_DATA", "-": "SPANNING_EDIT"}.get(base, "AMBIGUOUS_BASE"), []
        lo, hi = max(0, s - 1), min(length, e + 1)
        interior_lo, interior_hi = (s, e) if s < e else (lo, hi)
        if multiple_prefix[si, interior_hi] != multiple_prefix[si, interior_lo]:
            return None, "MULTIPLE_OCCURRENCES", []
        if ambiguity_prefix[si, interior_hi] != ambiguity_prefix[si, interior_lo]:
            return None, "AMBIGUOUS_ALIGNMENT", []
        cs, ce = chain_index[si]
        ci = bisect_right(cs, lo) - 1
        if ci < 0 or ce[ci] < hi or prefix[si, interior_hi] != prefix[si, interior_lo]:
            return None, "NO_DATA", []
        local = set().union(
            *(all_bins[si].get(k, set()) for k in range(lo // 1000, hi // 1000 + 1))
        )
        local_all = local
        focal = Edit(s, e, "")
        local = sorted(x for x in local if intersects(focal, x) and (s != e or x.structural))
        if any(x.start < s or x.end > e for x in local):
            return None, "OVERLAPPING_EDIT", local
        flanks = [v for v in (s - 1, e) if 0 <= v < length]
        if any(set(x.alt) - DNA for x in local_all for v in flanks if x.start <= v < x.end):
            return None, "UNCERTAIN_BREAKPOINT", local
        try:
            observed = replay(
                reference[s:e], [Edit(x.start - s, x.end - s, x.alt, x.complex) for x in local]
            )
        except ValueError:
            return None, "OVERLAPPING_EDIT", local
        if set(observed) - DNA:
            return None, "PARTIAL_SEQUENCE", local
        return observed, "PASS", local

    header = [
        "##fileformat=VCFv4.2",
        f"##source=PathAllele-{__version__}",
        f"##contig=<ID={args.reference.read_text().splitlines()[0][1:].split()[0]},length={length}>",
        '##INFO=<ID=PID,Number=.,Type=String,Description="Supporting native path source contexts, not nested parent loci">',
        '##INFO=<ID=SVLEN,Number=A,Type=Integer,Description="ALT minus REF length">',
        '##INFO=<ID=NATIVE,Number=A,Type=String,Description="Native orientation class; SEQUENCE, INV or COMPLEX">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Haploid genotype">',
        '##FORMAT=<ID=FT,Number=1,Type=String,Description="PASS; PARTIAL_SEQUENCE; AMBIGUOUS_ALIGNMENT indicates alternative concrete base mappings; other statuses are unresolved">',
    ]
    from .variant_schema import VCF_HEADERS, classification_info

    header.extend(VCF_HEADERS)
    counts = Counter()
    positive = [set() for _ in samples]
    catalogue = JointCatalogue(samples, args.output)
    with (args.work / "Graph.raw.vcf").open("w") as f:
        f.write(
            "\n".join(header)
            + "\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
            + "\t".join(samples)
            + "\n"
        )
        seen = set()
        # Independent components precede wider complete alleles. This ordering
        # permits REF-backbone projection only after every residual edit has a
        # concrete emitted ALT witness; it never imputes a missing genotype.
        ordered_loci = sorted(origins, key=lambda z: (z[2] != "SNP", z[1] - z[0], z))
        locus_ids = {key: f"v{i}" for i, key in enumerate(sorted(origins))}
        for vi, (s, e, kind) in enumerate(ordered_loci):
            ref = reference[s:e]
            observations, states = [], []
            for si in range(n):
                observed, ft, local = observe(si, s, e, kind)
                observations.append((observed, ft))
                states.append(
                    (observed, ft)
                    if observed is None or kind == "SNP"
                    else factor_state(reference, s, e, observed, local, positive[si])
                )
            catalogue.add(locus_ids[s, e, kind], s, e, ref, observations, states, kind)
            alts = sorted(
                {a for a, ft in states if a is not None and a != ref}, key=lambda a: (len(a), a)
            )
            if not alts:
                continue
            lookup = {a: i for i, a in enumerate([ref] + alts)}
            pos, vr, va = vcf_alleles(reference, 0, s, e, alts, kind)
            for alt in va:
                assert (pos, vr, alt) not in seen
                seen.add((pos, vr, alt))
            info = dict(
                NATIVE=",".join(
                    "+".join(sorted(event_kinds.get((s, e, aa), {"SEQUENCE"}))) for aa in alts
                ),
                PID=",".join(sorted(origins[s, e, kind])),
                SVLEN=",".join((str(len(a) - len(vr)) for a in va)),
            )
            for tag, values in classification_info(vr, va).items():
                info[tag] = ",".join(map(str, values)) if isinstance(values, tuple) else values
            cells = []
            for si, (a, ft) in enumerate(states):
                cells.append(("." if a is None else str(lookup[a])) + ":" + ft)
                counts["GT_" + ft] += 1
                if a is not None and a != ref:
                    positive[si].update(
                        (x for x in (Edit(s, e, a), Edit(s, e, a, True)) if x in models[si])
                    )
            f.write(
                "\t".join(
                    [
                        args.reference.read_text().splitlines()[0][1:].split()[0],
                        str(pos),
                        locus_ids[s, e, kind],
                        vr,
                        ",".join(va),
                        ".",
                        ".",
                        ";".join((f"{k}={v}" for k, v in info.items())),
                        "GT:FT",
                    ]
                    + cells
                )
                + "\n"
            )
            counts["records"] += 1
            counts["ALT_alleles"] += len(alts)
            if (vi + 1) % 5000 == 0:
                print(f"Projected {vi + 1}/{len(origins)} sites", flush=True)
    catalogue.close()
    with (args.work / "projection.pkl").open("wb") as f:
        pickle.dump((models, positive, callable_mask), f, protocol=5)
    summary = dict(method="native_path_projection", samples=n, counts=dict(counts))
    (args.output / "summary.projection.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
