#!/usr/bin/env python3
"""Align native allele intervals and project one shared haploid variant catalogue."""

import argparse
from collections import Counter, defaultdict
from bisect import bisect_right
from concurrent.futures import ThreadPoolExecutor
import csv
import gzip
import json
import os
from pathlib import Path
import pickle
import subprocess
import numpy as np
from . import __version__
from .graph_io import reverse_complement
from .allele_model import DNA, Edit, atomic_edits, intersects, replay, vcf_alleles
from .graph_anchors import anchored_gaps, unique_runs
from .normalize_alleles import (
    normalize,
    directional_disagreements,
    preserve_complex_gaps,
    cohere_complex_models,
)


def align_alleles(parents, work, jobs, executable):
    gaps = {
        (pi, ai): anchored_gaps(p["alleles"][0], a, p["anchors"][ai])
        for pi, p in enumerate(parents)
        for ai, a in enumerate(p["alleles"][1:], 1)
        if p.get("native_kind", ["REFERENCE"] + ["SEQUENCE"] * (len(p["alleles"]) - 1))[ai]
        == "SEQUENCE"
    }
    tasks = [
        (len(ref) + len(alt), pi, ai, segment)
        for (pi, ai), parts in gaps.items()
        for segment, (start, ref, alt) in enumerate(parts)
        if ref and alt and not (set(ref + alt) - DNA)
    ]
    buckets, loads = ([[] for _ in range(jobs)], [0] * jobs)
    for weight, pi, ai, segment in sorted(tasks, reverse=True):
        j = min(range(jobs), key=loads.__getitem__)
        buckets[j].append((pi, ai, segment))
        loads[j] += weight
    header = '##fileformat=VCFv4.2\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Synthetic allele carrier; ignored in inference">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tallele\n'
    for j, bucket in enumerate(buckets):
        with (work / f"{j}.alleles.vcf").open("w") as out:
            out.write(header)
            for pi, ai, segment in bucket:
                p = parents[pi]
                start, ref, alt = gaps[pi, ai][segment]
                out.write(
                    "\t".join(
                        [
                            p["chrom"],
                            str(p["pos"] + start),
                            f"a{pi}_{ai}_{segment}",
                            ref,
                            alt,
                            ".",
                            ".",
                            ".",
                            "GT",
                            "1",
                        ]
                    )
                    + "\n"
                )

    def run(j):
        with (work / f"{j}.edits.tsv").open("w") as out, (work / f"{j}.log").open("w") as log:
            subprocess.run(
                [executable, str(work / f"{j}.alleles.vcf")],
                stdout=out,
                stderr=log,
                check=True,
                env=dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1"),
            )
        print(f"Allele alignments: batch {j + 1}/{jobs} complete", flush=True)

    with ThreadPoolExecutor(max_workers=jobs) as pool:
        list(pool.map(run, range(jobs)))
    edits = [[[] for _ in p["alleles"]] for p in parents]
    for (pi, ai), parts in gaps.items():
        for start, ref, alt in parts:
            if set(ref + alt) - DNA:
                # An assembly gap is not a measured deletion length. Preserve
                # the bounded sequence, but do not align through unknown bases.
                edits[pi][ai].append(Edit(start, start + len(ref), alt, True))
            elif not ref or not alt:
                edits[pi][ai].extend(atomic_edits(start, ref, alt))
    for j in range(jobs):
        with (work / f"{j}.edits.tsv").open() as handle:
            for line in handle:
                if line.startswith("#"):
                    continue
                x = line.rstrip().split("\t")
                pi, ai, segment = (int(v) for v in x[0][1:].split("_"))
                p = parents[pi]
                start = int(x[1]) - p["pos"]
                ref, alt = ("" if v == "-" else v for v in x[2:4])
                gap_start, gap_ref, gap_alt = gaps[pi, ai][segment]
                assert gap_start <= start <= start + len(ref) <= gap_start + len(gap_ref)
                if p["alleles"][0][start : start + len(ref)] != ref:
                    raise ValueError("Alignment REF mismatch: " + p["pid"])
                edits[pi][ai].extend(atomic_edits(start, ref, alt, x[4] == "1"))
    for p, ee in zip(parents, edits):
        for ai, kind in enumerate(p.get("native_kind", [])):
            if kind in ("INV", "COMPLEX"):
                ee[ai] = [Edit(1, len(p["alleles"][0]) - 1, p["alleles"][ai][1:-1], True)]
        for ai, (allele, edit_list) in enumerate(zip(p["alleles"], ee)):
            if replay(p["alleles"][0], edit_list) != allele:
                raise ValueError(f"Sequence decomposition failed: {p['pid']} allele {ai}")
    return edits


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--blocks", type=Path, required=True)
    p.add_argument("--coverage", type=Path, required=True)
    p.add_argument("--reference", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--work", type=Path, required=True)
    p.add_argument("--jobs", type=int, default=8)
    p.add_argument("--aligner", required=True)
    args = p.parse_args(argv)
    args.work.mkdir(parents=True, exist_ok=False)
    samples, parents = pickle.load(args.blocks.open("rb"))
    coverage, uncertain, chains = pickle.load(args.coverage.open("rb"))
    reference = "".join(
        (l for l in args.reference.read_text().splitlines() if not l.startswith(">"))
    )
    n = len(samples)
    length = len(reference)
    for b in parents:
        b["anchors"] = [
            unique_runs(b["alleles"][0], a, anchors)
            for a, anchors in zip(b["alleles"], b["anchors"])
        ]
    anchored_work = args.work / "anchored"
    anchored_work.mkdir()
    forward = align_alleles(parents, anchored_work, args.jobs, args.aligner)
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
    reverse_work = args.work / "reverse"
    reverse_work.mkdir()
    reverse_raw = align_alleles(reverse_parents, reverse_work, args.jobs, args.aligner)
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
    ambiguous = np.zeros((n, length), dtype=bool)
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
    compound_loci = set()
    known_alleles = defaultdict(set)
    with gzip.open(args.output / "Graph.alleles.tsv.gz", "wt") as f:
        out = csv.writer(f, delimiter="\t", lineterminator="\n")
        out.writerow(["block_id", "allele_index", "ref_start0", "ref_end0", "sequence"])
        for b, ee in zip(parents, edits):
            offset = b["pos"] - 1
            for ai, (seq, es) in enumerate(zip(b["alleles"], ee)):
                out.writerow([b["gid"], ai, offset, offset + len(b["alleles"][0]), seq])
                absolute = {Edit(e.start + offset, e.end + offset, e.alt, e.complex) for e in es}
                assert all((0 <= e.start <= e.end <= length for e in absolute)), b["gid"]
                for si in b["carriers"][ai]:
                    models[si].update(absolute)
                for e in absolute:
                    if e.structural:
                        known_alleles[e.start, e.end].add(e.alt)
                    if e.complex:
                        compound_loci.add((e.start, e.end))
                        native = b.get("native_kind", [])[ai]
                        event_kinds[e.start, e.end, e.alt].add(
                            "COMPLEX" if native == "SEQUENCE" else native
                        )
                    origins[e.start, e.end, "EDIT" if e.structural else "SNP"].add(b["gid"])
    promoted = cohere_complex_models(models, reference)
    print(
        f"Unified equivalent complete allele representations: {len(promoted)} sample/locus pairs",
        flush=True,
    )
    with (args.work / "representation_changes.json").open("w") as f:
        json.dump(promoted, f)
    known = np.zeros((n, length), dtype=np.uint8)
    for si, intervals in enumerate(coverage):
        for s, e in intervals:
            known[si, s:e] = np.minimum(known[si, s:e] + 1, 2)
        for s, e in uncertain[si]:
            known[si, s:e] = 2
    callable_mask = known == 1
    bases = np.broadcast_to(np.frombuffer(reference.encode(), dtype=np.uint8), (n, length)).copy()
    bases[~callable_mask] = ord("?")
    structural = []
    bins = []
    all_bins = []
    for si, mm in enumerate(models):
        ss = sorted((e for e in mm if e.structural))
        structural.append(ss)
        bb = defaultdict(set)
        all_bb = defaultdict(set)
        for e in sorted(mm):
            for k in range(e.start // 1000, e.end // 1000 + 1):
                all_bb[k].add(e)
            if e.structural:
                bases[si, e.start : e.end] = ord("-")
                for k in range(e.start // 1000, e.end // 1000 + 1):
                    bb[k].add(e)
            else:
                bases[si, e.start : e.end] = np.frombuffer(e.alt.encode(), dtype=np.uint8)
        bins.append(bb)
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

    def state(si, s, e, kind):
        if kind == "SNP":
            if known[si, s] >= 2:
                return (None, "MULTIPLE_OCCURRENCES")
            if ambiguous[si, s]:
                return (None, "AMBIGUOUS_ALIGNMENT")
            base = chr(int(bases[si, s]))
            return (
                (base, "PASS")
                if base in DNA
                else (None, {"?": "NO_DATA", "-": "SPANNING_EDIT"}.get(base, "AMBIGUOUS_BASE"))
            )
        if multiple_prefix[si, min(length, e + 1)] != multiple_prefix[si, max(0, s - 1)]:
            return (None, "MULTIPLE_OCCURRENCES")
        if ambiguity_prefix[si, min(length, e + 1)] != ambiguity_prefix[si, max(0, s - 1)]:
            return (None, "AMBIGUOUS_ALIGNMENT")
        cs, ce = chain_index[si]
        ci = bisect_right(cs, max(0, s - 1)) - 1
        if ci < 0 or ce[ci] < min(length, e + 1):
            return (None, "NO_DATA")
        if prefix[si, e] != prefix[si, s]:
            return (None, "NO_DATA")
        if (s, e) in compound_loci:
            local = set().union(
                *(
                    all_bins[si].get(k, set())
                    for k in range(max(0, s - 1) // 1000, (e + 1) // 1000 + 1)
                )
            )
            local = [x for x in local if intersects(Edit(s, e, ""), x)]
            if any(x.start < s or x.end > e for x in local):
                return (None, "OVERLAPPING_EDIT")
            if any(chr(int(bases[si, v])) not in DNA for v in (s - 1, e) if 0 <= v < length):
                return (None, "UNCERTAIN_BREAKPOINT")
            try:
                observed = replay(
                    reference[s:e], [Edit(x.start - s, x.end - s, x.alt, x.complex) for x in local]
                )
            except ValueError:
                return (None, "OVERLAPPING_EDIT")
            if set(observed) - DNA:
                return (None, "PARTIAL_SEQUENCE")
            # Match complete observed sequence, regardless of which primitive
            # decomposition represented it in this particular sample.
            if observed == reference[s:e]:
                return (observed, "PASS")
            if observed in known_alleles[s, e]:
                if not any((x.start, x.end, x.alt) == (s, e, observed) for x in local):
                    return (None, "AMBIGUOUS_REPRESENTATION")
                return (observed, "PASS")
            return (None, "UNREPRESENTED_ALLELE")
        local = set().union(
            *(bins[si].get(k, set()) for k in range(max(0, s - 1) // 1000, (e + 1) // 1000 + 1))
        )
        focal = Edit(s, e, "")
        overlap = [x for x in local if intersects(focal, x)]
        exact = [x for x in overlap if (x.start, x.end) == (s, e)]
        if len(exact) > 1:
            return (None, "MULTIPLE_OCCURRENCES")
        if exact:
            if any((x != exact[0] for x in overlap)):
                return (None, "OVERLAPPING_EDIT")
            junctions = [v for v in [s - 1, e] if 0 <= v < length]
            if any((chr(int(bases[si, v])) not in DNA for v in junctions)):
                return (None, "UNCERTAIN_BREAKPOINT")
            alt = exact[0].alt
            if set(alt) - DNA:
                return (None, "PARTIAL_SEQUENCE")
            return (alt, "PASS")
        if overlap:
            return (None, "OVERLAPPING_EDIT")
        junctions = set((v for v in [s - 1, s, e - 1, e] if 0 <= v < length))
        if any((chr(int(bases[si, v])) not in DNA for v in junctions)):
            return (None, "UNCERTAIN_BREAKPOINT")
        return (reference[s:e], "PASS")

    header = [
        "##fileformat=VCFv4.2",
        f"##source=PathAllele-{__version__}",
        f"##contig=<ID={args.reference.read_text().splitlines()[0][1:].split()[0]},length={length}>",
        '##INFO=<ID=PID,Number=.,Type=String,Description="Supporting native path intervals">',
        '##INFO=<ID=SOURCE,Number=1,Type=String,Description="PATH_INTERVAL">',
        '##INFO=<ID=TYPE,Number=1,Type=String,Description="SNP or EDIT">',
        '##INFO=<ID=SVLEN,Number=A,Type=Integer,Description="ALT minus REF length">',
        '##INFO=<ID=NATIVE,Number=A,Type=String,Description="Native orientation class; SEQUENCE, INV or COMPLEX">',
        '##INFO=<ID=NBASES,Number=A,Type=Integer,Description="Unknown bases in ALT">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Haploid genotype">',
        '##FORMAT=<ID=FT,Number=1,Type=String,Description="PASS; PARTIAL_SEQUENCE; AMBIGUOUS_ALIGNMENT indicates alternative concrete base mappings; other statuses are unresolved">',
    ]
    counts = Counter()
    positive = [set() for _ in samples]
    with (args.output / "Graph.vcf").open("w") as f:
        f.write(
            "\n".join(header)
            + "\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
            + "\t".join(samples)
            + "\n"
        )
        seen = set()
        for vi, (s, e, kind) in enumerate(sorted(origins)):
            states = [state(si, s, e, kind) for si in range(n)]
            ref = reference[s:e]
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
                SOURCE="PATH_INTERVAL",
                TYPE=kind,
                SVLEN=",".join((str(len(a) - len(vr)) for a in va)),
                NBASES=",".join((str(sum((b not in DNA for b in a))) for a in va)),
            )
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
                        f"v{vi}",
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
    with (args.work / "projection.pkl").open("wb") as f:
        pickle.dump((models, positive, callable_mask), f, protocol=5)
    summary = dict(method="native_path_projection", samples=n, counts=dict(counts))
    (args.output / "summary.projection.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
