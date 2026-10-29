"""Stable local projection for source alleles crossing unresolved territory.

Complete allele boundaries and their edit types are established independently.
This projection uses no known-site list or external callset.
"""

import os
import subprocess
from concurrent.futures import ThreadPoolExecutor
from .allele_model import DNA, Edit, atomic_edits, replay
from .graph_anchors import anchored_gaps
from .graph_io import reverse_complement
from .normalize_alleles import normalize, directional_disagreements


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


def partial_projections(
    parents, complete, known, reference, reference_path, work, jobs, executable
):
    requested = {}
    for pi, (parent, alleles) in enumerate(zip(parents, complete)):
        off = parent["pos"] - 1
        for ai, edits in enumerate(alleles):
            selected = []
            for si in parent["carriers"][ai]:
                if any(
                    not (known[si, off + e.start : off + e.end] == 1).all()
                    for e in edits
                    if e.end > e.start
                ):
                    selected.append(si)
            if selected:
                requested[pi, ai] = set(selected)
    if not requested:
        return {}, {}
    keys = sorted(requested)
    subset = []
    for pi, ai in keys:
        p = parents[pi]
        subset.append(
            dict(
                p,
                alleles=[p["alleles"][0], p["alleles"][ai]],
                anchors=[p["anchors"][0], p["anchors"][ai]],
                native_kind=["REFERENCE", p["native_kind"][ai]],
            )
        )
    work.mkdir()
    fw = work / "forward"
    fw.mkdir()
    rv = work / "reverse"
    rv.mkdir()
    forward = align_alleles(subset, fw, jobs, executable)
    rc = reverse_complement
    reversed_parents = [
        dict(
            p,
            alleles=[rc(s) for s in p["alleles"]],
            anchors=[
                [
                    (len(p["alleles"][0]) - r - n, len(s) - q - n, n, node)
                    for r, q, n, node in anchors
                ][::-1]
                for s, anchors in zip(p["alleles"], p["anchors"])
            ],
        )
        for p in subset
    ]
    reverse_raw = align_alleles(reversed_parents, rv, jobs, executable)
    reverse = [
        [
            [
                Edit(
                    len(p["alleles"][0]) - e.end,
                    len(p["alleles"][0]) - e.start,
                    rc(e.alt),
                    e.complex,
                )
                for e in es
            ]
            for es in alleles
        ]
        for p, alleles in zip(subset, reverse_raw)
    ]
    forward, _ = normalize(subset, forward, reference, reference_path, work / "normalized_forward")
    reverse, _ = normalize(subset, reverse, reference, reference_path, work / "normalized_reverse")
    result = {}
    for key, p, fw, rv in zip(keys, subset, forward, reverse):
        positions = directional_disagreements(reference, p["pos"] - 1, fw[1], rv[1])
        # A rejected normalization already returns the sequence-valid source
        # edits. Only actual forward/reverse disagreements are uncertain; a
        # conflict in one repeat must not mask stable positions elsewhere.
        result[key] = (fw[1], positions)
    return requested, result
