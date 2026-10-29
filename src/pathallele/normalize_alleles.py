"""Use bcftools for canonical edit representation, retaining source replay."""

import subprocess
from pathlib import Path
from .allele_model import Edit, atomic_edits, replay, vcf_alleles


def bounded_indel_normalization(reference, edits):
    """Left-align simple indels inside this source's physical flanking bases."""
    result = []
    for e in edits:
        s, t, alt = e.start, e.end, e.alt
        if s == t and alt:
            while s > 1 and reference[s - 1] == alt[-1]:
                alt = alt[-1] + alt[:-1]
                s -= 1
                t -= 1
        elif not alt and s < t:
            while s > 1 and reference[s - 1] == reference[t - 1]:
                s -= 1
                t -= 1
        result.append(Edit(s, t, alt, e.complex))
    return result


def normalize(parents, edits, reference, reference_path, work):
    work.mkdir()
    source = work / "input.vcf"
    output = work / "normalized.vcf"
    chrom = Path(reference_path).read_text().splitlines()[0][1:].split()[0]
    with source.open("w") as f:
        f.write(
            f"##fileformat=VCFv4.2\n##contig=<ID={chrom},length={len(reference)}>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        )
        for pi, (b, ee) in enumerate(zip(parents, edits)):
            for ai, es in enumerate(ee):
                for ei, e in enumerate(es):
                    pos, ref, alts = vcf_alleles(
                        reference,
                        b["pos"] - 1,
                        e.start,
                        e.end,
                        [e.alt],
                        "EDIT" if e.structural else "SNP",
                    )
                    f.write(f"{chrom}\t{pos}\te{pi}_{ai}_{ei}\t{ref}\t{alts[0]}\t.\t.\t.\n")
    with (work / "normalization.log").open("w") as log:
        subprocess.run(
            [
                "bcftools",
                "norm",
                "-c",
                "e",
                "-f",
                str(reference_path),
                str(source),
                "-Ov",
                "-o",
                str(output),
            ],
            stderr=log,
            check=True,
        )
    result = [[[] for a in b["alleles"]] for b in parents]
    with output.open() as f:
        for line in f:
            if line.startswith("#"):
                continue
            x = line.rstrip().split("\t")
            pi, ai, ei = map(int, x[2][1:].split("_"))
            b = parents[pi]
            result[pi][ai].extend(
                atomic_edits(int(x[1]) - b["pos"], x[3], x[4], edits[pi][ai][ei].complex)
            )
    failed = set()
    changed = 0
    for pi, (b, ee) in enumerate(zip(parents, result)):
        offset = b["pos"] - 1
        for ai, (allele, es) in enumerate(zip(b["alleles"], ee)):
            lo = min([offset] + [offset + e.start for e in es])
            hi = max([offset + len(b["alleles"][0])] + [offset + e.end for e in es])
            expected = reference[lo:offset] + allele + reference[offset + len(b["alleles"][0]) : hi]
            try:
                valid = (
                    replay(
                        reference[lo:hi],
                        [
                            Edit(offset + e.start - lo, offset + e.end - lo, e.alt, e.complex)
                            for e in es
                        ],
                    )
                    == expected
                )
            except ValueError:
                valid = False
            # Independent source alleles must not normalize across their
            # physical boundary and collapse a second repeat-copy change.
            bounded = all(1 <= e.start <= e.end <= len(b["alleles"][0]) - 1 for e in es)
            if valid and not bounded:
                candidate = bounded_indel_normalization(b["alleles"][0], edits[pi][ai])
                try:
                    valid = replay(b["alleles"][0], candidate) == allele
                except ValueError:
                    valid = False
                if valid:
                    result[pi][ai] = candidate
                    continue
            if not valid:
                failed.add((pi, ai))
                result[pi][ai] = edits[pi][ai]
            else:
                changed += set(es) != set(edits[pi][ai])
    print(
        f"bcftools normalization: {changed} allele decompositions changed; {len(failed)} coupled conflicts remain unresolved",
        flush=True,
    )
    return result, failed


def directional_disagreements(reference, offset, left, right):
    """Exact reference-coordinate disagreements after canonical normalization."""
    if set(left) == set(right):
        return set()
    all_edits = set(left) | set(right)
    lo = min(offset + e.start for e in all_edits)
    hi = max(offset + e.end for e in all_edits)
    lo = max(0, lo - 1)
    hi = min(len(reference), hi + 1)

    def bases(es):
        out = list(reference[lo:hi])
        for e in es:
            s, t = offset + e.start - lo, offset + e.end - lo
            out[s:t] = ["?"] * (t - s) if e.structural else list(e.alt)
        return out

    aa, bb = bases(left), bases(right)
    positions = {lo + i for i, (x, y) in enumerate(zip(aa, bb)) if x != y}
    for e in set(left) ^ set(right):
        if e.structural:
            positions.update(
                range(max(0, offset + e.start - 1), min(len(reference), offset + e.end + 1))
            )
    return positions


def preserve_complex_gaps(
    parents, forward, reverse, normal_forward, normal_reverse, failed, reference
):
    """Keep a complete bounded sequence where primitive homology is unstable.

    Exact flanks still define the observed allele. A gap with multiple valid
    primitive representations is one compound allele, not a set of certain SNP
    positions. Stable adjacent gaps remain separately decomposed.
    """
    from .graph_anchors import anchored_gaps

    out_f = [[list(es) for es in ee] for ee in forward]
    out_r = [[list(es) for es in ee] for ee in reverse]
    changed = []
    for pi, b in enumerate(parents):
        offset = b["pos"] - 1
        for ai in range(1, len(b["alleles"])):
            pos = directional_disagreements(
                reference, offset, normal_forward[pi][ai], normal_reverse[pi][ai]
            )
            if not pos and (pi, ai) not in failed:
                continue
            raw_pos = directional_disagreements(reference, offset, forward[pi][ai], reverse[pi][ai])
            components = []
            for start, ref, alt in anchored_gaps(
                b["alleles"][0], b["alleles"][ai], b["anchors"][ai]
            ):
                if (pi, ai) in failed or any(
                    offset + start - 1 <= p <= offset + start + len(ref) for p in raw_pos
                ):
                    components.append((start, ref, alt))
            assert components, (b["gid"], ai)
            for source in (out_f, out_r):
                keep = source[pi][ai]
                for start, ref, alt in components:
                    end = start + len(ref)
                    keep = [e for e in keep if not (start <= e.start and e.end <= end)]
                    # Pure indels and single-base changes already have simple
                    # semantics; only a mixed replacement becomes compound.
                    compound = bool(ref and alt and not (len(ref) == len(alt) == 1))
                    keep.extend(atomic_edits(start, ref, alt, compound))
                assert replay(b["alleles"][0], keep) == b["alleles"][ai], b["gid"]
                source[pi][ai] = keep
            changed.append((pi, ai, len(components)))
    return out_f, out_r, changed


def cohere_complex_models(models, reference, known=None, complete_alleles=()):
    """Use one complete representation for identical observed complex alleles.

    Promote an equivalent primitive haplotype only on exact sequence replay.
    Broader or interacting edits are left unresolved. Descending spans make
    contained compound alternatives deterministic without overlapping ALT calls.
    """
    from collections import defaultdict
    from .allele_model import intersects

    sequences = defaultdict(set)
    compound_loci = {(e.start, e.end) for model in models for e in model if e.complex}
    for model in models:
        for e in model:
            if (
                e.structural
                and (e.start, e.end) in compound_loci
                and not (set(e.alt) - set("ACGT"))
            ):
                sequences[e.start, e.end].add(e.alt)
    for e in complete_alleles:
        if e.complex and not (set(e.alt) - set("ACGT")):
            sequences[e.start, e.end].add(e.alt)
    loci = sorted(sequences, key=lambda z: (-(z[1] - z[0]), z))
    promoted = []
    for si, model in enumerate(models):
        bins = defaultdict(set)
        for e in model:
            for k in range(e.start // 1000, e.end // 1000 + 1):
                bins[k].add(e)
        for s, t in loci:
            if known is not None and not (known[si, s:t] == 1).all():
                continue
            focal = Edit(s, t, "")
            local = set().union(
                *(bins.get(k, set()) for k in range(max(0, s - 1) // 1000, (t + 1) // 1000 + 1))
            )
            local = {e for e in local if intersects(focal, e)}
            if not local or any(e.start < s or e.end > t for e in local):
                continue
            if len(local) == 1 and next(iter(local)).complex:
                continue
            try:
                observed = replay(
                    reference[s:t], [Edit(e.start - s, e.end - s, e.alt, e.complex) for e in local]
                )
            except ValueError:
                continue
            if observed not in sequences[s, t]:
                continue
            replacement = Edit(s, t, observed, True)
            for e in local:
                model.remove(e)
                for k in range(e.start // 1000, e.end // 1000 + 1):
                    bins[k].discard(e)
            model.add(replacement)
            for k in range(s // 1000, t // 1000 + 1):
                bins[k].add(replacement)
            promoted.append((si, s, t, len(local)))
    return promoted


def close_touching_units(reference, edits):
    """A normalized boundary with no retained homology cannot split two units."""
    groups = []
    for edit in sorted(set(edits)):
        if groups and edit.start <= max(e.end for e in groups[-1]):
            groups[-1].append(edit)
        else:
            groups.append([edit])
    result = []
    for group in groups:
        if len(group) == 1:
            result.extend(group)
            continue
        s = min(e.start for e in group)
        t = max(e.end for e in group)
        alt = replay(
            reference[s:t], [Edit(e.start - s, e.end - s, e.alt, e.complex) for e in group]
        )
        result.extend(atomic_edits(s, reference[s:t], alt, keep_complex=True))
    assert replay(reference, result) == replay(reference, edits)
    return result
