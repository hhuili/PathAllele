"""Check joint source models against raw physical alleles before VCF export."""

from collections import Counter, defaultdict
from .allele_model import Edit, intersects, replay


def physical_models(parents, sample_count):
    models = [set() for _ in range(sample_count)]
    for parent in parents:
        start = parent["pos"]
        end = parent["pos"] + len(parent["alleles"][0]) - 2
        for ai, allele in enumerate(parent["alleles"][1:], 1):
            edit = Edit(start, end, allele[1:-1], True)
            for si in parent["carriers"][ai]:
                models[si].add(edit)
    return models


def indexed(edits):
    bins = defaultdict(set)
    for edit in edits:
        for i in range(edit.start // 1000, edit.end // 1000 + 1):
            bins[i].add(edit)
    return bins


def local(index, start, end):
    candidates = set().union(
        *(index.get(i, set()) for i in range(max(0, start - 1) // 1000, (end + 1) // 1000 + 1))
    )
    return {edit for edit in candidates if intersects(Edit(start, end, ""), edit)}


def validate_models(parents, models, reference, known, ambiguous):
    counts = Counter()
    physical = physical_models(parents, len(models))
    for si, (raw, model) in enumerate(zip(physical, models)):
        left, right = indexed(raw), indexed(model)
        done = set()
        for unit in raw:
            s, t = unit.start, unit.end
            while True:
                original = local(left, s, t)
                current = local(right, s, t)
                both = original | current
                lo = min([s] + [e.start for e in both])
                hi = max([t] + [e.end for e in both])
                if (lo, hi) == (s, t):
                    break
                s, t = lo, hi
            if (s, t) in done:
                continue
            done.add((s, t))
            lo, hi = (s, t) if s < t else (max(0, s - 1), min(len(reference), t + 1))
            if not (known[si, lo:hi] == 1).all() or ambiguous[si, lo:hi].any():
                counts["source_windows_with_unresolved_territory"] += 1
                continue

            def sequence(edits):
                return replay(
                    reference[s:t], [Edit(e.start - s, e.end - s, e.alt, e.complex) for e in edits]
                )

            expected, observed = sequence(original), sequence(current)
            if expected != observed:
                raise ValueError(
                    f"Joint model changed physical sequence: sample index {si}, interval {s}:{t}, lengths {len(expected)}:{len(observed)}"
                )
            counts["joint_source_windows_replayed"] += 1
    return dict(counts)
