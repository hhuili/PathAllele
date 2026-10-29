"""Exact, unique internal homology inside an already resolved physical allele.

This module sees complete sequences, never a CIGAR, a variant catalogue or a
population label. It cannot extend a physical source interval or resolve a
competing occurrence. Tied collinear chains retain only obligatory anchors.
"""

from collections import defaultdict

from .graph_io import reverse_complement


def obligatory_chain(runs):
    """Anchors present in every maximum-base chain; both axes may overlap."""
    runs = sorted(set(runs))
    if not runs:
        return []

    def precedes(a, b):
        return a[0] + a[2] <= b[0] and a[1] + a[2] <= b[1]

    def scores(items):
        values = []
        for i, run in enumerate(items):
            candidates = [values[j] for j in range(i) if precedes(items[j], run)]
            best = max([0] + [v[0] for v in candidates])
            count = sum(v[1] for v in candidates if v[0] == best) if best else 1
            values.append((best + run[2], count))
        best = max(v[0] for v in values)
        return values, (best, sum(v[1] for v in values if v[0] == best))

    forward, total = scores(runs)
    reflected = sorted((-r - n, -q - n, n) for r, q, n in runs)
    reverse_scores, _ = scores(reflected)
    reverse_lookup = dict(zip(reflected, reverse_scores))
    backward = [reverse_lookup[-r - n, -q - n, n] for r, q, n in runs]
    return [
        run
        for run, f, b in zip(runs, forward, backward)
        if f[0] + b[0] - run[2] == total[0] and f[1] * b[1] == total[1]
    ]


def exact_runs(ref, alt, minimum=31):
    """Maximal exact matches seeded uniquely in both complete gap sequences."""
    if minimum < 1:
        raise ValueError("Sequence anchor length must be positive")
    if min(len(ref), len(alt)) < minimum or set(ref + alt) - set("ACGT"):
        return []
    indices = []
    for sequence in (ref, alt):
        index = defaultdict(list)
        for i in range(len(sequence) - minimum + 1):
            index[sequence[i : i + minimum]].append(i)
        indices.append(index)
    seeds = sorted(
        (qp[0], indices[0][seq][0])
        for seq, qp in indices[1].items()
        if len(qp) == 1 and len(indices[0].get(seq, [])) == 1
    )
    # Coalesce overlapping seeds on one diagonal before extension. Long exact
    # runs are extended once, rather than once per k-mer.
    merged = []
    for q, r in seeds:
        if merged and r - q == merged[-1][0] - merged[-1][1] and q <= merged[-1][1] + merged[-1][2]:
            merged[-1][2] = max(merged[-1][2], q + minimum - merged[-1][1])
        else:
            merged.append([r, q, minimum])
    matches = set()
    rc_ref, rc_alt = reverse_complement(ref), reverse_complement(alt)
    for r, q, n in merged:
        while r > 0 and q > 0 and ref[r - 1] == alt[q - 1]:
            r -= 1
            q -= 1
            n += 1
        while r + n < len(ref) and q + n < len(alt) and ref[r + n] == alt[q + n]:
            n += 1
        sequence = ref[r : r + n]
        if sequence in rc_ref or sequence in rc_alt:
            continue
        if ref.find(sequence) != r or ref.find(sequence, r + 1) >= 0:
            continue
        if alt.find(sequence) != q or alt.find(sequence, q + 1) >= 0:
            continue
        matches.add((r, q, n))
    return obligatory_chain(matches)


def refine_anchors(ref, alt, anchors, minimum=31):
    """Add internal sequence anchors without crossing physical graph anchors."""
    ordered = [(0, 0, 1, None)] + list(anchors) + [(len(ref) - 1, len(alt) - 1, 1, None)]
    result = list(anchors)
    re = qe = 0
    for r, q, n, node in ordered:
        if r < re or q < qe or ref[r : r + n] != alt[q : q + n]:
            raise ValueError("Invalid physical anchors")
        if ref[re:r] != alt[qe:q]:
            for rr, qq, length in exact_runs(ref[re:r], alt[qe:q], minimum):
                result.append((re + rr, qe + qq, length, None))
        re, qe = r + n, q + n
    return sorted(result, key=lambda x: (x[0], x[1]))
