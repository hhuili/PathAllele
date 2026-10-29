"""Occurrence-scoped anchors; no global node-to-reference base projection."""

from bisect import bisect_left, bisect_right
from collections import Counter, defaultdict


def mandatory_chain(runs):
    """Anchors present in every maximum exact-base-weight collinear chain.

    Runs are (reference offset, query offset, length, payload). Source intervals
    must not overlap. Reference intervals may overlap, hence end-aware queries.
    Exact integer path counts preserve ties instead of choosing a first match.
    """
    if not runs:
        return []
    runs = sorted(runs, key=lambda x: (x[1], x[0]))
    assert all(a[1] + a[2] <= b[1] for a, b in zip(runs, runs[1:])), "Overlapping query anchors"
    coordinates = sorted({x[0] for x in runs} | {x[0] + x[2] for x in runs})

    def combine(a, b):
        if a[0] != b[0]:
            return a if a[0] > b[0] else b
        return (a[0], a[1] + b[1])

    def scores(items):
        tree = [(0, 0)] * (len(coordinates) + 1)
        out = []

        def query(position):
            value = (0, 1)
            k = bisect_right(coordinates, position)
            while k:
                value = combine(value, tree[k])
                k -= k & -k
            return value

        for s, e, w in items:
            score, count = query(s)
            value = (score + w, count)
            out.append(value)
            k = bisect_left(coordinates, e) + 1
            while k < len(tree):
                tree[k] = combine(tree[k], value)
                k += k & -k
        return out, query(coordinates[-1])

    forward, total = scores([(r, r + n, n) for r, q, n, p in runs])
    limit = coordinates[-1] + coordinates[0]
    reverse_coordinates = sorted(limit - x for x in coordinates)
    # Reflected endpoints preserve the same coordinate set.
    old = coordinates
    coordinates = reverse_coordinates
    backward, _ = scores([(limit - (r + n), limit - r, n) for r, q, n, p in runs[::-1]])
    backward = backward[::-1]
    coordinates = old
    return [
        x
        for x, f, b in zip(runs, forward, backward)
        if f[0] + b[0] - x[2] == total[0] and f[1] * b[1] == total[1]
    ]


def context_runs(reference, query, anchors, min_length=1, sequence_unique=True):
    """Merge exact adjacent occurrences, require unique sequence, then chain."""
    runs = []
    for r, q, n, node in sorted(anchors, key=lambda x: (x[1], x[0])):
        assert reference[r : r + n] == query[q : q + n]
        if runs and r == runs[-1][0] + runs[-1][2] and q == runs[-1][1] + runs[-1][2]:
            runs[-1][2] += n
            runs[-1][3].append((r, q, n, node))
        else:
            runs.append([r, q, n, [(r, q, n, node)]])
    unique = []
    for r, q, n, payload in runs:
        if n < min_length:
            continue
        seq = reference[r : r + n]
        if set(seq) - set("ACGT"):
            continue
        if not sequence_unique or (
            reference.find(seq) == r
            and reference.find(seq, r + 1) < 0
            and query.find(seq) == q
            and query.find(seq, q + 1) < 0
        ):
            unique.append((r, q, n, payload))
    return mandatory_chain(unique)


def bounded_anchors(
    reference, query, rw, qw, nodes, rbase, qbase, context, ref_offset, query_offset
):
    """Local exact matches restricted to one physical, already bounded allele."""
    ri = defaultdict(list)
    qi = defaultdict(list)
    r = 0
    q = 0
    for node in rw:
        ri[node].append(r)
        r += len(nodes[abs(node)])
    for node in qw:
        qi[node].append(q)
        q += len(nodes[abs(node)])
    anchors = []
    for node, positions in qi.items():
        if len(positions) == 1 and len(ri.get(node, [])) == 1:
            anchors.append((ri[node][0] + rbase, positions[0] + qbase, len(nodes[abs(node)]), node))
    supported = context_supported(
        [(r + ref_offset, q + query_offset, n, node) for r, q, n, node in anchors], context
    )
    anchors = [(r - ref_offset, q - query_offset, n, node) for r, q, n, node in supported]
    return [
        item
        for run in context_runs(reference, query, anchors, sequence_unique=False)
        for item in run[3]
    ]


def context_supported(anchors, index):
    """Keep graph seeds with an exact, uniquely placed mapping overlap.

    Alignment and graph-node boundaries need not coincide. A long graph node
    requires at least 31 exact aligned bases at the same offset; a short node
    requires its complete sequence. Global anchor runs still require 31 bases.
    """
    if index is None:
        return []
    runs, starts, ends = index
    out = []
    for r, q, n, node in anchors:
        need = min(31, n)
        i = bisect_right(starts, q + n - need) - 1
        locations = set()
        while i >= 0 and ends[i] >= q + need:
            qs, qe, rs = runs[i]
            if min(q + n, qe) - max(q, qs) >= need:
                locations.add(rs + q - qs)
            i -= 1
        if locations == {r}:
            out.append((r, q, n, node))
    return out


def occurrence_seeds(walk, query_offsets, reference_occurrences, nodes, index):
    """Resolve graph-node occurrences by mapped physical offsets, not node ID.

    Repeated nodes may form a sequence-unique long context. They must not be
    discarded before contexts are merged. Each seed still requires one exact
    mature-alignment offset and an actual matching reference path occurrence.
    """
    if index is None:
        return []
    runs, starts, ends = index
    out = []
    for qi, node in enumerate(walk):
        q = query_offsets[qi]
        n = len(nodes[abs(node)])
        need = min(31, n)
        i = bisect_right(starts, q + n - need) - 1
        locations = set()
        while i >= 0 and ends[i] >= q + need:
            qs, qe, rs = runs[i]
            if min(q + n, qe) - max(q, qs) >= need:
                locations.add(rs + q - qs)
            i -= 1
        if len(locations) != 1:
            continue
        r = next(iter(locations))
        ri = reference_occurrences.get((node, r))
        if ri is not None:
            out.append((r, q, n, (node, ri, qi)))
    return out


def extend_terminal_context(reference, query, graph_runs, index, min_length=31):
    """Bound otherwise unplaced ends by unique exact sequence context.

    Graph construction can put homologous sequence in different nodes. Use
    mature mapping only outside the existing graph-bounded interval; do not
    replace its internal homology. Overlapping mappings must agree on offset.
    These are projection anchors, never direct PAF variant calls.
    """
    if index is None:
        return graph_runs
    events = defaultdict(list)
    for q, e, r in index[0]:
        events[q].append((r - q, 1))
        events[e].append((r - q, -1))
    active = Counter()
    segments = []
    previous = None
    for q in sorted(events):
        if previous is not None and q > previous and len(active) == 1:
            offset = next(iter(active))
            r = previous + offset
            n = q - previous
            if (
                segments
                and segments[-1][0] + segments[-1][2] == r
                and segments[-1][1] + segments[-1][2] == previous
            ):
                old = segments[-1]
                segments[-1] = (old[0], old[1], old[2] + n, None)
            else:
                segments.append((r, previous, n, None))
        for offset, delta in events[q]:
            active[offset] += delta
            if not active[offset]:
                del active[offset]
        previous = q
    if not graph_runs:
        return context_runs(reference, query, segments, min_length)
    first, last = graph_runs[0], graph_runs[-1]
    before = []
    after = []
    for r, q, n, payload in segments:
        length = min(n, first[0] - r, first[1] - q)
        if length > 0:
            before.append((r, q, length, None))
        skip = max(0, last[0] + last[2] - r, last[1] + last[2] - q)
        if skip < n:
            after.append((r + skip, q + skip, n - skip, None))
    return (
        context_runs(reference, query, before, min_length)
        + graph_runs
        + context_runs(reference, query, after, min_length)
    )


def mapped_copy_intervals(paf, samples, length, min_mapq=20, path_samples=None):
    """Mask distinct physical mapping components that cover the same locus.

    Overlapping query alignments belong to one physical component, including
    split/inversion mappings of that same copy. Disjoint query intervals or
    different paths are separate copies. This only removes callability; it
    cannot introduce an ALT or a reference genotype.
    """
    import numpy as np

    by_path = defaultdict(list)
    by_sample = defaultdict(list)
    with open(paf) as f:
        for line in f:
            x = line.rstrip().split("\t")
            tags = {v[:2]: v[5:] for v in x[12:]}
            if int(x[11]) >= min_mapq and tags.get("tp") == "P":
                by_path[x[0]].append((int(x[2]), int(x[3]), int(x[7]), int(x[8])))
    component = 0
    for path, records in by_path.items():
        end = -1
        for qs, qe, rs, re in sorted(records):
            if qs >= end:
                component += 1
            end = max(end, qe)
            sample = path_samples[path] if path_samples is not None else path.split("#")[0]
            by_sample[sample].append((rs, re, component))
    result = []
    for sample in samples:
        seen = np.zeros(length, dtype=np.int32)
        conflict = np.zeros(length, dtype=bool)
        for s, e, identity in by_sample[sample]:
            old = seen[s:e]
            conflict[s:e] |= (old != 0) & (old != identity)
            seen[s:e] = np.where(old == 0, identity, old)
        edges = np.flatnonzero(np.diff(np.pad(conflict.astype(np.int8), (1, 1))))
        result.append(list(zip(edges[::2].tolist(), edges[1::2].tolist())))
    return result
