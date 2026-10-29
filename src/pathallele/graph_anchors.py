"""Conservative native graph anchors and explicit alignment ambiguity."""


def unique_runs(ref, alt, anchors):
    from .path_context import context_runs

    # These anchors already have physical graph-occurrence and independent
    # whole-path alignment context from path_blocks, including short nodes.
    return [
        item for run in context_runs(ref, alt, anchors, sequence_unique=False) for item in run[3]
    ]


def anchored_gaps(ref, alt, anchors):
    """Return differing gaps between fixed graph anchors, with REF offsets."""
    assert ref[0] == alt[0] and ref[-1] == alt[-1]
    ordered = (
        [(0, 0, 1)] + [(r, q, n) for r, q, n, node in anchors] + [(len(ref) - 1, len(alt) - 1, 1)]
    )
    r_end = q_end = 0
    gaps = []
    for r, q, n in ordered:
        if r < r_end or q < q_end or n <= 0:
            raise ValueError("Overlapping/noncollinear native graph anchors")
        if ref[r : r + n] != alt[q : q + n] or len(ref[r : r + n]) != n:
            raise ValueError("Native graph anchor sequence mismatch")
        a, b = ref[r_end:r], alt[q_end:q]
        if a != b:
            gaps.append((r_end, a, b))
        r_end, q_end = r + n, q + n
    assert r_end == len(ref) and q_end == len(alt)
    return gaps


def aligned_bases(ref, edits):
    """Concrete reference-position bases; spanning structures are unresolved."""
    bases = list(ref)
    for e in edits:
        if e.structural:
            bases[e.start : e.end] = ["?"] * (e.end - e.start)
        else:
            bases[e.start : e.end] = list(e.alt)
    return bases
