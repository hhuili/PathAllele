"""Reference-relative edits and haploid genotypes derived from complete alleles.

Coordinates are zero-based, half-open, relative to the unnormalized local REF.
VCF padding is a serialization detail and never participates in overlap tests.
"""

from dataclasses import dataclass

DNA = frozenset("ACGT")


@dataclass(frozen=True, order=True)
class Edit:
    start: int
    end: int
    alt: str
    complex: bool = False

    @property
    def structural(self):
        return self.complex or self.end - self.start != len(self.alt)


def trim_edit(start, ref, alt):
    while ref and alt and ref[0] == alt[0]:
        start += 1
        ref, alt = ref[1:], alt[1:]
    while ref and alt and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    return Edit(start, start + len(ref), alt)


def atomic_edits(start, ref, alt, keep_complex=False):
    edit = trim_edit(start, ref, alt)
    if edit.start == edit.end and not edit.alt:
        return []
    if edit.end - edit.start == len(edit.alt) == 1:
        return [Edit(edit.start, edit.end, edit.alt)]
    if keep_complex:
        return [Edit(edit.start, edit.end, edit.alt, True)]
    # All length-preserving substitutions use the same nucleotide-level units.
    if edit.end - edit.start == len(edit.alt):
        offset = edit.start - start
        return [
            Edit(edit.start + i, edit.start + i + 1, a)
            for i, a in enumerate(edit.alt)
            if a != ref[offset + i]
        ]
    return [edit]


def replay(ref, edits):
    cursor, parts = 0, []
    previous_insertion = None
    for e in sorted(set(edits)):
        if e.start < cursor or e.start < 0 or e.end > len(ref):
            raise ValueError(f"Overlapping/out-of-bounds edits: {e}, cursor={cursor}")
        if e.start == e.end and previous_insertion == e.start:
            raise ValueError("Multiple insertions at one boundary have no defined order")
        parts.extend((ref[cursor : e.start], e.alt))
        cursor = e.end
        previous_insertion = e.start if e.start == e.end else None
    return "".join(parts) + ref[cursor:]


def intersects(a, b):
    """Unpadded reference footprints, including insertion boundaries."""
    if a.start == a.end:
        if b.start == b.end:
            return a.start == b.start
        return b.start <= a.start <= b.end
    if b.start == b.end:
        return a.start <= b.start <= a.end
    return a.start < b.end and b.start < a.end


def vcf_alleles(reference, parent_pos0, start, end, alts, kind):
    s, e = parent_pos0 + start, parent_pos0 + end
    if kind == "SNP":
        return s + 1, reference[s:e], alts
    if s > 0:
        anchor = reference[s - 1]
        return s, reference[s - 1 : e], [anchor + a for a in alts]
    if e < len(reference):
        anchor = reference[e]
        return 1, reference[s : e + 1], [a + anchor for a in alts]
    raise ValueError("Cannot pad an event covering the entire reference")
