"""Allele-wise sequence classification, independent of event size and graph hierarchy.

Lengths exclude common VCF padding. These describe edit representations, not
mutational mechanisms; the NATIVE field independently records orientation evidence.
"""

from dataclasses import dataclass

from .allele_model import trim_edit

SV_MIN_BP = 50
SCHEMA_VERSION = 1


@dataclass(frozen=True)
class AlleleClass:
    type: str
    sizeclass: str
    vclass: str
    reflen: int
    altlen: int


def classify_allele(ref, alt):
    """Classify one explicit sequence ALT, even in a mixed multiallelic record."""
    if not ref or not alt or set(ref + alt) - set("ACGT"):
        raise ValueError("Classification requires nonempty, explicit ACGT VCF alleles")
    edit = trim_edit(0, ref, alt)
    r, a = edit.end - edit.start, len(edit.alt)
    if r == a == 0:
        raise ValueError("ALT is identical to REF")
    kind = (
        "INS"
        if r == 0
        else "DEL"
        if a == 0
        else "SNV"
        if r == a == 1
        else "MNV"
        if r == a
        else "REPLACEMENT"
    )
    size = "SV" if max(r, a) >= SV_MIN_BP else "SMALL"
    group = "SV" if size == "SV" else kind if kind in ("SNV", "MNV") else "INDEL"
    return AlleleClass(kind, size, group, r, a)


def classification_info(ref, alts):
    rows = [classify_allele(ref, alt) for alt in alts]
    return {
        "TYPE": tuple(x.type for x in rows),
        "SIZECLASS": tuple(x.sizeclass for x in rows),
    }


VCF_HEADERS = [
    f"##pathallele_variant_schema={SCHEMA_VERSION}",
    '##pathallele_hierarchy="Unassessed: flat reference projection does not assign nested parents; PID links source contexts in Graph.paths.tsv.gz and Graph.alleles.tsv.gz"',
    '##INFO=<ID=TYPE,Number=A,Type=String,Description="Sequence edit type after trimming common padding: SNV, MNV, INS, DEL or REPLACEMENT; not a mutational mechanism">',
    f'##INFO=<ID=SIZECLASS,Number=A,Type=String,Description="SV if the larger REF or ALT edit length after trimming common padding is at least {SV_MIN_BP} bp, otherwise SMALL; independent of TYPE and net length change">',
]
