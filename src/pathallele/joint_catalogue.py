"""Complete sequence catalogue and joint genotypes, before VCF factorization.

The VCF can represent a complete observed haplotype with several independent
records. A reference backbone at one record is allowed only when all residual
edits have already been emitted and reconstruct that observed sequence exactly.
No missing sequence is turned into reference evidence.
"""

import csv
import gzip
from pathlib import Path
import numpy as np

from .allele_model import Edit, replay


def factor_state(reference, start, end, observed, local, emitted):
    """Project a jointly observed sequence onto one compact event record."""
    if observed is None:
        raise ValueError("Factorization requires a resolved observed sequence")
    ref = reference[start:end]
    if observed == ref:
        return ref, "PASS"
    # One footprint may contain both substitution and length-changing alleles.
    # A component already emitted at a more specific record must not be called
    # again merely because this locus has a differently padded REF.
    if local and all(edit in emitted for edit in local):
        actual = replay(
            ref, [Edit(x.start - start, x.end - start, x.alt, x.complex) for x in local]
        )
        if actual != observed:
            raise ValueError("Independent VCF components changed the observed allele")
        return ref, "PASS"
    if len(local) == 1 and (local[0].start, local[0].end, local[0].alt) == (start, end, observed):
        return observed, "PASS"
    return None, "UNRESOLVED_COMPONENT"


class JointCatalogue:
    def __init__(self, samples, output):
        self.samples = samples
        self.output = Path(output)
        self.ids, self.bounds, self.genotypes, self.statuses, self.exports = [], [], [], [], []
        self.kinds, self.projected = [], []
        self.status_codes = {"PASS": 0}
        self.stream = gzip.open(self.output / "Graph.loci.tsv.gz", "wt")
        self.writer = csv.writer(self.stream, delimiter="\t", lineterminator="\n")
        self.writer.writerow(["locus_id", "allele_index", "ref_start0", "ref_end0", "sequence"])

    def codes(self, values):
        for value in values:
            if value not in self.status_codes:
                self.status_codes[value] = len(self.status_codes)
        return np.array([self.status_codes[value] for value in values], dtype=np.uint8)

    def add(self, ident, start, end, reference, observations, exports, kind="EDIT"):
        alleles = [reference] + sorted(
            {seq for seq, ft in observations if seq is not None and seq != reference},
            key=lambda seq: (len(seq), seq),
        )
        index = {seq: i for i, seq in enumerate(alleles)}
        for ai, sequence in enumerate(alleles):
            self.writer.writerow([ident, ai, start, end, sequence])
        self.ids.append(ident)
        self.bounds.append((start, end))
        self.kinds.append(kind)
        self.genotypes.append(
            np.array(
                [index[seq] if seq is not None else -1 for seq, ft in observations], dtype=np.int32
            )
        )
        self.statuses.append(self.codes([ft for seq, ft in observations]))
        self.exports.append(self.codes([ft for seq, ft in exports]))
        self.projected.append(
            np.array([index[seq] if seq is not None else -1 for seq, ft in exports], dtype=np.int32)
        )

    def close(self):
        self.stream.close()
        shape = (0, len(self.samples))
        np.savez_compressed(
            self.output / "Graph.locus_genotypes.npz",
            samples=np.array(self.samples),
            loci=np.array(self.ids),
            kinds=np.array(self.kinds),
            bounds=np.array(self.bounds, dtype=np.int64).reshape(-1, 2),
            gt=np.stack(self.genotypes) if self.genotypes else np.empty(shape, dtype=np.int32),
            projected_gt=np.stack(self.projected)
            if self.projected
            else np.empty(shape, dtype=np.int32),
            source_status=np.stack(self.statuses)
            if self.statuses
            else np.empty(shape, dtype=np.uint8),
            export_status=np.stack(self.exports)
            if self.exports
            else np.empty(shape, dtype=np.uint8),
            status_labels=np.array(list(self.status_codes)),
        )
