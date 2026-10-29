# Compact VCF schema (version 1)

The output follows VCF 4.2, uses explicit ACGT REF/ALT, a single haploid GT, sorted coordinates, BGZF compression and a tabix index. Missing calls stay `.`. `QUAL` and record `FILTER` stay `.` because no probabilistic quality or additional site filter is claimed.

| Field | Number | Meaning |
| --- | --- | --- |
| TYPE | A | SNV (1→1), MNV (equal lengths >1), INS (0→n), DEL (n→0), REPLACEMENT (unequal nonzero lengths), after trimming common padding |
| SIZECLASS | A | SV if max(trimmed REF length, trimmed ALT length) ≥50 bp, otherwise SMALL |
| SVLEN | A | Signed ALT minus REF length; zero does not exclude an SV-sized MNV |
| NATIVE | A | Native orientation evidence: SEQUENCE, INV, COMPLEX or combined support; separate from sequence-edit type |
| PID | . | Supporting source contexts in Graph.paths.tsv.gz / Graph.alleles.tsv.gz; record-level provenance, not graph parents |
| AC, AN, AF | A, 1, A | Observed haploid allele counts, called allele denominator and frequency |
| GT, FT | 1, 1 | Haploid allele index and genotype status; resolved calls PASS, unresolved calls retain their specific reason |

Per-ALT arrays follow exactly the ALT order, including records mixing edit types or sizes. `SV` is a size class; it is not a mutually exclusive edit mechanism. Small INS/DEL/REPLACEMENT may be grouped as INDEL by consumers. Large MNV/INS/DEL/REPLACEMENT belong to the SV size class. Adjacent SNVs or separate DEL/INS components are not automatically merged into a compound haplotype allele.

No repeated SOURCE/NBASES/NS, derivable REFLEN/ALTLEN/VCLASS or per-record UNASSESSED hierarchy fields are written. Sample-level FT is retained because missingness reasons are scientific information. Detailed source sequences and physical intervals remain in compressed sidecar tables. Graph hierarchy is explicitly unassessed once in `pathallele_hierarchy`; an independent native nested catalogue must carry its own coordinates/parents.

Schema validation recomputes TYPE and SIZECLASS from every REF/ALT, checks Number=A, reference agreement, haploid GT/FT, AC/AN, sample order, IDs and sorting, and replays source models. Format reference: https://samtools.github.io/hts-specs/VCFv4.2.pdf .
