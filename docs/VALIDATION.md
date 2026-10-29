# Validation contract

A successful run means all implemented invariants pass; it does not mean every biological variant is assessable.

* Every selected native path base belongs to exactly one MATCH, bounded ALLELE/PARTIAL_SEQUENCE, or UNBOUNDED interval, in original query coordinates.
* Every bounded REF/ALT sequence is checked against its actual graph-path occurrence, including reverse orientation. Carrier identities agree with these observations.
* Every normalized source decomposition replays the complete source allele, including reference sequence needed after shifting an indel.
* Every final ALT genotype has a sequence-derived model witness; its serialization agrees with the positive projection model. Confident genotypes require resolved coverage.
* VCF sample order, reference bases, sorting, record/allele uniqueness, haploid GT indices and AC/AN counts are checked. Unknown sequence does not become confident ALT or reference calls.
* Full source-to-final haplotype completeness is **not** implied. Some models remain unrepresented because of competing copies, interacting events, missing flanks or unresolved homology. Their counts and FT reasons are retained. Synthetic integration tests additionally compare final reconstructed complete sequences with independent known truth.

The reproducible synthetic suite covers cross-PAF false bridges, supported within-PAF long deletion, graph-only tail recovery, true unbounded ends, physical overlap, opposite graph steps, unknown sequence, overlapping-word uniqueness, reverse CIGAR coordinates, invalid graph inputs, SNV/INDEL reconstruction, an all-reference empty callset, and restart integrity. Real-cohort diagnostics are maintained outside the software package; no diagnostic target list enters production calling.

When reporting empirical performance, use an independent truth set with a declared callable region and matching representation policy. Allele frequency and catalogue overlap alone are not sensitivity or specificity estimates. Establish organism/cohort-specific validity before extending beyond the tested scope.
