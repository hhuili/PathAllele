# Implementation-to-figure map

Revised 2026-11-11; implementation baseline `057d45a9`.

| Figure element | Implementation | Meaning |
|---|---|---|
| Graph paths | `src/graph.rs::Graph::read/spell` | Read existing oriented walks and retain path offsets; no exhaustive graph traversal |
| minimap2 stream | `src/source.rs::build`, `src/context.rs::read_alignments` | Map whole path sequences to the matching reference; accept primary records at the MAPQ threshold |
| Node visits at r/q | `src/graph.rs::ref_occurrences`, `src/anchors.rs::occurrence_seeds` | Distinguish a graph node from one of its path occurrences; check oriented node and mapped reference position |
| Evidence convergence | `src/context.rs::plan`, `src/anchors.rs::context_runs/obligatory` | Source/uniqueness filters and anchors common to all optimal eligible collinear chains |
| Complete REF/ALT | `src/source.rs::Collector::accept`, `src/alleles.rs::prepare` | Keep bounded complete sequence and physical source evidence |
| Crossed RC blocks / INV_COMPLEX | `src/classify.rs::Classifier::annotate`, `Annotation::assign` | Generic local reverse-homology classification; separate from native source class |
| Allele calls | `src/catalogue.rs`, `src/project.rs` | Complete-sequence GT/FT, separate from event projection |
| Normalization / replay / event projection | `src/alleles.rs`, `src/variant.rs`, `src/project.rs` | Normalize edits, reconstruct observed sequence, preserve unresolved states |

The synthetic core is AGTCCCGA → TCGGGACT with TT→AC beside it. Input branches are forward-oriented; the classification stage shows sequence homology, not an asserted reverse graph traversal. Black flanks abstract longer eligible anchors. The GGA example is a short occurrence glyph, not a claim that a three-base node alone satisfies default anchoring thresholds.

Not drawn: FASTA/reference-path equality check; graph-only supplementation with uniqueness/orientation/gap restrictions; conditional integrated vcflib/WFA2; output housekeeping; SV allele clustering. These remain in the algorithm and legend. Actual oriP evidence is in `ORIP_CASE_STUDY.zh-CN.md`; the figure is not an oriP benchmark or a literal run on the eight-base drawing.
