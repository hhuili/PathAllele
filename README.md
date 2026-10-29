# PathAllele

**Source-resolved haploid variants from pangenome graph paths.**

PathAllele reads sequence-resolved GFA1 paths and a reference FASTA, extracts physically bounded alleles, and projects a joint VCF with explicit missing-genotype reasons. Every emitted ALT is checked against a sequence-derived source model. The caller uses minimap2 to resolve physical occurrence context. Graph identity and unique internal sequence anchors define complete allele units before variant classification; a gapped CIGAR does not determine their boundaries or edit types. A separate, direction-checked BiWFA projection retains assessable small changes when a physical source allele crosses unresolved reference territory. It does not run `vg deconstruct` or `vcfwave`, and is not an implementation of the HPRC calling workflow.

**Status: 0.2.0rc1, research pre-release.** The supported scope is one reference, haploid samples, and compact genomes. Synthetic tests establish specified algorithmic behavior; biological sensitivity and specificity require independent validation. The current dense projection uses memory proportional to sample count × reference length. Whole human cohorts, diploid genotypes, GFA W records, unresolved segment sequences, and overlapping graph segments are outside this release's scope.

## Install on Linux

```bash
git clone https://github.com/hhuili/PathAllele.git
cd PathAllele
conda env create -f environments/calling.yml
conda activate pathallele
python -m pip install .
bash environments/build_native.sh
pathallele doctor --aligner "$PWD/tools/pathallele-align"
```

The small BiWFA adapter is used for local projection in partially assessable source contexts. It does not define the complete allele catalogue or classify MNV/replacement alleles. The build helper verifies pinned WFA2 headers; dependencies are separately installed. The full calling workflow is tested on Linux; unit tests also run on macOS.

## Synthetic example

```bash
python examples/make_example.py demo/input
pathallele call --graph demo/input/graph.gfa --reference demo/input/reference.fa \
  --samples demo/input/samples.txt --output demo/output \
  --aligner "$PWD/tools/pathallele-align" --threads 2 --jobs 2
pathallele validate demo/output
```

The example has independent known SNV, insertion and deletion sequences, an unchanged sample, and a reverse-oriented path. The integration test reconstructs each complete sample sequence from the final VCF and compares it with the known input path.

## Real inputs

```bash
pathallele call --graph cohort.gfa --reference reference.fa --samples samples.txt \
  --output calls --aligner /path/to/pathallele-align --threads 16 --jobs 8
```

* GFA1 `S`, `L`, and `P` records must contain explicit sequences and zero-overlap links. Arbitrary segment identifiers are supported. Each selected path edge and the spelled reference are validated.
* `reference.fa` contains exactly one fully resolved A/C/G/T sequence. Select a specific reference path with `--reference-path` if needed.
* `samples.txt` lists distinct, non-reference haploid sample IDs, one per line. By default a PanSN path `sample#haplotype#contig` belongs to `sample`; multiple haplotypes for one sample are rejected. Multiple retained contigs are allowed.
* For other path names, provide `--path-map mapping.tsv` with `path` and `sample` columns. An explicit map declares haploid ownership; the user is responsible for assigning biological sample identities correctly.
* Sample paths may contain IUPAC ambiguity codes. Bounded unresolved sequence is retained as source evidence; it cannot produce a confident genotype through those bases.

`--sequence-anchor` (unique internal homology within an accepted source interval) defaults to 31 bp. `--minimum-anchor` defaults to 31 bp, `--mapping-k` to 19 bp, and `--minimum-mapq` to 20. The graph-only fallback permits gaps up to 1,000 bp and requires unique anchors in both orientations. This bound applies to unmapped graph contexts, not deletions supported inside one primary alignment. See [the algorithm](docs/ALGORITHM.md) for their distinct roles.

## Outputs and restart

| File | Meaning |
| --- | --- |
| `Graph.vcf.gz` / `.tbi` | Joint, source-bounded normalized haploid VCF; GT/FT; AC/AN/AF; per-ALT TYPE and SIZECLASS |
| `Graph.paths.tsv.gz` | Physical path intervals, alignment component, strand, source class and source allele |
| `Graph.alleles.tsv.gz` | Complete bounded source sequences |
| `Graph.anchors.tsv.gz` | Source-relative graph and unique sequence homology anchors |
| `Graph.loci.tsv.gz` | Complete observed REF/ALT sequence catalogue at each joint locus |
| `Graph.locus_genotypes.npz` | Actual joint sequence indices, compact projection indices and their separate missingness reasons |
| `Graph.alignment_ambiguity.tsv.gz` | Direction-dependent unresolved source positions |
| `summary.*.json` | Source and projection counts |
| `run.json`, `validation.json`, `run.log` | Input/code/tool hashes, stage status, invariants, diagnostics |
| `work/` | Auditable mapping and alignment intermediates; required for full validation and resume |

Use a fresh output directory. To resume an interrupted run, repeat the identical command with `--resume`. Completed stages are reused only after checking input, code, options, executable identity and stage-output hashes. A changed dependency or input is an error, not an implicit rerun. An output lock prevents concurrent writers. Failed stage evidence remains available in the log; that stage's own temporary projection work is replaced on an authorized resume.

Unknown coverage, competing physical copies, ambiguous alignment, overlapping edits and partial sequence remain missing. Complete mixed units retain their full sequence; separate graph-anchored SNVs/INDELs remain independent. Compact VCF records can overlap in coordinate span, but called ALT edits cannot overlap within a sample. A reference backbone at a compound record is allowed only when other explicitly called edits reconstruct that sample's complete observed allele. The sidecar always retains the actual complete sequence and joint genotype; a nonreference haplotype is never erased from that catalogue. Exact catalogue counts alone therefore do not measure detection accuracy. See [validation scope](docs/VALIDATION.md).

## Variant classification

`TYPE` is per ALT (`Number=A`): `SNV`, `MNV`, `INS`, `DEL`, or `REPLACEMENT`, after removing common padding. `SIZECLASS` independently records `SMALL` or `SV` (larger edited REF/ALT length ≥50 bp). A 50-bp MNV is therefore `TYPE=MNV;SIZECLASS=SV`; net length change alone does not define size.

The compact VCF retains `SVLEN`, `NATIVE` orientation evidence, and `PID` links to source-context tables, alongside standard counts and genotypes. Derivable lengths/classes and repeated constant fields are omitted. Uncompressed serialization stays in `work/`; the deliverable is BGZF-compressed and tabix-indexed. See [the VCF schema](docs/VCF.md).

Graph hierarchy is a separate dimension: this reference-projected caller does not assign nested parents. This limitation is recorded once in the header. `PID` does not mean a biological parent variant. Independent graph-nested catalogues can overlap sequence carried inside a VCF ALT; do not add their counts without occurrence-aware deduplication.

## Development

```bash
python -m pip install -e .
python -m unittest discover -s tests -v
```

Set `PATHALLELE_ALIGNER` to the compiled adapter to run integration tests with minimap2, bcftools and samtools. CI tests installation, unit behavior, full sequence reconstruction, complete joint alleles, missingness and empty-catalogue handling. See [contribution guidance](CONTRIBUTING.md). Original code is [MIT licensed](LICENSE); dependencies retain their [own notices](THIRD_PARTY.md). Cite the software version and commit with [CITATION.cff](CITATION.cff).
