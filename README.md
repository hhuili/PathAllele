# PathAllele

**Source-resolved haploid variants from pangenome graph paths.**

PathAllele reads sequence-resolved GFA1 paths and a reference FASTA, extracts physically bounded alleles, and projects a joint VCF with explicit missing-genotype reasons. Every emitted ALT is checked against a sequence-derived source model. The caller uses minimap2 for occurrence context and vcflib's BiWFA engine for allele alignment. It does not run `vg deconstruct` or `vcfwave`, and is not an implementation of the HPRC calling workflow.

**Status: 0.1.0rc1, research pre-release.** The supported scope is one reference, haploid samples, and compact genomes. Synthetic tests establish specified algorithmic behavior; biological sensitivity and specificity require independent validation. The current dense projection uses memory proportional to sample count × reference length. Whole human cohorts, diploid genotypes, GFA W records, unresolved segment sequences, and overlapping graph segments are outside this release's scope.

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

The small native adapter must be built against the installed vcflib/WFA2 libraries. The build helper downloads a pinned, checksum-verified WFA2 header tree; it does not modify third-party sources. To supply existing dependencies manually, use `pathallele build-aligner --help`. Python-only tests also run on macOS; native calling is currently tested on Linux.

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

`--minimum-anchor` defaults to 31 bp, `--mapping-k` to 19 bp, and `--minimum-mapq` to 20. The graph-only fallback permits gaps up to 1,000 bp and requires unique anchors in both orientations. This bound applies to unmapped graph contexts, not deletions supported inside one primary alignment. See [the algorithm](docs/ALGORITHM.md) for their distinct roles.

## Outputs and restart

| File | Meaning |
| --- | --- |
| `Graph.vcf.gz` / `.tbi` | Joint, reference-normalized haploid VCF; GT/FT and AC/AN/AF/NS |
| `Graph.paths.tsv.gz` | Physical path intervals, alignment component, strand, source class and parent allele |
| `Graph.alleles.tsv.gz` | Complete bounded source sequences |
| `Graph.anchors.tsv.gz` | Source-relative retained homology anchors |
| `Graph.alignment_ambiguity.tsv.gz` | Direction-dependent unresolved source positions |
| `summary.*.json` | Source and projection counts |
| `run.json`, `validation.json`, `run.log` | Input/code/tool hashes, stage status, invariants, diagnostics |
| `work/` | Auditable mapping and alignment intermediates; required for full validation and resume |

Use a fresh output directory. To resume an interrupted run, repeat the identical command with `--resume`. Completed stages are reused only after checking input, code, options, executable identity and stage-output hashes. A changed dependency or input is an error, not an implicit rerun. An output lock prevents concurrent writers. Failed stage evidence remains available in the log; that stage's own temporary projection work is replaced on an authorized resume.

Unknown coverage, competing physical copies, ambiguous alignment, overlapping edits and partial sequence remain missing. Equal-length substitutions can be decomposed into SNVs; unstable compound alleles retain their full sequence. Exact catalogue counts alone therefore do not measure detection accuracy. See [validation scope](docs/VALIDATION.md).

## Development

```bash
python -m pip install -e .
python -m unittest discover -s tests -v
PATHALLELE_ALIGNER="$PWD/tools/pathallele-align" python -m unittest discover -s tests -v
```

The second command runs the native end-to-end tests as well. CI tests installation, unit behavior, native compilation, sequence reconstruction and empty-catalogue handling. See [contribution guidance](CONTRIBUTING.md). Original code is [MIT licensed](LICENSE); dependencies retain their [own notices](THIRD_PARTY.md). Cite the software version and commit with [CITATION.cff](CITATION.cff).
