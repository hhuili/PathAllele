# PathAllele

PathAllele calls haploid variants from assembled pangenome graph paths. It retains
complete observed alleles, resolves their physical source intervals and exports
exact variants and structural allele groups as indexed VCFs.

## Getting started

Build with Rust >=1.82. The native calling dependencies are pinned below; the
complete workflow is tested on Linux.

```sh
git clone https://github.com/hhuili/PathAllele.git
cd PathAllele
conda env create -f environments/calling.yml
conda activate pathallele
cargo build --release --locked
export PATH="$PWD/target/release:$PATH"
bash environments/build_native.sh

# Call the bundled known-sequence example
pathallele call -g examples/demo/graph.gfa -r examples/demo/reference.fa \
    -s examples/demo/samples.txt -o calls --aligner tools/pathallele-align -t 2 -j 2
pathallele validate calls

# Inspect options and output formats
pathallele call -h
man ./pathallele.1
```

The caller and validators are Rust. Mapping uses minimap2; VCF processing uses
bcftools and samtools; a small native adapter calls vcflib/BiWFA. See
[THIRD_PARTY.md](THIRD_PARTY.md) for dependency notices.

## Usage

```sh
pathallele call -g graph.gfa -r ref.fa -s samples.txt -o out \
    --aligner pathallele-align -t 16 -j 8
pathallele validate out --outputs-only
pathallele export out -r ref.fa -o loci.vcf -l v2383
```

Input is a GFA1 graph with explicit `S` sequences, `P` paths and `0M` links, a
single A/C/G/T reference FASTA matching a graph path, and one non-reference sample
ID per line. Paths normally use `sample#haplotype#contig` names; `--path-map`
accepts a TSV with `path` and `sample` columns. Multiple contigs per sample are
allowed. Use a fresh output directory, or repeat an identical run with `--resume`.

## Output

| File | Contents |
| --- | --- |
| `Graph.vcf.gz` | Exact event alleles; haploid GT and missingness reason FT. |
| `Graph.loci.vcf.gz` | Complete-locus GT/FT and event-projection PG/PF, including reference-only loci. |
| `Graph.loci.tsv.gz` | Complete REF/ALT sequence dictionary. |
| `Graph.sv.vcf.gz` | Structural-group presence/absence/missing calls. |
| `Graph.sv_members.tsv.gz` | Exact complete/event alleles assigned to each structural group. |
| `Graph.paths.tsv.gz`, `Graph.alleles.tsv.gz` | Physical source intervals and allele evidence. |
| `run.json`, `validation.json` | Input/tool identities, execution metrics and validation results. |

All three VCFs are BGZF-compressed and indexed. `TYPE` distinguishes SNV, MNV,
INS, DEL, REPLACEMENT, INV, INV_COMPLEX and COMPLEX; `SEQTYPE` records sequence
shape separately. `SIZECLASS=SV` uses a 50-bp threshold. Structural groups allow
compatible sequence subtypes under a fixed 95% similarity rule.

Use complete-locus GT for sequence identity. Event GT=0 can contain changes
represented at linked records; overlapping complete loci are separate sequence
views. Structural groups are analysis units, not inferred ancestral events.
Unresolved observations remain missing. File semantics and grouping rules are
in [pathallele.1](pathallele.1).

Keep `work/` for full source validation and resume. `validate --outputs-only`
checks the portable output bundle without the original graph. Supported inputs
are single-reference haploid assemblies; diploid/polyploid samples, GFA `W` walks,
segment overlaps and circular reference coordinates are not supported. Source
consistency is not an independent estimate of biological accuracy.

## Development

```sh
cargo test --locked
cargo clippy --locked --all-targets -- -D warnings
PATHALLELE_ALIGNER="$PWD/tools/pathallele-align" \
    cargo test --release --locked --test integration -- --ignored --test-threads=1
```

Report reproducible problems through GitHub Issues. Cite the version used
([CITATION.cff](CITATION.cff)). MIT license.
