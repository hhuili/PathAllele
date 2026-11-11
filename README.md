# PathAllele

PathAllele calls haploid variants from assembly pangenome graphs. It currently
targets small genomes, primarily EBV graphs built with PGGB. Other viral genomes
are intended applications; large eukaryotic pangenomes have not been benchmarked.

A path in a cyclic graph may revisit a node. PathAllele tracks each **occurrence**:
one visit to an oriented node on a path, with its step and sequence offset.
Mapping and path context resolve reference correspondence; node identity or visit
order alone is insufficient. The method combines occurrence evidence with complete
alleles, sequence replay and conservative missing calls. These latter operations
are also available in other assembly-based callers.
Within each complete source allele, local reverse-complement evidence constrains
decomposition boundaries before variant classification.

## Installation

Published Linux x86-64 build:

```sh
conda create -n pathallele --override-channels -c conda-forge -c bioconda -c hhuili pathallele=0.2.0.dev20261110
conda activate pathallele
```

To build from source, use Rust >=1.82 and the pinned native dependencies:

```sh
git clone https://github.com/hhuili/PathAllele.git
cd PathAllele
conda env create -f environments/calling.yml
conda activate pathallele
bash environments/build_native.sh integrated
export PATH="$PWD/target/release:$PATH"
```

## Usage

```sh
pathallele -g graph.gfa -r ref.fa -o calls -t 16 -j 8
pathallele -h
```

Input is GFA1 with explicit `S` sequences, `P` paths and `0M` links, plus one
A/C/G/T reference FASTA matching a complete graph path. Use `--reference-path`
when necessary. All non-reference samples are selected by default; use
`-s samples.txt` to select/order samples or `--path-map` for custom path names.
Multiple contigs per sample are allowed. Use a fresh output directory or
`--resume` for an identical run.

minimap2 maps path sequences to the reference; only primary alignments passing
MAPQ filtering are retained. Graph evidence provides bounded supplementation,
so missed or incorrect mappings can affect calls. The mapper has not been
established as optimal. Local projection uses integrated vcflib/WFA2;
normalization uses bcftools. See [dependency notices](THIRD_PARTY.md).

## Output

| File | Contents |
| --- | --- |
| `Graph.vcf.gz` | Exact events with haploid GT/FT. |
| `Graph.loci.vcf.gz`, `Graph.loci.tsv.gz` | Complete alleles, sequence GT/FT and event-projection PG/PF. |
| `Graph.sv.vcf.gz`, `Graph.sv_members.tsv.gz` | Structural allele clusters and membership. |
| `Graph.paths.tsv.gz`, `Graph.alleles.tsv.gz` | Source intervals and allele evidence. |
| `run.json` | Input, tool and output identities. |

Sequence GT and event GT describe different views; missing does not mean reference.
Replay checks sequence consistency, not biological accuracy. SV clusters are not
inferred ancestral mutations. Keep `work/` for resume.

The caller is single-reference and haploid. GFA `W` walks, segment overlaps,
diploid/polyploid genotyping and circular reference coordinates are unsupported.
Graph cycles and repeated-node visits are supported separately.

## Development

```sh
cargo test --locked
cargo clippy --locked --all-targets -- -D warnings
```

See [native integration tests](recipe/README.md). Cite the version used
([CITATION.cff](CITATION.cff)). MIT license.
