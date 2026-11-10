# PathAllele

PathAllele calls haploid variants from assembled pangenome graph paths. It retains
complete observed alleles, resolves their physical source intervals and exports
exact variants and structural allele groups as indexed VCFs.

## Getting started

Install the Linux x86-64 development build from the personal Conda channel:

```sh
conda create -n pathallele --override-channels -c conda-forge -c hhuili -c bioconda pathallele=0.2.0.dev20261110
conda activate pathallele
pathallele -g graph.gfa -r ref.fa -o out
```

For a source build:

Build with Rust >=1.82. The native calling dependencies are pinned below; the
complete workflow is tested on Linux.

```sh
git clone https://github.com/hhuili/PathAllele.git
cd PathAllele
conda env create -f environments/calling.yml
conda activate pathallele
bash environments/build_native.sh integrated
export PATH="$PWD/target/release:$PATH"

# Call the bundled known-sequence example
pathallele -g examples/demo/graph.gfa -r examples/demo/reference.fa \
    -o calls -t 2 -j 2

# Inspect options
pathallele --help
```

The calling algorithm is implemented in Rust. Mapping uses minimap2; VCF processing uses
bcftools and samtools. The integrated build embeds the native vcflib/WFA2
adapter in the PathAllele executable. It uses isolated internal workers and
does not require a separate `pathallele-align` executable. Native libraries
are still runtime dependencies; this is not a fully static build. See
[THIRD_PARTY.md](THIRD_PARTY.md) for dependency notices.

## Usage

```sh
pathallele -g graph.gfa -r ref.fa -o out -t 16 -j 8
```

PathAllele has one public operation and no subcommands. Select or filter its
standard VCF outputs with bcftools. The calling engine is built into PathAllele.
There is no `--aligner` option,
external-adapter discovery or `PATHALLELE_ALIGNER` override. The executable and
loaded libraries remain recorded in `run.json`. Bioconda builds enable
`native-aligner`; compilers are build dependencies, not user prerequisites.

Calling requires the `native-aligner` build feature. Manual native
builds use `PATHALLELE_NATIVE_PREFIX` (or `CONDA_PREFIX`). The dependency environment
provides vcflib/WFA2 libraries and their matching headers; no third-party headers
are bundled in this repository. Use a C++ compiler compatible with the installed
vcflib ABI; `PATHALLELE_CXX11_ABI` defaults to `1` for the pinned Conda libraries.
The adapter embeds the existing vcflib/WFA2 backend, not a Rust reimplementation
of BiWFA. The isolated internal worker protects the main process from upstream
library exit paths.

Bioconda recipe preparation and local package tests are described in
[recipe/README.md](recipe/README.md).

Changing a binary or backend intentionally invalidates an existing run's strict
resume identity. Use a fresh output directory for comparisons.

Input is a GFA1 graph with explicit `S` sequences, `P` paths and `0M` links, a
single A/C/G/T reference FASTA matching a graph path, and optionally one non-reference sample
ID per line. Paths normally use `sample#haplotype#contig` names; `--path-map`
accepts a TSV with `path` and `sample` columns. Multiple contigs per sample are
allowed. Without `-s`, all non-reference samples are selected in first-path
order (or first-row order with `--path-map`). A plain path name is treated as
one sample ID. Use `-s` to select a subset or define its order.
Use a fresh output directory, or repeat an identical run with `--resume`.

## Output

| File | Contents |
| --- | --- |
| `Graph.vcf.gz` | Exact event alleles; haploid GT and missingness reason FT. |
| `Graph.loci.vcf.gz` | Complete-locus GT/FT and event-projection PG/PF, including reference-only loci. |
| `Graph.loci.tsv.gz` | Complete REF/ALT sequence dictionary. |
| `Graph.sv.vcf.gz` | Structural-group presence/absence/missing calls. |
| `Graph.sv_members.tsv.gz` | Exact complete/event alleles assigned to each structural group. |
| `Graph.paths.tsv.gz`, `Graph.alleles.tsv.gz` | Physical source intervals and allele evidence. |
| `run.json` | Input/tool identities, execution metrics and output identities for resume. |

All three VCFs are BGZF-compressed and indexed. `TYPE` distinguishes SNV, MNV,
INS, DEL, REPLACEMENT, INV, INV_COMPLEX and COMPLEX; `SEQTYPE` records sequence
shape separately. `SIZECLASS=SV` uses a 50-bp threshold. Structural groups allow
compatible sequence subtypes under a fixed 95% similarity rule.

Use complete-locus GT for sequence identity. Event GT=0 can contain changes
represented at linked records; overlapping complete loci are separate sequence
views. Structural groups are analysis units, not inferred ancestral events.
Unresolved observations remain missing. VCF positions are 1-based; TSV intervals
are 0-based and half-open. All genotypes are haploid; `.` means missing, and
`FT` explains the source status. `AC/AN/AF` exclude missing calls.

Complete-locus `GT` and `PG` index the sequence dictionary in `Graph.loci.tsv.gz`
(`locus_id`, `allele_index`, with 0 for REF), not the event VCF's ALT order.
Structural-group `GT` is 1 for presence, 0 for supported absence, or `.` for
insufficient/conflicting evidence. Absence requires complete sequence coverage
of the group; event `GT=0` alone is insufficient. In `Graph.sv_members.tsv.gz`,
`SEED` rows identify emitted SV alleles and `COMPLETE` rows identify additional
compatible complete alleles, which may lack an event allele index.

Allele classification uses trimmed sequences: `REPLACEMENT` has nonempty,
unequal-length REF and ALT, and `SIZECLASS=SV` requires a larger trimmed allele
of at least 50 bp. `INV` requires an exact same-interval reverse complement;
`INV_COMPLEX` retains an inverted component with additional changes or shifted
boundaries. Sequence-only inversion evidence requires a unique informative
reverse core of at least 31 bp.

Structural groups use complete-link clustering with compatible classes and
inversion axes, endpoint distances within 5% of the shorter span (bounded to
1–50 bp), and at least 95% reference-span, ALT-length and sequence similarity.
Sequences are compared in a common reference frame; deletions compare deleted
sequence. Approximate non-deletion groups must be distinguishable from the
reference beyond that distance budget. Insertions at different positions
require a closed common frame to merge. `COMPLEX` receives exact deduplication
only, and an additional complete allele must uniquely match every group member.

Keep `work/` for resume. Independent sequence fixtures and frozen-output comparisons run in the
integration tests; production calls do not repeat a whole-run audit. Supported inputs
are single-reference haploid assemblies; diploid/polyploid samples, GFA `W` walks,
segment overlaps and circular reference coordinates are not supported. Source
consistency is not an independent estimate of biological accuracy.

## Development

```sh
cargo test --locked
cargo clippy --locked --all-targets -- -D warnings
cargo test --release --locked --features native-aligner --test integration -- \
    --ignored --test-threads=1 --skip embedded_and_external_calls_have_identical_output_contracts
```

Report reproducible problems through GitHub Issues. Cite the version used
([CITATION.cff](CITATION.cff)). MIT license.
