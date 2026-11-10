#!/usr/bin/env bash
set -euo pipefail
! command -v pathallele-align >/dev/null 2>&1
pathallele --help > help.txt
example="$PREFIX/share/pathallele/examples/demo"
PATHALLELE_ALIGNER=/absent/obsolete-adapter pathallele \
  -g "$example/graph.gfa" -r "$example/reference.fa" \
  -o smoke-calls -t 2 -j 1 > smoke.log 2>&1
bcftools query -l smoke-calls/Graph.loci.vcf.gz > samples.txt
diff "$example/samples.txt" samples.txt
for name in Graph.vcf.gz Graph.loci.vcf.gz Graph.sv.vcf.gz; do
  bcftools view -Ov -o /dev/null "smoke-calls/$name"
  test -s "smoke-calls/$name.tbi"
done
