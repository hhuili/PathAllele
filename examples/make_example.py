#!/usr/bin/env python3
"""Write a deterministic synthetic GFA and independent SNV/indel truth table."""

import argparse
from pathlib import Path
import random


def make_example(output, all_reference=False):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    rng = random.Random(71029)
    reference = "".join(rng.choice("ACGT") for _ in range(6000))
    cuts = [0, 900, 901, 2000, 3000, 3003, 4000, 4010, 6000]
    nodes = {f"ref_{i}": reference[s:e] for i, (s, e) in enumerate(zip(cuts, cuts[1:]))}
    keys = list(nodes)
    refwalk = [(n, "+") for n in keys]
    changed = next(b for b in "ACGT" if b != reference[900])
    nodes.update(snv=changed, insertion="GGACT", unknown="NNNNNNNNNN")
    altwalk = [("snv", "+") if i == 1 else (n, "+") for i, n in enumerate(keys) if i != 4]
    altwalk.insert(3, ("insertion", "+"))
    walks = {
        "ref#0#chrom": refwalk,
        "sample0#1#contig": refwalk,
        "sample1#1#contig": refwalk if all_reference else altwalk,
        "sample2#1#contig": refwalk if all_reference else altwalk,
        "reverse#1#contig": [
            (n, "-" if sign == "+" else "+")
            for n, sign in reversed(refwalk if all_reference else altwalk)
        ],
    }
    links = {(*a, *b) for walk in walks.values() for a, b in zip(walk, walk[1:])}
    with (output / "graph.gfa").open("w") as f:
        f.write("H\tVN:Z:1.0\n")
        for name, seq in nodes.items():
            f.write(f"S\t{name}\t{seq}\n")
        for a, x, b, y in sorted(links):
            f.write(f"L\t{a}\t{x}\t{b}\t{y}\t0M\n")
        for name, walk in walks.items():
            f.write("P\t" + name + "\t" + ",".join(n + s for n, s in walk) + "\t*\n")
    (output / "reference.fa").write_text(">ref\n" + reference + "\n")
    (output / "samples.txt").write_text("sample0\nsample1\nsample2\nreverse\n")
    (output / "truth.tsv").write_text(
        "start0\tend0\talt\n"
        + ("" if all_reference else f"900\t901\t{changed}\n2000\t2000\tGGACT\n3000\t3003\t\n")
    )
    return output


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output")
    parser.add_argument("--all-reference", action="store_true")
    args = parser.parse_args()
    make_example(args.output, args.all_reference)
