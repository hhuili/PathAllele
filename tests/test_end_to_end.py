import importlib.util
import shutil
import os
import random
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import pysam
from pathallele.graph_io import Graph
from pathallele.allele_model import replay, trim_edit

spec = importlib.util.spec_from_file_location(
    "example", Path(__file__).parents[1] / "examples/make_example.py"
)
example = importlib.util.module_from_spec(spec)
spec.loader.exec_module(example)


class InputTests(unittest.TestCase):
    def test_noninteger_ids_and_reference_integrity(self):
        with tempfile.TemporaryDirectory() as temp:
            data = example.make_example(temp)
            graph = Graph(data / "graph.gfa", data / "reference.fa", data / "samples.txt")
            self.assertEqual(len(graph.samples), 4)
            (data / "reference.fa").write_text(">ref\nACGT\n")
            with self.assertRaisesRegex(ValueError, "differs"):
                Graph(data / "graph.gfa", data / "reference.fa", data / "samples.txt")

    def test_missing_edge_and_duplicate_samples_fail(self):
        with tempfile.TemporaryDirectory() as temp:
            data = example.make_example(temp)
            original = (data / "graph.gfa").read_text()
            (data / "graph.gfa").write_text(
                "\n".join(l for l in original.splitlines() if not l.startswith("L\t"))
            )
            with self.assertRaisesRegex(ValueError, "edge"):
                Graph(data / "graph.gfa", data / "reference.fa", data / "samples.txt")
            (data / "graph.gfa").write_text(original)
            (data / "samples.txt").write_text("sample0\nsample0\n")
            with self.assertRaisesRegex(ValueError, "unique"):
                Graph(data / "graph.gfa", data / "reference.fa", data / "samples.txt")


@unittest.skipUnless(
    os.environ.get("PATHALLELE_ALIGNER")
    and all(shutil.which(tool) for tool in ("minimap2", "bcftools", "samtools")),
    "Integration needs PATHALLELE_ALIGNER, minimap2, bcftools and samtools",
)
class IntegrationTests(unittest.TestCase):
    def run_fixture(self, all_reference=False):
        with tempfile.TemporaryDirectory(prefix="pathallele test ") as temp:
            root = Path(temp)
            data = example.make_example(root / "inputs with spaces", all_reference)
            output = root / "result with spaces"
            command = [
                sys.executable,
                "-m",
                "pathallele.cli",
                "call",
                "--graph",
                str(data / "graph.gfa"),
                "--reference",
                str(data / "reference.fa"),
                "--samples",
                str(data / "samples.txt"),
                "--output",
                str(output),
                "--aligner",
                os.environ["PATHALLELE_ALIGNER"],
                "--threads",
                "2",
            ]
            run = subprocess.run(command, text=True, capture_output=True)
            if run.returncode:
                self.fail(
                    run.stderr
                    + "\n"
                    + ((output / "run.log").read_text() if (output / "run.log").exists() else "")
                )
            reference = "".join((data / "reference.fa").read_text().splitlines()[1:])
            graph = Graph(data / "graph.gfa", data / "reference.fa", data / "samples.txt")
            edits = {s: [] for s in graph.samples}
            count = 0
            with pysam.VariantFile(output / "Graph.vcf.gz") as vcf:
                for record in vcf:
                    count += 1
                    for s in graph.samples:
                        gt = record.samples[s]["GT"][0]
                        self.assertIsNotNone(gt)
                        if gt:
                            edits[s].append(
                                trim_edit(record.start, record.ref, record.alts[gt - 1])
                            )
            self.assertEqual(count, 0 if all_reference else 3)
            from pathallele.graph_io import reverse_complement

            for s in graph.samples:
                expected = graph.spell(graph.paths[s][0])
                if s == "reverse":
                    expected = reverse_complement(expected)
                self.assertEqual(replay(reference, edits[s]), expected)
            repeated = subprocess.run(command, text=True, capture_output=True)
            self.assertEqual(repeated.returncode, 2)
            resumed = subprocess.run(command + ["--resume"], text=True, capture_output=True)
            self.assertEqual(resumed.returncode, 0, resumed.stderr)
            changed = subprocess.run(
                command + ["--resume", "--minimum-anchor", "35"], text=True, capture_output=True
            )
            self.assertEqual(changed.returncode, 2)
            with (output / "Graph.vcf.gz").open("ab") as stream:
                stream.write(b"corrupt")
            corrupt = subprocess.run(command + ["--resume"], text=True, capture_output=True)
            self.assertEqual(corrupt.returncode, 2)

    def test_snv_insertion_deletion_and_reverse_paths(self):
        self.run_fixture()

    def test_all_reference_empty_catalogue(self):
        self.run_fixture(True)

    def test_complete_alleles_and_independent_snv_backbones(self):
        with tempfile.TemporaryDirectory(prefix="pathallele joint ") as temp:
            root = Path(temp)
            rng = random.Random(912)
            flank = lambda: "".join(rng.choice("ACGT") for _ in range(2000))
            left, right = flank(), flank()
            ref, alt = "A" * 60 + "CG", "CG" + "T" * 59 + "C"
            snv = ref[:30] + "G" + ref[31:]
            partial = ref[:10] + "N" * 10 + ref[20:]
            nodes = dict(left=left, right=right, ref=ref, alt=alt, snv=snv, partial=partial)
            paths = {"ref#0#chrom": ["left+", "ref+", "right+"]}
            for name, mid in [
                ("same", "ref"),
                ("changed", "alt"),
                ("small", "snv"),
                ("unknown", "partial"),
            ]:
                paths[name + "#1#contig"] = ["left+", mid + "+", "right+"]
            paths["reverse#1#contig"] = ["right-", "alt-", "left-"]
            links = {(a, b) for walk in paths.values() for a, b in zip(walk, walk[1:])}
            lines = ["H\tVN:Z:1.0"] + [f"S\t{n}\t{s}" for n, s in nodes.items()]
            lines += [f"L\t{a[:-1]}\t{a[-1]}\t{b[:-1]}\t{b[-1]}\t0M" for a, b in sorted(links)]
            lines += ["P\t" + n + "\t" + ",".join(w) + "\t*" for n, w in paths.items()]
            (root / "graph.gfa").write_text("\n".join(lines) + "\n")
            (root / "reference.fa").write_text(">ref\n" + left + ref + right + "\n")
            samples = ["same", "changed", "small", "unknown", "reverse"]
            (root / "samples.txt").write_text("\n".join(samples) + "\n")
            out = root / "output"
            cmd = [
                sys.executable,
                "-m",
                "pathallele.cli",
                "call",
                "--graph",
                str(root / "graph.gfa"),
                "--reference",
                str(root / "reference.fa"),
                "--samples",
                str(root / "samples.txt"),
                "--output",
                str(out),
                "--aligner",
                os.environ["PATHALLELE_ALIGNER"],
                "--threads",
                "2",
            ]
            result = subprocess.run(cmd, text=True, capture_output=True)
            self.assertEqual(
                result.returncode, 0, result.stderr + "\n" + (out / "run.log").read_text()
            )
            edits = {s: [] for s in samples}
            with pysam.VariantFile(out / "Graph.vcf.gz") as v:
                records = [r.copy() for r in v]
            compound = next(r for r in records if "MNV" in r.info["TYPE"])
            self.assertEqual(compound.samples["small"]["GT"], (0,))
            self.assertEqual(compound.samples["unknown"]["GT"], (None,))
            for r in records:
                for s in samples:
                    gt = r.samples[s]["GT"][0]
                    if gt:
                        edits[s].append(trim_edit(r.start, r.ref, r.alts[gt - 1]))
            for sample, sequence in [
                ("same", ref),
                ("small", snv),
                ("changed", alt),
                ("reverse", alt),
            ]:
                self.assertEqual(replay(left + ref + right, edits[sample]), left + sequence + right)
            import numpy as np

            joint = np.load(out / "Graph.locus_genotypes.npz")
            i = list(joint["loci"]).index(compound.id)
            self.assertGreater(joint["gt"][i, 2], 0)
            self.assertEqual(joint["projected_gt"][i, 2], 0)
            self.assertEqual(joint["gt"][i, 3], -1)


if __name__ == "__main__":
    unittest.main()
