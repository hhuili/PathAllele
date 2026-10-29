import importlib.util
import os
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
    os.environ.get("PATHALLELE_ALIGNER"),
    "Native integration needs PATHALLELE_ALIGNER and calling dependencies",
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
                "--threads",
                "2",
                "--jobs",
                "2",
                "--aligner",
                os.environ["PATHALLELE_ALIGNER"],
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


if __name__ == "__main__":
    unittest.main()
