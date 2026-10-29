import random
import tempfile
from pathlib import Path
from types import SimpleNamespace
import unittest
from collections import defaultdict
from pathallele.components import (
    Alignment,
    graph_chains,
    plan_path,
    read_alignments,
    unique_oriented,
)
from pathallele.graph_io import offsets, reverse_complement
from pathallele.allele_model import Edit, replay, atomic_edits, vcf_alleles


def fake_graph(sequences, refwalk, querywalk):
    nodes = {i + 1: s for i, s in enumerate(sequences)}

    def spell(walk):
        return "".join(nodes[n] if n > 0 else reverse_complement(nodes[-n]) for n in walk)

    places = defaultdict(list)
    for i, n in enumerate(refwalk):
        places[abs(n)].append(i)
    ro = offsets(nodes, refwalk)
    return SimpleNamespace(
        nodes=nodes,
        refwalk=refwalk,
        ref_offsets=ro,
        ref_places=places,
        ref_occurrences={(n, ro[i]): i for i, n in enumerate(refwalk)},
        walks={"sample": querywalk},
        reference=spell(refwalk),
        spell=lambda name: spell(querywalk),
    )


class SourceTests(unittest.TestCase):
    def setUp(self):
        rng = random.Random(401)
        self.seq = ["".join(rng.choice("ACGT") for _ in range(n)) for n in (100, 2000, 100, 100, 1)]

    def test_separate_pafs_do_not_create_large_deletion(self):
        graph = fake_graph(self.seq, [1, 2, 3], [1, 3])
        records = [
            Alignment("A", "sample", "+", 0, 100, 0, 100, ((0, 100, 0),)),
            Alignment("B", "sample", "+", 100, 200, 2100, 2200, ((100, 200, 2100),)),
        ]
        blocks, missing, _ = plan_path(graph, "sample", records)
        self.assertEqual(len(blocks), 2)
        self.assertFalse(missing)
        self.assertTrue(all(len(b["runs"]) == 1 for b in blocks))

    def test_single_paf_preserves_supported_large_deletion(self):
        graph = fake_graph(self.seq, [1, 2, 3], [1, 3])
        record = Alignment("A", "sample", "+", 0, 200, 0, 2200, ((0, 100, 0), (100, 200, 2100)))
        blocks, missing, _ = plan_path(graph, "sample", [record])
        self.assertEqual(len(blocks), 1)
        self.assertEqual(len(blocks[0]["runs"]), 2)
        self.assertFalse(missing)

    def test_graph_context_recovers_unmapped_tail(self):
        graph = fake_graph(self.seq, [1, 3, 4], [1, 3, 5, 4])
        record = Alignment("A", "sample", "+", 0, 100, 0, 100, ((0, 100, 0),))
        blocks, missing, _ = plan_path(graph, "sample", [record])
        extra = [b for b in blocks if b["kind"] == "graph"]
        self.assertEqual(len(extra), 1)
        self.assertEqual(len(extra[0]["runs"]), 2)
        self.assertFalse(missing)

    def test_unbounded_tail_stays_missing(self):
        graph = fake_graph(self.seq, [1, 3], [1, 5])
        blocks, missing, _ = plan_path(graph, "sample", [])
        self.assertEqual(missing, [(100, 101)])

    def test_overlapping_physical_ownership_is_not_resolved_by_order(self):
        graph = fake_graph(self.seq, [1, 3], [1, 3])
        first = Alignment("A", "sample", "+", 0, 200, 0, 200, ((0, 200, 0),))
        second = Alignment("B", "sample", "+", 50, 150, 50, 150, ((50, 150, 50),))
        a = plan_path(graph, "sample", [first, second])
        b = plan_path(graph, "sample", [second, first])
        self.assertEqual(a[1], b[1])
        self.assertTrue(a[1])

    def test_repeated_and_reverse_steps_break_graph_chains(self):
        graph = fake_graph(self.seq, [1, 3, 4], [1, -3, 4])
        runs = graph_chains(
            graph.nodes,
            graph.refwalk,
            graph.ref_offsets,
            graph.walks["sample"],
            graph.spell("sample"),
            graph.reference,
        )
        self.assertTrue(all(len(chunk) == 1 for chunk in runs))

    def test_unknown_gap_breaks_graph_chain(self):
        seq = self.seq + ["NNNN"]
        graph = fake_graph(seq, [1, 3], [1, 6, 3])
        blocks, missing, _ = plan_path(graph, "sample", [])
        self.assertEqual(missing, [(100, 104)])
        self.assertTrue(all(len(b["runs"]) == 1 for b in blocks))

    def test_uniqueness_checks_overlapping_occurrences(self):
        self.assertFalse(unique_oriented("AAA", "AAAA", "AAA"))
        self.assertFalse(unique_oriented("ACGT", "ACGT", "ACGT"))

    def test_paf_record_identity_and_reverse_coordinates(self):
        with tempfile.TemporaryDirectory() as temp:
            p = Path(temp) / "a.paf"
            p.write_text(
                "s\t200\t20\t120\t-\tr\t500\t100\t201\t99\t101\t60\ttp:A:P\tcg:Z:50=1D50=\n"
            )
            r = read_alignments(p)["s"][0]
            self.assertEqual((r.query_start, r.query_end), (80, 180))
            self.assertEqual(r.runs, ((80, 130, 100), (130, 180, 151)))
            p.write_text(p.read_text().replace("50=1D50=", "100M"))
            with self.assertRaisesRegex(ValueError, "Exact"):
                read_alignments(p)

    def test_iupac_reverse_complement(self):
        s = "ACGTRYSWKMBDHVN"
        self.assertEqual(reverse_complement(reverse_complement(s)), s)

    def test_edit_replay_rejects_overlap(self):
        with self.assertRaises(ValueError):
            replay("ACGT", [Edit(1, 3, ""), Edit(2, 3, "T")])
        self.assertEqual(replay("ACGT", atomic_edits(0, "ACGT", "ATGT")), "ATGT")

    def test_origin_insertion_right_padding(self):
        self.assertEqual(vcf_alleles("ACGT", 0, 0, 0, ["G"], "EDIT"), (1, "A", ["GA"]))


if __name__ == "__main__":
    unittest.main()
