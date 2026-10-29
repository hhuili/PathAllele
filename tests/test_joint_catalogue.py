import gzip
import random
import tempfile
from pathlib import Path
import unittest
import numpy as np

from pathallele.allele_model import Edit, replay
from pathallele.graph_io import reverse_complement
from pathallele.joint_catalogue import JointCatalogue, factor_state
from pathallele.sequence_anchors import exact_runs, obligatory_chain, refine_anchors


class JointTests(unittest.TestCase):
    def test_reference_backbone_requires_emitted_sequence_witness(self):
        ref = "ACGTCAGT"
        es = [Edit(2, 3, "A"), Edit(5, 6, "T")]
        sequence = replay(ref, es)
        self.assertEqual(
            factor_state(ref, 0, 8, sequence, es, set()), (None, "UNRESOLVED_COMPONENT")
        )
        self.assertEqual(
            factor_state(ref, 0, 8, sequence, es, {es[0]}), (None, "UNRESOLVED_COMPONENT")
        )
        self.assertEqual(factor_state(ref, 0, 8, sequence, es, set(es)), (ref, "PASS"))
        full = Edit(0, 8, sequence, True)
        self.assertEqual(factor_state(ref, 0, 8, sequence, [full], set()), (sequence, "PASS"))
        with self.assertRaises(ValueError):
            factor_state(ref, 0, 8, None, es, set(es))
        with self.assertRaises(ValueError):
            factor_state(ref, 0, 8, "TTTTTTTT", es, set(es))

    def test_complete_catalogue_keeps_reference_like_haplotypes(self):
        with tempfile.TemporaryDirectory() as tmp:
            c = JointCatalogue(["reference", "small", "compound", "missing"], tmp)
            observations = [("ACGT", "PASS"), ("ATGT", "PASS"), ("TGCA", "PASS"), (None, "NO_DATA")]
            projected = [("ACGT", "PASS"), ("ACGT", "PASS"), ("TGCA", "PASS"), (None, "NO_DATA")]
            c.add("v1", 4, 8, "ACGT", observations, projected)
            c.close()
            data = np.load(Path(tmp) / "Graph.locus_genotypes.npz")
            self.assertEqual(data["gt"].tolist(), [[0, 1, 2, -1]])
            with gzip.open(Path(tmp) / "Graph.loci.tsv.gz", "rt") as f:
                self.assertIn("ATGT", f.read())

    def test_obligatory_chain_preserves_ambiguous_placement(self):
        self.assertEqual(obligatory_chain([(0, 0, 5), (10, 0, 5)]), [])
        self.assertEqual(obligatory_chain([(0, 0, 5), (10, 0, 5), (20, 20, 10)]), [(20, 20, 10)])
        self.assertEqual(obligatory_chain([(0, 0, 100), (2, 2, 20)]), [(0, 0, 100)])

    def test_sequence_boundaries_do_not_need_edit_alignment(self):
        rng = random.Random(1307)
        middle = "".join(rng.choice("ACGT") for _ in range(96))
        ref = "A" * 60 + "CG" + middle + "C"
        alt = "CG" + "T" * 59 + "C" + middle + "A"
        for k in (19, 31, 51):
            runs = exact_runs(ref, alt, k)
            self.assertIn((62, 62, 96), runs)
            reverse = exact_runs(reverse_complement(ref), reverse_complement(alt), k)
            back = sorted((len(ref) - r - n, len(alt) - q - n, n) for r, q, n in reverse)
            self.assertEqual(runs, back)

    def test_no_crossing_of_physical_anchor_or_unknown_sequence(self):
        rng = random.Random(21)
        seq = "".join(rng.choice("ACGT") for _ in range(70))
        self.assertEqual(exact_runs(seq, seq[:35] + "N" + seq[36:]), [])
        ref = "G" + seq + "A" * 20 + "T"
        alt = "G" + seq + "C" * 20 + "T"
        anchors = [(1, 1, len(seq), 7)]
        self.assertEqual(refine_anchors(ref, alt, anchors), anchors)


class CompleteUnitTests(unittest.TestCase):
    def test_same_footprint_snv_is_not_emitted_twice(self):
        ref = "ACGTA"
        edit = Edit(2, 3, "A")
        self.assertEqual(factor_state(ref, 2, 3, "A", [edit], set()), ("A", "PASS"))
        self.assertEqual(factor_state(ref, 2, 3, "A", [edit], {edit}), ("G", "PASS"))

    def test_rejected_normalization_does_not_mask_stable_distant_snv(self):
        from pathallele.normalize_alleles import directional_disagreements

        ref = "A" * 80
        forward = [Edit(1, 2, "G"), Edit(40, 50, "")]
        reverse = [Edit(1, 2, "G"), Edit(41, 51, "")]
        uncertain = directional_disagreements(ref, 0, forward, reverse)
        self.assertTrue(uncertain)
        self.assertNotIn(1, uncertain)
        self.assertTrue(all(p >= 39 for p in uncertain))

    def test_indel_normalization_respects_physical_parent_boundary(self):
        from pathallele.normalize_alleles import bounded_indel_normalization

        ref = "AAAAA"
        self.assertEqual(bounded_indel_normalization(ref, [Edit(4, 4, "A")]), [Edit(1, 1, "A")])
        self.assertEqual(bounded_indel_normalization(ref, [Edit(3, 4, "")]), [Edit(1, 2, "")])

    def test_joint_model_replay_catches_collapsed_repeat_copy(self):
        from pathallele.validate_models import validate_models
        import numpy as np

        ref = "GAAAAAC"
        parents = [
            dict(pos=1, alleles=["GAA", "GAAA"], carriers=[[], [0]]),
            dict(pos=4, alleles=["AAA", "AAAA"], carriers=[[], [0]]),
        ]
        known = np.ones((1, len(ref)), dtype=np.uint8)
        ambiguous = np.zeros_like(known, dtype=bool)
        valid = [{Edit(1, 1, "A"), Edit(4, 4, "A")}]
        self.assertGreater(
            validate_models(parents, valid, ref, known, ambiguous)["joint_source_windows_replayed"],
            0,
        )
        with self.assertRaisesRegex(ValueError, "physical sequence"):
            validate_models(parents, [{Edit(1, 1, "A")}], ref, known, ambiguous)

    def test_touching_normalized_units_are_one_sequence(self):
        from pathallele.normalize_alleles import close_touching_units

        ref = "AACCGGTT"
        edits = [Edit(2, 4, "G", True), Edit(4, 6, "")]
        closed = close_touching_units(ref, edits)
        self.assertEqual(len(closed), 1)
        self.assertEqual(replay(ref, closed), replay(ref, edits))

    def test_joint_coherence_does_not_absorb_unique_bases_into_ambiguous_span(self):
        from pathallele.normalize_alleles import cohere_complex_models
        import numpy as np

        complete = Edit(1, 4, "TGG", True)
        primitive = {Edit(1, 2, "T"), Edit(2, 3, "G"), Edit(3, 4, "G")}
        models = [set(primitive), set(primitive)]
        known = np.ones((2, 5), dtype=np.uint8)
        known[0, 3] = 2
        cohere_complex_models(models, "AAAAA", known, [complete])
        self.assertEqual(models[0], primitive)
        self.assertEqual(models[1], {complete})

    def test_mixed_units_do_not_depend_on_alignment_gaps(self):
        from pathallele.project_variants import decompose_alleles

        parent = dict(
            pid="synthetic",
            alleles=["GACGTAC", "GTGGTTC", "GACGTTAC"],
            anchors=[[], [], []],
            native_kind=["REFERENCE", "SEQUENCE", "SEQUENCE"],
        )
        result = decompose_alleles([parent])[0]
        self.assertEqual(len(result[1]), 1)
        self.assertTrue(result[1][0].complex)
        self.assertEqual(result[1][0].end - result[1][0].start, len(result[1][0].alt))
        self.assertEqual(replay(parent["alleles"][0], result[2]), parent["alleles"][2])


if __name__ == "__main__":
    unittest.main()
