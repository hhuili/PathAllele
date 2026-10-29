"""Independent boundary and multiallelic examples for public VCF classification."""

import unittest

from pathallele.variant_schema import classify_allele, classification_info


class VariantSchemaTests(unittest.TestCase):
    def check(self, ref, alt, expected):
        row = classify_allele(ref, alt)
        self.assertEqual((row.type, row.sizeclass, row.vclass, row.reflen, row.altlen), expected)

    def test_single_base_substitution_with_padding(self):
        self.check("AAC", "AGC", ("SNV", "SMALL", "SNV", 1, 1))

    def test_mnv_and_equal_length_large_variant(self):
        self.check("GACG", "GGTG", ("MNV", "SMALL", "MNV", 2, 2))
        self.check("G" + "A" * 50 + "G", "G" + "C" * 50 + "G", ("MNV", "SV", "SV", 50, 50))

    def test_indel_size_boundary_excludes_padding(self):
        for size, group in [(49, "SMALL"), (50, "SV")]:
            cat = "INDEL" if group == "SMALL" else "SV"
            self.check("G", "G" + "A" * size, ("INS", group, cat, 0, size))
            self.check("G" + "A" * size, "G", ("DEL", group, cat, size, 0))

    def test_replacement_size_is_not_absolute_svlen(self):
        self.check("G" + "A" * 50, "G" + "C" * 49, ("REPLACEMENT", "SV", "SV", 50, 49))
        self.check("GAA", "GCCC", ("REPLACEMENT", "SMALL", "INDEL", 2, 3))

    def test_mixed_multiallelic_tags_follow_alt_order(self):
        tags = classification_info("AC", ["AT", "A", "AC" + "G" * 50])
        self.assertEqual(tags["TYPE"], ("SNV", "DEL", "INS"))
        self.assertEqual(tags["SIZECLASS"], ("SMALL", "SMALL", "SV"))
        self.assertEqual(set(tags), {"TYPE", "SIZECLASS"})

    def test_reject_non_sequence_and_reference_alleles(self):
        for ref, alt in [("A", "A"), ("A", "<DEL>"), ("N", "A"), ("A", "*"), ("A", "")]:
            with self.subTest(ref=ref, alt=alt), self.assertRaises(ValueError):
                classify_allele(ref, alt)


if __name__ == "__main__":
    unittest.main()
