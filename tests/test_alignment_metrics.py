"""Regression checks for GeneWise intron and gap columns."""
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))
from alignment_metrics import (aligned_query_length, flank_identities,
                               mutation_query_fraction, query_coverage,
                               sequence_identity)


class AlignmentMetricsTests(unittest.TestCase):
    def test_spaces_and_gaps_do_not_inflate_coverage(self):
        query = "AB    CD--E"
        target = "AB    CX--E"
        self.assertEqual(aligned_query_length(query, target), 5)
        self.assertEqual(query_coverage(query, target, 5), 1.0)
        self.assertEqual(sequence_identity(query, target), .8)

    def test_query_coverage_cannot_exceed_one(self):
        self.assertEqual(query_coverage("AAAAA", "AAAAA", 4), 1.0)

    def test_flanks_use_amino_acid_positions(self):
        query = "AB    CDE"
        target = "AB    XDF"
        self.assertEqual(flank_identities(query, target, 6, 2), (1.0, .5))

    def test_mutation_fraction_ignores_intron_columns(self):
        self.assertEqual(mutation_query_fraction("AB    CDE", "AB    XDE", 6), (.6, 3, 5))


if __name__ == "__main__":
    unittest.main()
