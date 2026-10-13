#!/usr/bin/env python3
import importlib.util
import tempfile
import threading
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class CatalogUtilsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.utils = load_module("catalog_utils", ROOT / "lib" / "catalog_utils.py")

    def test_file_sha256_is_stable(self):
        import hashlib
        with tempfile.TemporaryDirectory() as td:
            path=Path(td)/"x";path.write_bytes(b"abc")
            self.assertEqual(self.utils.file_sha256(path),hashlib.sha256(b"abc").hexdigest())

    def test_merge_intervals_uses_half_open_union(self):
        observed = self.utils.merge_half_open_intervals([(10, 20), (15, 25), (30, 40), (40, 45)])
        self.assertEqual(observed, [(10, 25), (30, 45)])

    def test_interval_index_does_not_double_count_overlapping_cds(self):
        idx = self.utils.IntervalIndex({"Chr1": [(100, 200), (150, 250)]})
        self.assertEqual(idx.overlap_bp("Chr1", 100, 250), 150)
        self.assertAlmostEqual(idx.overlap_fraction("Chr1", 100, 300), 0.75)

    def test_union_coverage_excludes_internal_gaps(self):
        cov = self.utils.union_coverage_fraction([(1, 30), (51, 80), (70, 90)], query_length=100, one_based_inclusive=True)
        self.assertAlmostEqual(cov, 0.70)

    def test_generic_gene_id_parsing_prefers_gene_field(self):
        desc = "AT5G16970.1 pep chromosome:TAIR10:5:1:100:+ gene:AT5G16970 transcript:AT5G16970.1 gene_symbol:AER"
        self.assertEqual(self.utils.infer_gene_id("AT5G16970.1", desc), "AT5G16970")

    def test_generic_gene_id_parsing_supports_equals_header(self):
        desc = "XP_001 gene=GENE42 transcript_id=TX42"
        self.assertEqual(self.utils.infer_gene_id("XP_001", desc), "GENE42")

    def test_reverse_one_based_inclusive_coverage(self):
        cov = self.utils.union_coverage_fraction([(30, 1)], query_length=100, one_based_inclusive=True)
        self.assertAlmostEqual(cov, 0.30)

    def test_gene_id_parsing_supports_brackets_and_uniprot_gn(self):
        self.assertEqual(self.utils.infer_gene_id("XP1", "protein [gene=ABC1] [locus_tag=L1]"), "ABC1")
        self.assertEqual(self.utils.infer_gene_id("P1", "sp|P1|NAME GN=GENE9 PE=1"), "GENE9")

    def test_atomic_writer_uses_unique_temps_under_concurrency(self):
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / "shared.txt"
            errors = []
            def worker(text):
                try:
                    with self.utils.atomic_text_writer(path) as out:
                        out.write(text)
                except Exception as exc:
                    errors.append(exc)
            threads = [threading.Thread(target=worker, args=(f"value-{i}",)) for i in range(8)]
            for thread in threads: thread.start()
            for thread in threads: thread.join()
            self.assertFalse(errors)
            self.assertTrue(path.read_text().startswith("value-"))
            self.assertEqual(list(Path(td).glob("shared.txt.tmp.*")), [])

    def test_atomic_writer_replaces_destination(self):
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / "out.txt"
            with self.utils.atomic_text_writer(path) as out:
                out.write("complete\n")
            self.assertEqual(path.read_text(), "complete\n")
            self.assertFalse((Path(td) / "out.txt.tmp").exists())


if __name__ == "__main__":
    unittest.main(verbosity=2)
