"""Regression checks for output paths, GFF3 exports and offline report safety."""
import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from urllib.parse import unquote

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'lib'))
from publish_result_assets import gff3, html_table, publish_assets
from reference_cache import artifacts, begin
from run_layout import cache_file, work_directory


def locus(**overrides):
    row = dict(locus_id='cand1', query_id='q1', chrom='chr1', start0='0', end0='9', strand='+',
               evidence_level='high', evidence_method='miniprot', stop_count='1', frameshift_count='0')
    row.update(overrides)
    return row


class GffTests(unittest.TestCase):
    def test_one_based_conversion_and_reverse_strand(self):
        fields = gff3([locus(start0='9', end0='20', strand='-')]).splitlines()[-1].split('\t')
        self.assertEqual(fields[3:5], ['10', '20'])
        self.assertEqual(fields[6], '-')
        self.assertEqual(fields[2], 'pseudogene')

    def test_reserved_characters_are_escaped(self):
        row = locus(locus_id='a;b=c% d', query_id='q;= %', chrom='c;=%')
        fields = gff3([row]).splitlines()[-1].split('\t')
        attrs = {k: unquote(v) for k, v in (x.split('=', 1) for x in fields[8].split(';'))}
        self.assertEqual(attrs['ID'], row['locus_id'])
        self.assertEqual(attrs['reference_protein'], row['query_id'])
        self.assertEqual(unquote(fields[0]), row['chrom'])

    def test_empty_gff_is_valid_header(self):
        result = gff3([])
        self.assertTrue(result.startswith('##gff-version 3\n'))
        self.assertFalse(any(line and not line.startswith('#') for line in result.splitlines()))

    def test_duplicate_ids_rejected(self):
        with self.assertRaisesRegex(ValueError, 'duplicate'):
            gff3([locus(), locus()])

    def test_invalid_coordinates_and_strand_rejected(self):
        for overrides in [dict(start0='-1'), dict(end0='0'), dict(strand='x')]:
            with self.subTest(overrides=overrides), self.assertRaisesRegex(ValueError, 'Invalid'):
                gff3([locus(**overrides)])


class LayoutTests(unittest.TestCase):
    def test_current_and_legacy_cache_locations(self):
        with tempfile.TemporaryDirectory() as td:
            d = Path(td)
            self.assertEqual(work_directory(d), d)
            self.assertEqual(cache_file(d), d / 'reference_cache.json')
            (d / 'work').mkdir()
            self.assertEqual(work_directory(d), d / 'work')
            self.assertEqual(cache_file(d), d / 'metadata/reference_cache.json')

    def test_preflight_does_not_hide_legacy_stage_files(self):
        with tempfile.TemporaryDirectory() as td:
            d = Path(td)
            (d / 'work/input_check').mkdir(parents=True)
            (d / 'genewise').mkdir()
            self.assertEqual(work_directory(d), d)
            self.assertEqual(cache_file(d), d / 'reference_cache.json')

    def test_cached_report_changes_are_rejected(self):
        with tempfile.TemporaryDirectory() as td:
            d = Path(td)
            run, source = d / 'run', d / 'source'
            (run / 'work').mkdir(parents=True)
            source.mkdir()
            args = SimpleNamespace(directory=run, source=source, input=[], arguments=[], reuse_complete=True)
            self.assertFalse(begin(args))
            (run / 'report.html').write_text('original')
            cache = cache_file(run)
            data = json.loads(cache.read_text())
            data.update(status='complete', artifacts=artifacts(run))
            self.assertNotIn('metadata/reference_cache.json', data['artifacts'])
            self.assertIn('report.html', data['artifacts'])
            cache.write_text(json.dumps(data))
            original_cache = cache.read_bytes()
            self.assertTrue(begin(args))
            self.assertEqual(cache.read_bytes(), original_cache)
            (run / 'report.html').write_text('changed')
            with self.assertRaisesRegex(ValueError, 'missing or modified'):
                begin(args)
            self.assertEqual(cache.read_bytes(), original_cache)


class ReportTests(unittest.TestCase):
    def test_table_escapes_html(self):
        rendered = html_table([('id', '<name>')], [dict(id='<script>alert(1)</script>')])
        self.assertNotIn('<script>', rendered)
        self.assertIn('&lt;script&gt;', rendered)
        self.assertIn('&lt;name&gt;', rendered)

    def populate(self, d, rows):
        (d / 'work').mkdir()
        (d / 'results').mkdir()
        (d / 'metadata').mkdir()
        with (d / 'results/final_loci.tsv').open('w', newline='') as h:
            w = csv.DictWriter(h, fieldnames=list(locus()), delimiter='\t')
            w.writeheader()
            w.writerows(rows)
        (d / 'results/summary.json').write_text(json.dumps(dict(final_loci=len(rows), input_queries=1,
                                                selected_loci=len(rows), evidence_levels={'high': len(rows)})))
        (d / 'metadata/run_manifest.json').write_text(json.dumps(dict(species='<script>sample</script>',
                                                      comparison_mode='self', workflow='fast', files={})))

    def test_report_and_gff_written_for_current_layout(self):
        with tempfile.TemporaryDirectory() as td:
            d = Path(td)
            self.populate(d, [locus()])
            publish_assets(d)
            manifest = json.loads((d / 'metadata/run_manifest.json').read_text())
            self.assertEqual(manifest['layout_version'], 4)
            self.assertEqual(manifest['files']['final_gff3'], '../results/final_loci.gff3')
            self.assertEqual(manifest['files']['reference_cache'], 'reference_cache.json')
            report = (d / 'report.html').read_text()
            self.assertIn('results/final_loci.gff3', report)
            self.assertNotIn('<script>sample</script>', report)
            self.assertIn('&lt;script&gt;sample&lt;/script&gt;', report)

    def test_zero_candidates_have_report_and_header_only_gff(self):
        with tempfile.TemporaryDirectory() as td:
            d = Path(td)
            self.populate(d, [])
            publish_assets(d)
            self.assertTrue((d / 'report.html').is_file())
            self.assertEqual((d / 'results/final_loci.gff3').read_text(), gff3([]))

    def test_inconsistent_counts_fail_publication(self):
        with tempfile.TemporaryDirectory() as td:
            d = Path(td)
            self.populate(d, [])
            (d / 'results/summary.json').write_text(json.dumps(dict(final_loci=1)))
            with self.assertRaisesRegex(ValueError, 'counts differ'):
                publish_assets(d)
            self.assertFalse((d / 'report.html').exists())


if __name__ == '__main__':
    unittest.main()
