"""GeneWise viewer counts must match the codon-aware analysis parser."""
import json
import re
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'lib'))
from visualize_genewise import generate_html_report, parse_genewise_block


def pretty(codon='TAA', aa='X', name='query', query='MAM'):
    pad = ' ' * 20
    return (f'Query protein: {name}\nScore 10.0 bits\n'
            + f'{name[:15]:15s} {1:3d} {query}\n'
            + pad + 'M M\n' + pad + 'M' + aa + 'M\n'
            + f'{"chr":15s} {100:3d} a{codon[0].lower()}a\n'
            + pad + f't{codon[1].lower()}t\n'
            + pad + f'g{codon[2].lower()}g\n//\n')


class VisualizationTests(unittest.TestCase):
    def test_verified_stops_are_counted_and_highlighted(self):
        for codon in ['TAA', 'TAG', 'TGA']:
            row = parse_genewise_block(pretty(codon))
            self.assertEqual(row['stop_count'], 1)
            self.assertEqual(row['target_seq'], 'M*M')
            self.assertEqual(row['mutations'][0]['type'], 'stop')
            self.assertEqual(row['mutations'][0]['position'], 2)

    def test_ambiguous_and_sense_x_are_not_stops(self):
        for codon in ['NNN', 'TGG']:
            row = parse_genewise_block(pretty(codon))
            self.assertEqual(row['stop_count'], 0)
            self.assertEqual(row['target_seq'], 'M?M')
            self.assertEqual(row['mutations'], [])

    def test_frameshift_matches_analysis_count(self):
        row = parse_genewise_block(pretty('GCT', '!'))
        self.assertEqual((row['stop_count'], row['frameshift_count']), (0, 1))
        self.assertEqual(row['mutations'][0]['type'], 'frameshift')

    def test_metadata_symbols_are_not_mutations(self):
        row = parse_genewise_block(pretty('GCT', 'A', 'query*!'))
        self.assertEqual((row['stop_count'], row['frameshift_count']), (0, 0))

    def test_truncated_label_preserves_alignment(self):
        row = parse_genewise_block(pretty(name='very_long_query_name_123456789'))
        self.assertEqual(row['query_seq'], 'MAM')
        self.assertEqual(row['length'], 3)
        self.assertAlmostEqual(row['identity'], 200 / 3)

    def test_html_json_is_safe_and_has_offline_fallback(self):
        row = parse_genewise_block(pretty())
        row['transcript_id'] = '</script><script>unsafe</script>'
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'viewer.html'
            generate_html_report([row], path, title='<demo>')
            text = path.read_text()
            self.assertIn('&lt;demo&gt;', text)
            self.assertNotIn(row['transcript_id'], text)
            data = json.loads(re.search(r'const allPseudogenes\s*=\s*(\[.*?\]);', text, re.S).group(1))
            self.assertEqual(data[0]['transcript_id'], row['transcript_id'])
            self.assertIn("typeof Chart === 'undefined'", text)
            self.assertIn('escapeHtml(item.transcript_id)', text)


if __name__ == '__main__':
    unittest.main()
