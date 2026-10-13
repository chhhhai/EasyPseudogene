"""Regression checks for input exclusions and the shared result schema."""
import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'lib'))
from normalize_inputs import normalize_proteins, plain_genome
from publish_catalog_results import publish
from result_evidence import evidence, grade
from filter_self_mode_candidates import main as filter_self
from prepare_genewise_jobs import prepare
from reference_cache import begin


class EvidenceTests(unittest.TestCase):
    def test_thresholds(self):
        self.assertEqual(grade(1, 0, .8, .5, True)[0], 'high')
        self.assertEqual(grade(0, 1, .2, .3, True)[0], 'medium')
        self.assertEqual(grade(1, 0, .199, .3, True)[0], 'low')
        self.assertEqual(grade(0, 0, 1, 1, True)[0], 'low')

    def test_missing_and_invalid_measurements(self):
        for coverage, identity in [(None, .5), (.5, None), (1.1, .5), (.5, -1)]:
            self.assertEqual(grade(1, 0, coverage, identity, True)[0], 'unassessed')
        self.assertEqual(evidence(1, None, 1, 1, assessed=True)['evidence_level'], 'unassessed')


class NormalizationTests(unittest.TestCase):
    def test_invalid_protein_skip_report(self):
        with tempfile.TemporaryDirectory() as td:
            d = Path(td)
            (d / 'source').write_text('>ok description\nMAM*\n>bad\nMA*M\n')
            result = normalize_proteins(d / 'source', d / 'output', None, {}, d / 'mapping', 'skip', d / 'exclusions')
            self.assertEqual(result['proteins'], 1)
            self.assertEqual(result['excluded_proteins'], 1)
            self.assertIn('>ok description\nMAM\n', (d / 'output').read_text())
            excluded = list(csv.DictReader((d / 'exclusions').open(), delimiter='\t'))
            self.assertEqual(excluded[0]['query_id'], 'bad')
            self.assertEqual(excluded[0]['reason'], 'internal_stop')

    def test_duplicate_invalid_ids_are_not_silently_skipped(self):
        with tempfile.TemporaryDirectory() as td:
            d = Path(td)
            (d / 'source').write_text('>bad\nM*M\n>bad\nMAM\n')
            with self.assertRaisesRegex(ValueError, 'Duplicate protein FASTA ID'):
                normalize_proteins(d / 'source', d / 'output', None, {}, d / 'mapping', 'skip', d / 'exclusions')

    def test_plain_genome_is_local_and_repeated_staging_keeps_index(self):
        with tempfile.TemporaryDirectory() as td:
            d = Path(td)
            source, output = d / 'source.fa', d / 'stage.fa'
            source.write_text('>c\nATGC\n')
            self.assertEqual(plain_genome(source, output), output.resolve())
            self.assertNotEqual(source.resolve(), output.resolve())
            index = Path(str(output) + '.fai')
            index.write_text('existing index\n')
            plain_genome(source, output)
            self.assertTrue(index.exists())
            source.write_text('>c\nATGCAA\n')
            plain_genome(source, output)
            self.assertFalse(index.exists())


class CatalogPublisherTests(unittest.TestCase):
    def setup_files(self, d, structural='1', start='105', end='350', parent='q1'):
        proteins = d / 'proteins.fa'
        proteins.write_text('>q1\n' + 'M' * 100 + '\n')
        raw = dict(candidate_id='cand1', query_id='q1', chr='c', start='101', end='400', strand='+')
        validated = dict(raw, validator_status='validated_miniprot', alignment_source='miniprot',
                         structural_tags_complete=structural, selected_parent=parent,
                         validated_genomic_start=start, validated_genomic_end=end,
                         validated_strand='-', validated_stop_count='1', validated_frameshift_count='0',
                         validated_query_coverage='.9', validated_identity='.8', validated_aligned_length='90')
        paths = []
        for name, row in [('raw.tsv', raw), ('validated.tsv', validated)]:
            p = d / name
            with p.open('w', newline='') as h:
                w = csv.DictWriter(h, fieldnames=list(row), delimiter='\t')
                w.writeheader()
                w.writerow(row)
            paths.append(p)
        return proteins, paths[0], paths[1]

    def test_validated_disruption_uses_refined_coordinates(self):
        with tempfile.TemporaryDirectory() as td:
            d = Path(td)
            p, raw, validated = self.setup_files(d)
            result = publish(d / 'run', 'self', 'self', 'mmseqs-fast', p, raw, validated)
            self.assertEqual(result['final_loci'], 1)
            final = list(csv.DictReader((d / 'run/results/final_loci.tsv').open(), delimiter='\t'))[0]
            self.assertEqual((final['start0'], final['end0'], final['strand']), ('104', '350', '-'))
            self.assertEqual(final['evidence_level'], 'high')
            self.assertEqual(final['genewise_status'], 'not_run')
            self.assertEqual(final['evidence_query_id'], 'q1')

    def test_incomplete_structural_tags_are_unassessed(self):
        with tempfile.TemporaryDirectory() as td:
            d = Path(td)
            p, raw, validated = self.setup_files(d, structural='0')
            result = publish(d / 'run', 'self', 'self', 'mmseqs-fast', p, raw, validated)
            self.assertEqual(result['final_loci'], 0)
            self.assertEqual(result['evidence_levels'], {'unassessed': 1})
            locus = list(csv.DictReader((d / 'run/results/loci.tsv').open(), delimiter='\t'))[0]
            self.assertEqual(locus['stop_count'], '')

    def test_nonoverlapping_validation_is_unassessed(self):
        with tempfile.TemporaryDirectory() as td:
            d = Path(td)
            p, raw, validated = self.setup_files(d, start='500', end='700')
            result = publish(d / 'run', 'self', 'self', 'mmseqs-fast', p, raw, validated)
            self.assertEqual(result['evidence_levels'], {'unassessed': 1})

    def test_unknown_validated_parent_is_rejected(self):
        with tempfile.TemporaryDirectory() as td:
            d = Path(td)
            p, raw, validated = self.setup_files(d, parent='missing')
            with self.assertRaisesRegex(ValueError, 'Validated parent absent'):
                publish(d / 'run', 'self', 'self', 'mmseqs-fast', p, raw, validated)

    def test_unvalidated_catalog_is_not_final(self):
        with tempfile.TemporaryDirectory() as td:
            d = Path(td)
            p, raw, validated = self.setup_files(d)
            result = publish(d / 'run', 'self', 'self', 'mmseqs-fast', p, raw)
            self.assertEqual(result['final_loci'], 0)
            self.assertEqual(result['evidence_levels'], {'unassessed': 1})

    def test_empty_catalog_publishes_layout(self):
        with tempfile.TemporaryDirectory() as td:
            d = Path(td)
            p, raw, validated = self.setup_files(d)
            for path in [raw, validated]:
                path.write_text(path.read_text().splitlines()[0] + '\n')
            result = publish(d / 'run', 'self', 'self', 'mmseqs-fast', p, raw, validated)
            self.assertEqual(result['selected_loci'], 0)
            self.assertEqual(result['input_queries'], 1)
            self.assertTrue((d / 'run/results/final_loci.bed').is_file())
            self.assertEqual(json.loads((d / 'run/metadata/run_manifest.json').read_text())['layout_version'], 3)


class SelfLocusTests(unittest.TestCase):
    def test_multiple_query_loci_keep_their_own_flanks(self):
        with tempfile.TemporaryDirectory() as td:
            d = Path(td)
            (d / 'hits').write_text('q\t+\tc1\t100\t130\t1000\t30\t0\nq\t+\tc1\t200\t230\t1000\t30\t0\nq\t-\tc2\t300\t330\t1000\t30\t0\n')
            (d / 'flanks').write_text('q\t+\tc1\t90\t140\nq\t+\tc1\t190\t240\nq\t-\tc2\t290\t340\n')
            (d / 'gff').write_text('##gff-version 3\nc1\tx\tgene\t101\t130\t.\t+\t.\tID=own\n')
            argv = ['filter']
            for flag, name in [('best-hits', 'hits'), ('best-hits-flanked', 'flanks'), ('gff', 'gff'),
                               ('out-best-hits', 'kept'), ('out-best-hits-flanked', 'kept_flanks'),
                               ('report', 'report'), ('bed', 'bed'), ('removed', 'removed')]:
                argv += ['--' + flag, str(d / name)]
            with patch.object(sys, 'argv', argv):
                filter_self()
            self.assertEqual((d / 'kept_flanks').read_text(), 'q\t+\tc1\t190\t240\nq\t-\tc2\t290\t340\n')
            (d / 'proteins').write_text('>q\nMAM\n')
            (d / 'paf').write_text('q\t3\n')
            prepare(d / 'proteins', d / 'kept', d / 'kept_flanks', d, d / 'jobs', d / 'paf')
            self.assertEqual(len((d / 'jobs').read_text().splitlines()), 2)
            self.assertEqual((d / 'jobs').read_text().splitlines()[0].split('\t')[3:5], ['191', '240'])

    def test_same_chromosome_wrong_flank_is_rejected(self):
        with tempfile.TemporaryDirectory() as td:
            d = Path(td)
            (d / 'proteins').write_text('>q\nMAM\n')
            (d / 'hits').write_text('q\t+\tc\t100\t130\n')
            (d / 'flanks').write_text('q\t+\tc\t200\t240\n')
            (d / 'paf').touch()
            with self.assertRaisesRegex(ValueError, 'does not contain'):
                prepare(d / 'proteins', d / 'hits', d / 'flanks', d, d / 'jobs', d / 'paf')


class CompletedCacheTests(unittest.TestCase):
    def test_reuse_retains_completed_cache(self):
        with tempfile.TemporaryDirectory() as td:
            d = Path(td)
            out, source = d / 'out', d / 'source'
            out.mkdir()
            source.mkdir()
            for name in ['lib', 'scripts', 'bin', 'config']:
                (source / name).mkdir()
            a = SimpleNamespace(directory=out, source=source, input=[], arguments=[], reuse_complete=True)
            self.assertFalse(begin(a))
            path = out / 'reference_cache.json'
            data = json.loads(path.read_text())
            data.update(status='complete', artifacts={})
            path.write_text(json.dumps(data))
            original = path.read_bytes()
            self.assertTrue(begin(a))
            self.assertEqual(path.read_bytes(), original)
            a.arguments = ['changed']
            with self.assertRaisesRegex(ValueError, 'changed'):
                begin(a)


if __name__ == '__main__':
    unittest.main()
