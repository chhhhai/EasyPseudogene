"""Regression tests for the September 2026 audit findings (no external tools required)."""
import csv
import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'lib'))
from genewise_alignment import mutation_counts,segments
from prefilter_proteins import parse_header
from prepare_genewise_jobs import prepare
from reference_cache import begin,artifacts
from types import SimpleNamespace

def module(name):
    spec=importlib.util.spec_from_file_location(name,ROOT/'scripts'/f'{name}.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m


def pretty(codon='TAA',aa='X',name='query',strand='forward'):
    pad=' '*20
    return f'Query protein: {name}\nStrand: {strand}\n'+f'{name[:15]:15s} {1:3d} '+'MAM\n'+pad+'M M\n'+pad+'M'+aa+'M\n'+f'{"chr":15s} {100:3d} '+'a'+codon[0].lower()+'a\n'+pad+'t'+codon[1].lower()+'t\n'+pad+'g'+codon[2].lower()+'g\n//\n'

class StopTests(unittest.TestCase):
    def test_genomic_label_runs_into_coordinate(self):
        prefix='chr15:10261975-57767 '
        pad=' '*len(prefix)
        block=(f'{"query":<16} {1:>3} MXM\n'+pad+'M M\n'+pad+'MXM\n'
               +prefix+'ata\n'+pad+'tga\n'+pad+'gag\n//\n')
        self.assertEqual(mutation_counts(block),(1,0))
        self.assertEqual(list(segments(block)),[('MXM','MXM')])
    def test_each_stop_codon(self):
        for stop in ['TAA','TAG','TGA']:self.assertEqual(mutation_counts(pretty(stop)),(1,0))
    def test_ambiguous_x_is_not_stop(self):
        self.assertEqual(mutation_counts(pretty('NNN')),(0,0));self.assertIn('?',list(segments(pretty('NNN')))[0][1])
    def test_sense_x_is_not_stop(self):self.assertEqual(mutation_counts(pretty('TGG')),(0,0))
    def test_query_metadata_markers_not_counted(self):self.assertEqual(mutation_counts(pretty('GCT',aa='A',name='q*!')),(0,0))
    def test_frameshift(self):self.assertEqual(mutation_counts(pretty('GCT',aa='!')),(0,1))
    def test_long_truncated_query_name(self):self.assertEqual(mutation_counts(pretty(name='long_query_name_abcdefghijk')),(1,0))
    def test_reverse_coding_rows(self):self.assertEqual(mutation_counts(pretty(strand='reverse')),(1,0))

class InputTests(unittest.TestCase):
    def test_plain_fasta(self):self.assertEqual(parse_header('>protein.1')[1],'sequence:protein.1')
    def test_gene_id_precedes_display_symbol(self):self.assertNotEqual(parse_header('>t1 gene=A gene_symbol=S')[1],parse_header('>t2 gene=B gene_symbol=S')[1])
    def test_isoforms_group_by_gene(self):self.assertEqual(parse_header('>t1 gene=A')[1],parse_header('>t2 gene=A')[1])
    def test_no_symbol_does_not_merge_plain_versions(self):self.assertNotEqual(parse_header('>p.1')[1],parse_header('>p.2')[1])

class FoldTests(unittest.TestCase):
    def test_label_and_model_folds_identical_across_species(self):
        a=module('label_candidates_from_matches');b=module('multispecies_benchmark')
        for family in ['F1','F2','orthogroup99']:
            self.assertEqual(a.stable_fold('mouse',family,5),b.stable_fold('rat',family,5))
    def test_loso_recomputes_inconsistent_input_folds(self):
        m=module('multispecies_benchmark')
        with tempfile.TemporaryDirectory() as d:
            d=Path(d);items=[]
            for species in ['a','b','c']:
                p=d/(species+'.tsv');p.write_text('candidate_id\tfamily_id\tfamily_fold\tx\ty\n'+''.join(f'{species}{i}\tF{i}\t{ord(species)%5}\t{i}\t{i%2}\n' for i in range(12)));items.append({'species':species,'labeled_tsv':str(p)})
            (d/'manifest.json').write_text(json.dumps({'datasets':items}));m.loso(d/'manifest.json',['x'],'y',d/'out.json',d/'pred')
            for p in (d/'pred').glob('*.tsv'):
                for r in csv.DictReader(p.read_text().splitlines(),delimiter='\t'):self.assertEqual(int(r['family_fold']),m.stable_fold('',r['family_id'],5))

class LocusTests(unittest.TestCase):
    def test_multi_locus_and_sanitization_collisions(self):
        with tempfile.TemporaryDirectory() as d:
            d=Path(d);(d/'proteins').write_text('>q/a\nMAM\n>q_a\nMAM\n')
            text='q/a\t+\tc1\t0\t9\nq/a\t-\tc2\t10\t19\nq_a\t+\tc1\t20\t29\n'
            (d/'hits').write_text(text);(d/'flanks').write_text(text);(d/'paf').write_text('q/a\t3\nq_a\t3\n')
            prepare(d/'proteins',d/'hits',d/'flanks',d,d/'jobs',d/'paf')
            records=list(csv.DictReader((d/'locus_map.tsv').read_text().splitlines(),delimiter='\t'))
            self.assertEqual(len({r['locus_id'] for r in records}),3);self.assertEqual(len(list(d.glob('*.pep.fa'))),3)
    def test_empty_jobs(self):
        with tempfile.TemporaryDirectory() as d:
            d=Path(d)
            for n in ['proteins','hits','flanks','paf']:(d/n).touch()
            prepare(d/'proteins',d/'hits',d/'flanks',d,d/'jobs',d/'paf');self.assertEqual((d/'jobs').read_text(),'')

class CoordinateTests(unittest.TestCase):
    def test_one_base_interval_and_both_strands(self):
        m=module('normalize_pgre_coordinates')
        with tempfile.TemporaryDirectory() as d:
            d=Path(d);(d/'input').write_text('candidate_id\tchr\tstart\tend\tstrand\nx\tc\t1\t1\t+\ny\tc\t5\t8\t-\n');m.convert(d/'input',d/'out')
            self.assertEqual((d/'out').read_text(),'c\t0\t1\tx\t0\t+\nc\t4\t8\ty\t0\t-\n')
    def test_refuses_zero_based_input(self):
        m=module('normalize_pgre_coordinates')
        with tempfile.TemporaryDirectory() as d:
            d=Path(d);(d/'i').write_text('chr\tstart\tend\nc\t0\t9\n')
            with self.assertRaises(ValueError):m.convert(d/'i',d/'o')

class CacheTests(unittest.TestCase):
    def test_stale_incomplete_and_modified_outputs_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            d=Path(d);source=d/'source';out=d/'out';source.mkdir();out.mkdir();inp=d/'input';inp.write_text('first')
            a=SimpleNamespace(source=source,directory=out,input=[inp],arguments=['run'])
            begin(a)
            with self.assertRaisesRegex(ValueError,'incomplete'):begin(a)
            p=out/'reference_cache.json';j=json.loads(p.read_text());j.update(status='complete',artifacts={});p.write_text(json.dumps(j));begin(a)
            j=json.loads(p.read_text());j.update(status='complete',artifacts={});p.write_text(json.dumps(j));inp.write_text('second')
            with self.assertRaisesRegex(ValueError,'changed'):begin(a)
            inp.write_text('first');(out/'regions').mkdir();(out/'regions/x').write_text('tampered')
            with self.assertRaisesRegex(ValueError,'missing or modified'):begin(a)

if __name__=='__main__':unittest.main()
