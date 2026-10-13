#!/usr/bin/env python3
import importlib.util,sys,tempfile,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'lib'))
spec=importlib.util.spec_from_file_location('unitary',ROOT/'lib/unitary_cross_species.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
class UnitaryTests(unittest.TestCase):
 def test_conservative_unitary_requires_synteny(self):
  self.assertEqual(m.classify_unitary(False,True,False,False),'unitary_like_needs_genome_scan')
  self.assertEqual(m.classify_unitary(False,True,True,True,False,False,False,True),'unitary_like_needs_outgroup')
  self.assertEqual(m.classify_unitary(False,True,True,True,True,True,False,True),'unitary_candidate')
  self.assertEqual(m.classify_unitary(True,True,True,True),'not_unitary_functional_copy_present')
 def test_proteome_hit_query_coverage(self):
  with tempfile.TemporaryDirectory() as td:
   p=Path(td)/'a.m8';p.write_text('Q1\tT1\t80\t80\t1\t80\t1\t80\t1e-20\t100\t100\t100\n')
   hit=m.read_best_proteome_hits(p)['Q1'];self.assertAlmostEqual(hit['query_coverage'],.8);self.assertAlmostEqual(hit['identity'],.8)
 def test_genome_copy_excludes_candidate_locus_but_detects_other_copy(self):
  row={'chr':'Chr1','start':'101','end':'200'};hits=[{'target':'Chr1','start':100,'end':200,'query_coverage':1,'identity':1,'intact':True}]
  self.assertIsNone(m.intact_genome_copy(hits,row));hits.append({'target':'Chr2','start':10,'end':100,'query_coverage':.9,'identity':.8,'intact':True});self.assertEqual(m.intact_genome_copy(hits,row)['target'],'Chr2')
if __name__=='__main__':unittest.main()
