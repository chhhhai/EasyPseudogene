#!/usr/bin/env python3
import importlib.util,sys,tempfile,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]; sys.path.insert(0,str(ROOT/'lib'))
def load():
 spec=importlib.util.spec_from_file_location('validator',ROOT/'lib'/'validate_candidate_disruptions.py'); m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m
class ValidatorTests(unittest.TestCase):
 @classmethod
 def setUpClass(cls): cls.m=load()
 def test_expected_query_target_pair_and_tags(self):
  with tempfile.TemporaryDirectory() as td:
   paf=Path(td)/'a.paf'; paf.write_text('Q1\t100\t10\t90\t+\tEPWIN_C1\t500\t20\t300\t70\t100\t60\tAS:i:120\tfs:i:2\tst:i:1\nQ2\t100\t0\t100\t+\tEPWIN_C1\t500\t20\t300\t90\t100\t60\tAS:i:200\tfs:i:0\tst:i:0\n')
   hits,s=self.m.read_best_family_paf(paf,{'EPWIN_C1':{'candidate_id':'C1','initial_query':'Q1','expected_queries':['Q1'],'sequence':'A'*500}})
   self.assertEqual(s['validated_candidates'],1); self.assertEqual(hits['C1']['frameshift_count'],2); self.assertAlmostEqual(hits['C1']['query_coverage'],.8)
 def test_alternative_parent_can_validate_candidate(self):
  with tempfile.TemporaryDirectory() as td:
   paf=Path(td)/'a.paf'; paf.write_text('QALT\t100\t0\t90\t+\tEPWIN_C1\t500\t20\t300\t80\t100\t60\tAS:i:180\tfs:i:1\tst:i:0\tcg:Z:40M100N50M\n')
   mapping={'EPWIN_C1':{'candidate_id':'C1','initial_query':'Q1','expected_queries':['Q1','QALT'],'sequence':'A'*500}}
   hits,_summary=self.m.read_best_family_paf(paf,mapping)
   self.assertEqual(hits['C1']['query_id'],'QALT'); self.assertTrue(hits['C1']['parent_changed']); self.assertEqual(hits['C1']['candidate_intron_count'],1)

 def test_subtype_rules(self):
  processed={'query_coverage':.8,'frameshift_count':1,'stop_count':0,'candidate_intron_count':0}
  duplicated={'query_coverage':.8,'frameshift_count':1,'stop_count':0,'candidate_intron_count':2}
  self.assertEqual(self.m.classify_subtype(processed,3),'processed_like')
  self.assertEqual(self.m.classify_subtype(duplicated,3),'duplicated_like')

 def test_target_site_duplication(self):
  length,sequence=self.m.target_site_duplication('GGATGCAAAAATGCTT',6,10,min_len=4,max_len=6)
  self.assertEqual((length,sequence),(4,'ATGC'))
 def test_repeat_classification_and_category_fractions(self):
  self.assertEqual(self.m.classify_repeat_label('LTR/Gypsy'),'retroelement');self.assertEqual(self.m.classify_repeat_label('DNA/Helitron'),'dna_transposon');self.assertEqual(self.m.classify_repeat_label('Low_complexity'),'low_complexity')
  index=self.m.RepeatEvidenceIndex({'retroelement':[('Chr1',0,50)],'low_complexity':[('Chr1',50,100)]});fractions=index.fractions('Chr1',0,100)
  self.assertAlmostEqual(fractions['retroelement'],.5);self.assertAlmostEqual(fractions['low_complexity'],.5);self.assertAlmostEqual(index.overlap_fraction('Chr1',0,100),1.0)
 def test_repeatmasker_out_whitespace_parser(self):
  import tempfile
  with tempfile.TemporaryDirectory() as td:
   p=Path(td)/'rm.out';p.write_text('  123  1.0  0.0  0.0  Chr1  11  50  (0)  +  Gypsy  LTR/Gypsy  1  40  (0)  1\n')
   idx=self.m.load_repeat_index(p);self.assertAlmostEqual(idx.fractions('Chr1',10,50)['retroelement'],1.0)
 def test_processed_and_duplicated_evidence_confidence(self):
  processed={'alignment_source':'miniprot','query_coverage':.8,'frameshift_count':1,'stop_count':0,'candidate_intron_count':0,'polyA_run':10,'tsd_length':5,'strand':'+'}
  pe=self.m.subtype_evidence(processed,4,10000,1,'+',{'retroelement':.2},0);self.assertEqual(pe['subtype'],'processed_like');self.assertEqual(pe['processed_confidence'],'high');self.assertEqual(pe['intron_loss_count'],3)
  duplicated={'alignment_source':'miniprot','query_coverage':.8,'frameshift_count':1,'stop_count':0,'candidate_intron_count':2,'polyA_run':0,'tsd_length':0,'strand':'+'}
  de=self.m.subtype_evidence(duplicated,3,5000,1,'+',{},.8);self.assertEqual(de['subtype'],'duplicated_like');self.assertEqual(de['duplicated_confidence'],'high')
 def test_evidence_levels(self):
  high=self.m.evidence_level({'frameshift_count':1,'stop_count':0,'query_coverage':.8,'identity':.8,'query_left_unaligned_fraction':.1,'query_right_unaligned_fraction':.1})
  low=self.m.evidence_level({'frameshift_count':0,'stop_count':0,'query_coverage':.9,'identity':.8,'query_left_unaligned_fraction':.05,'query_right_unaligned_fraction':.05})
  self.assertEqual(high[0],'high'); self.assertEqual(low[0],'low')
if __name__=='__main__': unittest.main(verbosity=2)
