#!/usr/bin/env python3
import importlib.util,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('evaluate',ROOT/'scripts'/'evaluate_catalog.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
bspec=importlib.util.spec_from_file_location('bootstrap',ROOT/'scripts'/'bootstrap_locus_metrics.py');b=importlib.util.module_from_spec(bspec);bspec.loader.exec_module(b)
class EvaluateTests(unittest.TestCase):
 def test_one_prediction_cannot_match_two_truth_loci(self):
  truth=[('Chr1',100,200,'T1'),('Chr1',200,300,'T2')];by={'Chr1':[(100,200,'T1'),(200,300,'T2')]};result=m.evaluate([('Chr1',100,300)],truth,by,None)
  self.assertEqual(result['true_positives'],1);self.assertEqual(result['precision'],1);self.assertEqual(result['recall'],.5)
 def test_maximum_matching_avoids_greedy_under_count(self):
  matching=m.maximum_cardinality_matching([[0,1],[0]],2);self.assertEqual(len(matching),2)
 def test_bootstrap_quantile_interpolates(self):
  self.assertAlmostEqual(b.quantile([0,10],.25),2.5)
 def test_matched_pair_is_auditable(self):
  truth=[('Chr1',100,200,'T1')];by={'Chr1':truth};metrics,pairs=m.match_loci([('Chr1',120,180)],truth,by,None,['C1'])
  self.assertEqual(metrics['matching_policy'],'one_to_one_maximum_cardinality');self.assertEqual(pairs[0]['prediction_id'],'C1');self.assertEqual(pairs[0]['truth_id'],'T1');self.assertEqual(pairs[0]['overlap_bp'],60);self.assertAlmostEqual(pairs[0]['prediction_overlap_fraction'],1.0)
if __name__=='__main__':unittest.main()
