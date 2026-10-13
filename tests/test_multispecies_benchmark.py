#!/usr/bin/env python3
import importlib.util,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];spec=importlib.util.spec_from_file_location("multi",ROOT/"scripts"/"multispecies_benchmark.py");m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
class MultiSpeciesTests(unittest.TestCase):
 def test_family_fold_is_stable(self):
  self.assertEqual(m.stable_fold("rice","F1",5),m.stable_fold("arabidopsis","F1",5));self.assertTrue(0<=m.stable_fold("rice","F1",5)<5)
 def test_logistic_learns_separable_fixture(self):
  x=[[-2],[-1],[1],[2]];y=[0,0,1,1];w=m.fit_logistic(x,y,epochs=1000);score=[m.sigmoid(w[0]+w[1]*r[0]) for r in x];self.assertGreater(m.auc(y,score),.99)
if __name__=="__main__":unittest.main()
