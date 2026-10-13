#!/usr/bin/env python3
import csv,importlib.util,tempfile,unittest,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"lib"))
spec=importlib.util.spec_from_file_location("tiers",ROOT/"lib"/"catalog_tiers.py");m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
class TierTests(unittest.TestCase):
 def test_explicit_tier_semantics(self):
  with tempfile.TemporaryDirectory() as td:
   root=Path(td);src=root/"v.tsv";fields=["candidate_id","alignment_source","disruption_evidence_level","pseudogene_subtype"]
   rows=[["A","miniprot","high","processed_like"],["B","miniprot","medium","fragmentary"],["C","mmseqs_fallback","medium","duplicated_like"]]
   with src.open("w",newline="") as h:
    w=csv.writer(h,delimiter="\t",lineterminator="\n");w.writerow(fields);w.writerows(rows)
   summary=m.write_tiers(src,root/"catalog",root/"summary.json")
   self.assertEqual(summary["tier_counts"],{"sensitive":3,"miniprot_supported":2,"strict_core_disruption":1,"subtype_candidates":2})
   self.assertTrue((root/"catalog.strict_core_disruption.tsv").exists())
if __name__=="__main__":unittest.main()
