#!/usr/bin/env python3
import importlib.util,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def load(name):
 spec=importlib.util.spec_from_file_location(name,ROOT/'scripts'/f'{name}.py');mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);return mod
class SyntenyHelperTests(unittest.TestCase):
 def test_paf_interval_merge_is_chromosome_specific(self):
  m=load('paf_to_synteny_bed');self.assertEqual(m.merge([('1',10,20),('1',22,30),('2',1,5)],2),[('1',10,30),('2',1,5)])
 def test_mcscan_block_merge(self):
  m=load('mcscanx_blocks_to_bed');self.assertEqual(m.merge([('1',5,10),('1',9,20),('1',30,40)]),[('1',5,20),('1',30,40)])
if __name__=='__main__':unittest.main()
