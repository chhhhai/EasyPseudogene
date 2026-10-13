#!/usr/bin/env python3
import importlib.util,sys,tempfile,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]; sys.path.insert(0,str(ROOT/'lib'))
def load():
 spec=importlib.util.spec_from_file_location('mmfast',ROOT/'lib'/'self_catalog_mmseqs_fast.py'); m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m
class MmseqsFastHelperTests(unittest.TestCase):
 @classmethod
 def setUpClass(cls): cls.m=load()
 def test_hardmask_merges_overlapping_cds(self):
  with tempfile.TemporaryDirectory() as td:
   td=Path(td); fa=td/'g.fa'; gff=td/'g.gff3'; out=td/'m.fa'
   fa.write_text('>Chr1\n'+'A'*300+'\n')
   gff.write_text('Chr1\tx\tCDS\t101\t200\t.\t+\t0\tParent=T1\nChr1\tx\tCDS\t151\t250\t.\t+\t0\tParent=T2\n')
   summary=self.m.hardmask_cds(fa,gff,out)
   seq=''.join(x.strip() for x in out.read_text().splitlines() if not x.startswith('>'))
   self.assertEqual(summary['merged_cds_intervals'],1); self.assertEqual(summary['cds_masked_bp'],150)
   self.assertEqual(seq[100:250],'N'*150); self.assertEqual(seq[:100],'A'*100)
 def test_ncbi_protein_id_gene_map(self):
  with tempfile.TemporaryDirectory() as td:
   g=Path(td)/'a.gff3';g.write_text('Chr1\tx\tCDS\t1\t30\t.\t+\t0\tprotein_id=XP_1.1;gene=LOC1\n')
   self.assertEqual(self.m.load_gff_protein_gene_map(g)['XP_1.1'],'LOC1')
 def test_longest_isoform_uses_generic_gene_field(self):
  with tempfile.TemporaryDirectory() as td:
   td=Path(td); fa=td/'p.fa'; out=td/'o.fa'
   fa.write_text('>tx1 gene=G1\nMAA\n>tx2 gene=G1\nMAAAAA\n>tx3 gene=G2\nMA\n')
   s=self.m.select_longest_isoforms(fa,out)
   self.assertEqual(s['longest_proteins'],2); text=out.read_text(); self.assertIn('>tx2 ',text); self.assertNotIn('>tx1 ',text)
if __name__=='__main__': unittest.main(verbosity=2)
