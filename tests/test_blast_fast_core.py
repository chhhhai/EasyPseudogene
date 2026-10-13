#!/usr/bin/env python3
import importlib.util, sys, tempfile, unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'lib'))
def load():
 spec=importlib.util.spec_from_file_location('blastfast',ROOT/'lib'/'self_catalog_blast_fast.py')
 m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m
class BlastFastCoreTests(unittest.TestCase):
 @classmethod
 def setUpClass(cls): cls.m=load()
 def test_gff_cds_overlap_is_union_and_gene_locus_loaded(self):
  with tempfile.TemporaryDirectory() as td:
   g=Path(td)/'a.gff3'; g.write_text('Chr1\tx\tgene\t101\t400\t.\t+\t.\tID=gene:G1;Name=G1\nChr1\tx\tCDS\t101\t200\t.\t+\t0\tParent=T1\nChr1\tx\tCDS\t151\t250\t.\t+\t0\tParent=T2\n')
   cds,genes=self.m.load_gff_indexes(g)
   self.assertEqual(cds.overlap_bp('Chr1',100,250),150)
   self.assertIn('G1',genes)
   self.assertEqual(genes['G1'][0][:3],('Chr1',100,400))
 def test_gtf_gene_and_transcript_attributes_are_indexed(self):
  with tempfile.TemporaryDirectory() as td:
   g=Path(td)/'a.gtf'; g.write_text('Chr1\tx\tgene\t101\t400\t.\t+\t.\tgene_id "G1"; gene_name "N1";\nChr1\tx\ttranscript\t101\t400\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n')
   _cds,genes=self.m.load_gff_indexes(g)
   self.assertIn('G1',genes); self.assertIn('T1',genes)

 def test_cluster_uses_true_nonredundant_query_coverage(self):
  rows=[dict(query='Q1',chrom='Chr1',start=100,end=200,strand='+',pident=50,aln_len=30,qstart=1,qend=30,qlen=100,evalue=1e-5,bitscore=60),dict(query='Q1',chrom='Chr1',start=250,end=350,strand='+',pident=50,aln_len=40,qstart=51,qend=90,qlen=100,evalue=1e-5,bitscore=70)]
  c=self.m.cluster_hsps(rows,1000)[0]
  self.assertAlmostEqual(c['query_coverage'],0.70)
 def test_dag_chain_splits_query_order_reset(self):
  base=dict(query='Q',chrom='Chr1',strand='+',pident=80,aln_len=50,qlen=100,evalue=1e-10,bitscore=100,qlen_supplied=True)
  rows=[dict(base,start=100,end=249,qstart=1,qend=50),dict(base,start=300,end=449,qstart=51,qend=100),dict(base,start=500,end=649,qstart=1,qend=50)]
  clusters=self.m.cluster_hsps(rows,1000,chain_mode='dag')
  self.assertEqual(len(clusters),2);self.assertEqual(sorted(x['cluster_size'] for x in clusters),[1,2]);self.assertTrue(all(x['chain_mode']=='dag' for x in clusters))
 def test_nonredundant_score_downweights_overlapping_query_hsps(self):
  rows=[dict(query='Q1',chrom='Chr1',start=100,end=200,strand='+',pident=50,aln_len=100,qstart=1,qend=100,qlen=120,evalue=1e-5,bitscore=100),dict(query='Q1',chrom='Chr1',start=120,end=220,strand='+',pident=50,aln_len=100,qstart=1,qend=100,qlen=120,evalue=1e-5,bitscore=90)]
  raw=self.m.cluster_hsps(rows,1000,score_mode='raw')[0];nonred=self.m.cluster_hsps(rows,1000,score_mode='nonredundant')[0]
  self.assertEqual(raw['cluster_score'],190);self.assertEqual(nonred['cluster_score'],100);self.assertEqual(nonred['total_aligned_aa'],100)
 def test_source_overlap_detected_by_query_gene_mapping(self):
  c=dict(query='Q1',chrom='Chr1',start=101,end=400)
  genes={'G1':[('Chr1',100,400,'+')]}; qmap={'Q1':{'gene_id':'G1','length':100}}
  self.assertGreaterEqual(self.m.source_overlap_fraction(c,genes,qmap),0.99)
 def test_merge_keeps_metrics_from_selected_parent(self):
  clusters=[dict(query='Q1',chrom='Chr1',strand='+',start=100,end=200,cluster_size=2,cluster_score=100,best_bitscore=90,total_aligned_aa=80,mean_identity=60,query_length=100,query_coverage=.8,evalue=1e-8,candidate_type='x',cds_overlap_frac=0,source_overlap_frac=0),dict(query='Q2',chrom='Chr1',strand='+',start=150,end=250,cluster_size=20,cluster_score=1000,best_bitscore=500,total_aligned_aa=900,mean_identity=99,query_length=1000,query_coverage=.1,evalue=1e-20,candidate_type='x',cds_overlap_frac=0,source_overlap_frac=0)]
  out=self.m.merge_loci(clusters,gap=0,same_strand=True,coordinate_mode='union',cds_index=self.m.IntervalIndex(),gene_loci={},query_meta={})[0]
  self.assertEqual(out['query'],'Q1'); self.assertEqual(out['best_bitscore'],90); self.assertEqual(out['total_aligned_aa'],80); self.assertAlmostEqual(out['mean_identity'],60)

 def test_merge_recomputes_cds_overlap_on_union(self):
  cds=self.m.IntervalIndex({'Chr1':[(180,260)]})
  clusters=[dict(query='Q1',chrom='Chr1',strand='+',start=101,end=200,cluster_size=1,cluster_score=60,best_bitscore=60,total_aligned_aa=30,mean_identity=50,query_coverage=0.3,evalue=1e-5,candidate_type='x',cds_overlap_frac=0,source_overlap_frac=0),dict(query='Q2',chrom='Chr1',strand='+',start=220,end=300,cluster_size=1,cluster_score=70,best_bitscore=70,total_aligned_aa=40,mean_identity=55,query_coverage=0.4,evalue=1e-6,candidate_type='x',cds_overlap_frac=0,source_overlap_frac=0)]
  out=self.m.merge_loci(clusters,gap=30,same_strand=True,coordinate_mode='union',cds_index=cds,gene_loci={},query_meta={})[0]
  self.assertEqual((out['start'],out['end']),(101,300))
  self.assertAlmostEqual(out['cds_overlap_frac'],61/181)
  self.assertAlmostEqual(out['cds_outer_overlap_frac'],80/200)
if __name__=='__main__': unittest.main(verbosity=2)
