#!/usr/bin/env python3
"""MMseqs2 translated-search front end for the EasyPseudogene locus catalog."""
from __future__ import annotations
import argparse, hashlib, json, os, shutil, subprocess, sys, tempfile
from collections import defaultdict
from pathlib import Path
HERE=Path(__file__).resolve().parent
if str(HERE) not in sys.path: sys.path.insert(0,str(HERE))
from catalog_utils import atomic_text_writer, infer_gene_id, merge_half_open_intervals, read_fasta, write_fasta
BLAST_FAST=HERE/'self_catalog_blast_fast.py'

def file_fingerprint(path,block=1024*1024):
 path=Path(path); h=hashlib.sha256()
 with path.open('rb') as fh:
  while True:
   chunk=fh.read(block)
   if not chunk: break
   h.update(chunk)
 return {'path':str(path.resolve()),'size':path.stat().st_size,'sha256':h.hexdigest()}

def tool_version(executable):
 try:
  p=subprocess.run([executable,'version'],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=30)
  return p.stdout.strip().splitlines()[0] if p.stdout.strip() else f'exit={p.returncode}'
 except Exception as exc: return f'unavailable:{exc}'

def load_gff_protein_gene_map(gff):
 mapping={}
 if not gff:return mapping
 with Path(gff).open(errors='replace') as handle:
  for line in handle:
   if not line.strip() or line.startswith('#'):continue
   fields=line.rstrip('\n').split('\t')
   if len(fields)<9 or fields[2]!='CDS':continue
   attrs={}
   for part in fields[8].split(';'):
    if '=' in part:
     key,value=part.split('=',1);attrs[key]=value
   protein=attrs.get('protein_id') or attrs.get('Name');gene=attrs.get('gene') or attrs.get('gene_id')
   if protein and gene:mapping[protein]=gene
 return mapping

def select_longest_isoforms(proteins,out_path,gene_map=None):
 best={}; n_in=0; ambiguous=0; gene_map=gene_map or {}
 for name,desc,seq in read_fasta(Path(proteins)):
  n_in+=1; gene=gene_map.get(name) or infer_gene_id(name,desc)
  if not gene: ambiguous+=1; gene=name
  current=best.get(gene)
  if current is None or len(seq)>len(current[2]): best[gene]=(name,f'{name} gene_symbol={gene} {desc}',seq)
 write_fasta(best.values(),Path(out_path))
 return {'input_proteins':n_in,'longest_proteins':len(best),'ambiguous_gene_ids':ambiguous,'path':str(out_path)}

def parse_cds_intervals(gff,chroms):
 cds=defaultdict(list)
 with Path(gff).open(errors='replace') as handle:
  for line in handle:
   if not line.strip() or line.startswith('#'): continue
   parts=line.rstrip('\n').split('\t')
   if len(parts)<5 or parts[2]!='CDS' or parts[0] not in chroms: continue
   try: start=int(parts[3])-1; end=int(parts[4])
   except ValueError: continue
   if end>start: cds[parts[0]].append((start,end))
 return {chrom:merge_half_open_intervals(v) for chrom,v in cds.items()}

def hardmask_cds(genome,gff,out_path):
 records=list(read_fasta(Path(genome))); cds=parse_cds_intervals(gff,{n for n,_,_ in records})
 output=[]; masked_bp=0
 for name,desc,seq in records:
  arr=bytearray(seq.encode('ascii'))
  for start,end in cds.get(name,[]):
   s=max(0,start); e=min(len(arr),end)
   if e>s:
    masked_bp += sum(base not in (ord('N'),ord('n')) for base in arr[s:e])
    arr[s:e]=b'N'*(e-s)
  output.append((name,desc,arr.decode('ascii')))
 write_fasta(output,Path(out_path))
 return {'input_genome':str(genome),'masked_genome':str(out_path),'genome_records':len(records),'merged_cds_intervals':sum(len(v) for v in cds.values()),'cds_masked_bp':masked_bp}

def run(cmd,log_path=None):
 if log_path:
  Path(log_path).parent.mkdir(parents=True,exist_ok=True)
  with Path(log_path).open('w') as log: p=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,text=True)
 else: p=subprocess.run(cmd)
 if p.returncode!=0: raise SystemExit(f'Command failed with exit {p.returncode}: {" ".join(cmd)}')

def parse_args(argv=None):
 ap=argparse.ArgumentParser(description=__doc__)
 for option in ['proteins','genome','workdir','out-tsv','out-gff3','summary-json']: ap.add_argument('--'+option,required=True)
 ap.add_argument('--gff'); ap.add_argument('--threads',type=int,default=16)
 ap.add_argument('--mmseqs-sensitivity',type=float,default=5.7); ap.add_argument('--mmseqs-max-seqs',type=int,default=1000); ap.add_argument('--mmseqs-evalue',type=float,default=0.1); ap.add_argument('--mmseqs-translation-mode',type=int,default=0); ap.add_argument('--mmseqs-search-type',type=int,default=2)
 ap.add_argument('--keep-longest-isoform',action='store_true'); ap.add_argument('--sqlite-stream',action='store_true'); ap.add_argument('--mask-cds-before-search',action='store_true'); ap.add_argument('--reuse-m8'); ap.add_argument('--no-cache',action='store_true')
 ap.add_argument('--max-evalue',type=float,default=0.1); ap.add_argument('--min-bitscore',type=float,default=40); ap.add_argument('--min-align-aa',type=int,default=30); ap.add_argument('--cluster-gap',type=int,default=1000); ap.add_argument('--enforce-query-order',action='store_true'); ap.add_argument('--query-order-tolerance-aa',type=int,default=20); ap.add_argument('--min-cluster-score',type=float,default=50); ap.add_argument('--score-mode',choices=['raw','nonredundant'],default='raw'); ap.add_argument('--chain-mode',choices=['proximity','dag'],default='proximity'); ap.add_argument('--dag-max-query-gap',type=int,default=500); ap.add_argument('--max-cds-overlap',type=float,default=1.01); ap.add_argument('--overlap-filter-basis',choices=['block','outer','both'],default='both'); ap.add_argument('--source-overlap-cutoff',type=float,default=1.01); ap.add_argument('--max-candidates-per-query',type=int,default=0)
 ap.add_argument('--merge-loci',action='store_true'); ap.add_argument('--merge-gap',type=int,default=2000); ap.add_argument('--merge-same-strand',action='store_true'); ap.add_argument('--merge-coordinate-mode',choices=['union','best','consensus'],default='union'); ap.add_argument('--max-merged-span',type=int,default=0)
 args=ap.parse_args(argv)
 if args.threads < 1: ap.error('--threads must be >= 1')
 if args.mmseqs_sensitivity <= 0: ap.error('--mmseqs-sensitivity must be > 0')
 if args.mmseqs_max_seqs < 1: ap.error('--mmseqs-max-seqs must be >= 1')
 if args.mmseqs_evalue <= 0 or args.max_evalue <= 0: ap.error('e-values must be > 0')
 for name in ['min_align_aa','cluster_gap','query_order_tolerance_aa','min_cluster_score','max_cds_overlap','source_overlap_cutoff','max_candidates_per_query','merge_gap','max_merged_span']:
  if getattr(args,name) < 0: ap.error(f'--{name.replace("_","-")} must be >= 0')
 return args

def main(argv=None):
 args=parse_args(argv); workdir=Path(args.workdir); workdir.mkdir(parents=True,exist_ok=True)
 proteins_in=Path(args.proteins); genome_in=Path(args.genome); gff=Path(args.gff) if args.gff else None
 if args.mask_cds_before_search and not gff: raise SystemExit('--mask-cds-before-search requires --gff')
 manifest_core={'pipeline_schema_version':2,'inputs':{'proteins':file_fingerprint(proteins_in),'genome':file_fingerprint(genome_in),'gff':file_fingerprint(gff) if gff else None},'mmseqs_version':tool_version('mmseqs'),'parameters':{'sensitivity':args.mmseqs_sensitivity,'max_seqs':args.mmseqs_max_seqs,'evalue':args.mmseqs_evalue,'translation_mode':args.mmseqs_translation_mode,'search_type':args.mmseqs_search_type,'keep_longest':args.keep_longest_isoform,'mask_cds':args.mask_cds_before_search}}
 cache_key=hashlib.sha256(json.dumps(manifest_core,sort_keys=True).encode()).hexdigest()[:16]
 manifest_path=workdir/f'search_manifest.{cache_key}.json'
 m8=Path(args.reuse_m8) if args.reuse_m8 else workdir/f'mmseqs_translated.{cache_key}.m8'
 cache_hit=False
 if not args.reuse_m8 and not args.no_cache and m8.exists() and manifest_path.exists():
  try:
   cached=json.loads(manifest_path.read_text())
   cache_hit=cached.get('search')==manifest_core and cached.get('artifact')==file_fingerprint(m8)
  except Exception: cache_hit=False
 prep={}; proteins=proteins_in; genome=genome_in
 if not cache_hit and not args.reuse_m8:
  if args.keep_longest_isoform:
   proteins=workdir/f'query.longest_isoform.{cache_key}.fa'; prep['longest_isoform']=select_longest_isoforms(proteins_in,proteins,load_gff_protein_gene_map(gff))
  if args.mask_cds_before_search:
   genome=workdir/f'target.CDS_hardmasked.{cache_key}.fa'; prep['cds_hardmask']=hardmask_cds(genome_in,gff,genome)
  tmp=Path(tempfile.mkdtemp(prefix=f'mmseqs_tmp.{cache_key}.',dir=str(workdir)))
  fd,tmp_name=tempfile.mkstemp(prefix=f'mmseqs_translated.{cache_key}.tmp.',dir=str(workdir)); os.close(fd); tmp_m8=Path(tmp_name); tmp_m8.unlink()
  cmd=['mmseqs','easy-search',str(proteins),str(genome),str(tmp_m8),str(tmp),'--search-type',str(args.mmseqs_search_type),'--translation-mode',str(args.mmseqs_translation_mode),'-s',str(args.mmseqs_sensitivity),'--threads',str(args.threads),'--max-seqs',str(args.mmseqs_max_seqs),'-e',str(args.mmseqs_evalue),'--format-output','query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qlen']
  try:
   run(cmd,workdir/f'mmseqs.easy-search.{cache_key}.log'); os.replace(tmp_m8,m8)
  finally:
   shutil.rmtree(tmp,ignore_errors=True)
   if tmp_m8.exists(): tmp_m8.unlink()
  artifact=file_fingerprint(m8)
  with atomic_text_writer(manifest_path) as out: out.write(json.dumps({'search':manifest_core,'command':cmd,'artifact':artifact},indent=2)+'\n')
 elif cache_hit:
  prep['cache']={'hit':True,'manifest':str(manifest_path)}
 # Original query FASTA contains the same IDs and reliable gene metadata; qlen is also embedded in new M8 output.
 blast_cmd=[sys.executable,str(BLAST_FAST),'--blast6',str(m8),'--query-fasta',str(proteins_in),'--max-evalue',str(args.max_evalue),'--min-bitscore',str(args.min_bitscore),'--min-align-aa',str(args.min_align_aa),'--cluster-gap',str(args.cluster_gap),'--query-order-tolerance-aa',str(args.query_order_tolerance_aa),'--min-cluster-score',str(args.min_cluster_score),'--score-mode',args.score_mode,'--chain-mode',args.chain_mode,'--dag-max-query-gap',str(args.dag_max_query_gap),'--max-cds-overlap',str(args.max_cds_overlap),'--overlap-filter-basis',args.overlap_filter_basis,'--source-overlap-cutoff',str(args.source_overlap_cutoff),'--max-candidates-per-query',str(args.max_candidates_per_query),'--merge-gap',str(args.merge_gap),'--merge-coordinate-mode',args.merge_coordinate_mode,'--max-merged-span',str(args.max_merged_span),'--out-tsv',args.out_tsv,'--out-gff3',args.out_gff3,'--summary-json',args.summary_json]
 if gff: blast_cmd += ['--gff',str(gff)]
 if args.sqlite_stream: blast_cmd += ['--sqlite-stream']
 if args.enforce_query_order: blast_cmd += ['--enforce-query-order']
 if args.merge_loci: blast_cmd += ['--merge-loci']
 if args.merge_same_strand: blast_cmd += ['--merge-same-strand']
 run(blast_cmd,workdir/'blast_fast_postprocess.log')
 summary_path=Path(args.summary_json); summary=json.loads(summary_path.read_text())
 summary['mmseqs_fast']={'m8':str(m8),'m8_bytes':m8.stat().st_size,'threads':args.threads,'cache_hit':cache_hit,'manifest':str(manifest_path),'search':manifest_core,'prep':prep}
 with atomic_text_writer(summary_path) as out: out.write(json.dumps(summary,indent=2,sort_keys=True)+'\n')
 print(json.dumps(summary,indent=2,sort_keys=True))
if __name__=='__main__': main()
