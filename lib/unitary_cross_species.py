#!/usr/bin/env python3
"""Conservative cross-species unitary-pseudogene evidence annotator.
Definitive candidates require no intact target copy, strong disruption, synteny,
and an intact outgroup ortholog. Missing gates produce explicit needs_* labels.
"""
from __future__ import annotations
import argparse,csv,json,os,subprocess,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
if str(HERE) not in sys.path:sys.path.insert(0,str(HERE))
from catalog_utils import IntervalIndex,atomic_text_writer

def read_best_proteome_hits(path):
 best={}
 with Path(path).open(errors='replace') as handle:
  for line in handle:
   if not line.strip() or line.startswith('#'):continue
   fields=line.rstrip('\n').split('\t')
   if len(fields)<12:continue
   try:
    query,target=fields[0],fields[1];identity=float(fields[2])/100;aligned=int(float(fields[3]));qstart=int(float(fields[4]));qend=int(float(fields[5]));evalue=float(fields[8]);bits=float(fields[9]);qlen=int(float(fields[10]));tlen=int(float(fields[11]));qcov=(abs(qend-qstart)+1)/max(1,qlen)
   except ValueError:continue
   hit={'target':target,'identity':identity,'aligned_aa':aligned,'query_coverage':qcov,'evalue':evalue,'bitscore':bits,'target_length':tlen};rank=(qcov,identity,bits)
   if query not in best or rank>best[query][0]:best[query]=(rank,hit)
 return {query:value[1] for query,value in best.items()}

def is_intact(hit,qcov,identity,evalue):return bool(hit and hit['query_coverage']>=qcov and hit['identity']>=identity and hit['evalue']<=evalue)
def classify_unitary(functional_copy_present,strong_disruption,synteny_supported,synteny_supplied,outgroup_intact=False,outgroup_supplied=False,genome_copy_present=False,genome_scan_supplied=False):
 if functional_copy_present or genome_copy_present:return 'not_unitary_functional_copy_present'
 if not strong_disruption:return 'unresolved_no_functional_copy_without_strong_disruption'
 if not genome_scan_supplied:return 'unitary_like_needs_genome_scan'
 if not synteny_supplied or not synteny_supported:return 'unitary_like_needs_synteny'
 if not outgroup_supplied or not outgroup_intact:return 'unitary_like_needs_outgroup'
 return 'unitary_candidate'
def load_synteny(path):
 if not path:return IntervalIndex(),False
 intervals={}
 with Path(path).open() as handle:
  for line in handle:
   if not line.strip() or line.startswith('#'):continue
   f=line.rstrip('\n').split('\t')
   try:intervals.setdefault(f[0],[]).append((int(f[1]),int(f[2])))
   except (ValueError,IndexError):continue
 return IntervalIndex(intervals),True
def run_search(reference,proteome,out_path,tmpdir,threads):
 temporary=out_path.with_name(out_path.name+f'.tmp.{os.getpid()}');cmd=['mmseqs','easy-search',str(reference),str(proteome),str(temporary),str(tmpdir),'--threads',str(threads),'-e','0.1','--format-output','query,target,pident,alnlen,qstart,qend,tstart,tend,evalue,bits,qlen,tlen'];process=subprocess.run(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
 if process.returncode:raise SystemExit(process.stderr[-2000:])
 os.replace(temporary,out_path)

def run_genome_search(reference,target_genome,out_paf,threads):
 temporary=out_paf.with_name(out_paf.name+f'.tmp.{os.getpid()}');cmd=['miniprot','-t',str(threads),'--outn','500','--outs','0.1',str(target_genome),str(reference)]
 with temporary.open('w') as out:process=subprocess.run(cmd,stdout=out,stderr=subprocess.PIPE,text=True)
 if process.returncode:raise SystemExit(process.stderr[-2000:])
 os.replace(temporary,out_paf)
def read_genome_hits(path):
 hits={}
 with Path(path).open(errors='replace') as handle:
  for line in handle:
   if not line.strip() or line.startswith('#'):continue
   f=line.rstrip('\n').split('\t')
   if len(f)<12:continue
   try:
    query,qlen,qstart,qend,target,tstart,tend=f[0],int(f[1]),int(f[2]),int(f[3]),f[5],int(f[7]),int(f[8]);matches,block=int(f[9]),int(f[10]);tags={x[:2]:x[5:] for x in f[12:] if len(x)>5 and x[2:5] in {':i:',':f:',':Z:'}};cg=tags.get('cg','');qcov=(qend-qstart)/max(1,qlen);identity=matches/max(1,block)
   except ValueError:continue
   hits.setdefault(query,[]).append({'target':target,'start':min(tstart,tend),'end':max(tstart,tend),'query_coverage':qcov,'identity':identity,'intact':qcov>=.8 and identity>=.5 and 'F' not in cg and 'S' not in cg})
 return hits
def intact_genome_copy(hits,row):
 chrom=row.get('chr','');start=int(row.get('start','1'))-1;end=int(row.get('end','0'))
 for hit in hits:
  overlap=max(0,min(end,hit['end'])-max(start,hit['start']));candidate_fraction=overlap/max(1,end-start);hit_fraction=overlap/max(1,hit['end']-hit['start'])
  if hit['target']==chrom and min(candidate_fraction,hit_fraction)>=.8:continue
  if hit['intact']:return hit
 return None

def main(argv=None):
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--candidates',required=True,type=Path);ap.add_argument('--reference-proteins',required=True,type=Path);ap.add_argument('--target-proteins',required=True,type=Path);ap.add_argument('--target-genome',type=Path);ap.add_argument('--reuse-target-genome-paf',type=Path);ap.add_argument('--outgroup-proteins',type=Path);ap.add_argument('--workdir',required=True,type=Path);ap.add_argument('--out-tsv',required=True,type=Path);ap.add_argument('--summary-json',required=True,type=Path);ap.add_argument('--synteny-bed',type=Path);ap.add_argument('--reuse-proteome-m8',type=Path);ap.add_argument('--reuse-outgroup-m8',type=Path);ap.add_argument('--threads',type=int,default=16);ap.add_argument('--functional-qcov',type=float,default=.8);ap.add_argument('--functional-identity',type=float,default=.5);ap.add_argument('--functional-evalue',type=float,default=1e-5);args=ap.parse_args(argv);args.workdir.mkdir(parents=True,exist_ok=True)
 target_m8=args.reuse_proteome_m8 or args.workdir/'reference_vs_target_proteome.m8'
 if not args.reuse_proteome_m8:run_search(args.reference_proteins,args.target_proteins,target_m8,args.workdir/'target_mmseqs_tmp',args.threads)
 target_hits=read_best_proteome_hits(target_m8);genome_scan_supplied=bool(args.target_genome or args.reuse_target_genome_paf);genome_paf=args.reuse_target_genome_paf;
 if args.target_genome and not genome_paf:genome_paf=args.workdir/'reference_vs_target_genome.paf';run_genome_search(args.reference_proteins,args.target_genome,genome_paf,args.threads)
 genome_hits=read_genome_hits(genome_paf) if genome_paf else {};outgroup_supplied=bool(args.outgroup_proteins or args.reuse_outgroup_m8);outgroup_hits={};outgroup_m8=args.reuse_outgroup_m8
 if args.outgroup_proteins and not outgroup_m8:outgroup_m8=args.workdir/'reference_vs_outgroup_proteome.m8';run_search(args.reference_proteins,args.outgroup_proteins,outgroup_m8,args.workdir/'outgroup_mmseqs_tmp',args.threads)
 if outgroup_m8:outgroup_hits=read_best_proteome_hits(outgroup_m8)
 synteny,synteny_supplied=load_synteny(args.synteny_bed);candidate_handle=args.candidates.open();reader=csv.DictReader(candidate_handle,delimiter='\t');rows=list(reader);candidate_handle.close();base_fields=reader.fieldnames or [];extra=['target_genome_functional_copy','target_genome_query_coverage','target_genome_identity','genome_scan_supplied','target_best_hit','target_functional_copy','target_functional_query_coverage','target_functional_identity','target_functional_evalue','outgroup_best_hit','outgroup_functional_copy','outgroup_functional_query_coverage','outgroup_functional_identity','outgroup_intact','synteny_supported','unitary_classification'];counts={}
 with atomic_text_writer(args.out_tsv) as out:
  writer=csv.DictWriter(out,fieldnames=base_fields+extra,delimiter='\t',lineterminator='\n');writer.writeheader()
  for row in rows:
   query=row.get('selected_parent') or row.get('query_id');target=target_hits.get(query);genome_hit=intact_genome_copy(genome_hits.get(query,[]),row);outgroup=outgroup_hits.get(query);present=is_intact(target,args.functional_qcov,args.functional_identity,args.functional_evalue);out_intact=is_intact(outgroup,args.functional_qcov,args.functional_identity,args.functional_evalue);syn=synteny.overlap_bp(row['chr'],int(row['start'])-1,int(row['end']))>0 if synteny_supplied else False;strong=row.get('disruption_evidence_level')=='high';label=classify_unitary(present,strong,syn,synteny_supplied,out_intact,outgroup_supplied,bool(genome_hit),genome_scan_supplied);counts[label]=counts.get(label,0)+1;writer.writerow({**row,'target_genome_functional_copy':genome_hit['target'] if genome_hit else '','target_genome_query_coverage':f'{genome_hit["query_coverage"]:.4f}' if genome_hit else '','target_genome_identity':f'{genome_hit["identity"]:.4f}' if genome_hit else '','genome_scan_supplied':int(genome_scan_supplied),'target_best_hit':target['target'] if target else '','target_functional_copy':target['target'] if present else '','target_functional_query_coverage':f'{target["query_coverage"]:.4f}' if target else '','target_functional_identity':f'{target["identity"]:.4f}' if target else '','target_functional_evalue':target['evalue'] if target else '','outgroup_best_hit':outgroup['target'] if outgroup else '','outgroup_functional_copy':outgroup['target'] if out_intact else '','outgroup_functional_query_coverage':f'{outgroup["query_coverage"]:.4f}' if outgroup else '','outgroup_functional_identity':f'{outgroup["identity"]:.4f}' if outgroup else '','outgroup_intact':int(out_intact),'synteny_supported':int(syn),'unitary_classification':label})
 summary={'candidates':len(rows),'classifications':counts,'target_proteome_hits':len(target_hits),'target_genome_queries':len(genome_hits),'genome_scan_supplied':genome_scan_supplied,'outgroup_hits':len(outgroup_hits),'synteny_supplied':synteny_supplied,'outgroup_supplied':outgroup_supplied,'thresholds':{'functional_qcov':args.functional_qcov,'functional_identity':args.functional_identity,'functional_evalue':args.functional_evalue}}
 with atomic_text_writer(args.summary_json) as out:out.write(json.dumps(summary,indent=2,sort_keys=True)+'\n')
 print(json.dumps(summary,indent=2,sort_keys=True))
if __name__=='__main__':main()
