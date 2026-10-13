#!/usr/bin/env python3
"""Build same-species pseudogene loci from BLAST6/MMseqs M8 HSPs."""
from __future__ import annotations
import argparse, csv, json, os, re, sqlite3, sys, tempfile
from urllib.parse import unquote
from collections import defaultdict
from pathlib import Path
HERE=Path(__file__).resolve().parent
if str(HERE) not in sys.path: sys.path.insert(0,str(HERE))
from catalog_utils import IntervalIndex, atomic_text_writer, load_protein_metadata, merge_half_open_intervals

def parse_args(argv=None):
 ap=argparse.ArgumentParser(description=__doc__)
 ap.add_argument('--blast6',type=Path,required=True); ap.add_argument('--sqlite-stream',action='store_true'); ap.add_argument('--keep-sqlite',action='store_true')
 ap.add_argument('--out-tsv',type=Path,required=True); ap.add_argument('--out-gff3',type=Path,required=True); ap.add_argument('--summary-json',type=Path,required=True)
 ap.add_argument('--gff',type=Path); ap.add_argument('--query-fasta',type=Path)
 ap.add_argument('--max-evalue',type=float,default=0.1); ap.add_argument('--min-bitscore',type=float,default=40); ap.add_argument('--min-align-aa',type=int,default=30)
 ap.add_argument('--cluster-gap',type=int,default=1000); ap.add_argument('--min-cluster-score',type=float,default=50); ap.add_argument('--score-mode',choices=['raw','nonredundant'],default='raw');ap.add_argument('--chain-mode',choices=['proximity','dag'],default='proximity');ap.add_argument('--dag-max-query-gap',type=int,default=500); ap.add_argument('--enforce-query-order',action='store_true'); ap.add_argument('--query-order-tolerance-aa',type=int,default=20)
 ap.add_argument('--max-cds-overlap',type=float,default=1.01); ap.add_argument('--overlap-filter-basis',choices=['block','outer','both'],default='both'); ap.add_argument('--source-overlap-cutoff',type=float,default=1.01)
 ap.add_argument('--max-candidates-per-query',type=int,default=0)
 ap.add_argument('--merge-loci',action='store_true'); ap.add_argument('--merge-gap',type=int,default=0); ap.add_argument('--merge-same-strand',action='store_true')
 ap.add_argument('--merge-coordinate-mode',choices=['union','best','consensus'],default='union'); ap.add_argument('--max-merged-span',type=int,default=0)
 args=ap.parse_args(argv)
 for name in ['min_align_aa','cluster_gap','query_order_tolerance_aa','merge_gap','max_candidates_per_query','max_merged_span']:
  if getattr(args,name)<0: ap.error(f'--{name.replace("_","-")} must be >= 0')
 for name in ['max_cds_overlap','source_overlap_cutoff']:
  if getattr(args,name)<0: ap.error(f'--{name.replace("_","-")} must be >= 0')
 return args

def parse_attrs(text):
 out={}
 for part in text.split(';'):
  part=part.strip()
  if not part: continue
  if '=' in part:
   k,v=part.split('=',1)
  else:
   match=re.match(r'([^\s]+)\s+["]?(.+?)["]?$',part)
   if not match: continue
   k,v=match.groups()
  out[k.strip()]=unquote(v.strip().strip('"'))
 return out

def aliases(value):
 out=[]
 for token in re.split(r'[,|]',value or ''):
  token=token.strip()
  if not token: continue
  out.extend([token,re.sub(r'^gene:','',token),re.sub(r'\.\d+$','',token),re.sub(r'\.\d+$','',re.sub(r'^gene:','',token))])
 return list(dict.fromkeys(x for x in out if x))

def load_gff_indexes(gff):
 cds=defaultdict(list); genes=defaultdict(list); records=[]
 if not gff: return IntervalIndex(),genes
 with Path(gff).open(errors='replace') as fh:
  for line in fh:
   if not line.strip() or line.startswith('#'): continue
   f=line.rstrip('\n').split('\t')
   if len(f)<9: continue
   try: start0=int(f[3])-1; end0=int(f[4])
   except ValueError: continue
   attrs=parse_attrs(f[8]); records.append((f[0],f[2],start0,end0,f[6],attrs))
   if f[2]=='CDS': cds[f[0]].append((start0,end0))
 # Index direct gene features first.
 gene_alias_to_loci=defaultdict(list)
 for chrom,feature,start0,end0,strand,attrs in records:
  if feature.lower() not in {'gene','pseudogene','ncrna_gene'}: continue
  loc=(chrom,start0,end0,strand); ids=[]
  for key in ['ID','Name','gene_id','gene','locus_tag']:
   ids.extend(aliases(attrs.get(key,'')))
  for key in ids:
   if loc not in gene_alias_to_loci[key]: gene_alias_to_loci[key].append(loc)
 # Add transcript aliases and Parent/gene_id mappings. This supports GFF3 and GTF.
 for chrom,feature,start0,end0,strand,attrs in records:
  if feature.lower() not in {'mrna','transcript','rna','lnc_rna','ncrna'}: continue
  loc=(chrom,start0,end0,strand); own=[]; parent=[]
  for key in ['ID','Name','transcript_id']:
   own.extend(aliases(attrs.get(key,'')))
  for key in ['Parent','gene_id','gene','locus_tag']:
   parent.extend(aliases(attrs.get(key,'')))
  parent_loci=[]
  for key in parent: parent_loci.extend(gene_alias_to_loci.get(key,[]))
  use_loci=parent_loci or [loc]
  for key in dict.fromkeys(own+parent):
   for source_loc in use_loci:
    if source_loc not in gene_alias_to_loci[key]: gene_alias_to_loci[key].append(source_loc)
 # NCBI-style CDS records often expose protein_id and gene directly.
 for chrom,feature,start0,end0,strand,attrs in records:
  if feature!='CDS': continue
  protein_ids=aliases(attrs.get('protein_id',''))+aliases(attrs.get('Name',''))
  parent_keys=aliases(attrs.get('gene',''))+aliases(attrs.get('gene_id',''))+aliases(attrs.get('Parent',''))
  parent_loci=[]
  for key in parent_keys: parent_loci.extend(gene_alias_to_loci.get(key,[]))
  for protein_id in dict.fromkeys(protein_ids):
   for source_loc in parent_loci:
    if source_loc not in gene_alias_to_loci[protein_id]: gene_alias_to_loci[protein_id].append(source_loc)
 return IntervalIndex(cds),gene_alias_to_loci

def read_blast6(path,max_evalue,min_bitscore,min_align_aa,query_meta=None):
 rows=[]; total=0; malformed=0
 query_meta=query_meta or {}
 with Path(path).open(errors='replace') as fh:
  for line in fh:
   if not line.strip() or line.startswith('#'): continue
   total+=1; f=line.rstrip('\n').split('\t')
   if len(f)<12: malformed+=1; continue
   try:
    q,s=f[0],f[1]; pid=float(f[2]); alen=int(float(f[3])); qs=int(float(f[6])); qe=int(float(f[7])); ss=int(float(f[8])); se=int(float(f[9])); ev=float(f[10]); bits=float(f[11])
    qlen_supplied=len(f)>12 and f[12] not in {'','.'}
    qlen=int(float(f[12])) if qlen_supplied else int(query_meta.get(q,{}).get('length') or max(qs,qe))
   except (ValueError,TypeError): malformed+=1; continue
   if ev>max_evalue or bits<min_bitscore or alen<min_align_aa: continue
   rows.append(dict(query=q,chrom=s,start=min(ss,se),end=max(ss,se),strand='+' if ss<=se else '-',pident=pid,aln_len=alen,qstart=min(qs,qe),qend=max(qs,qe),qlen=qlen,evalue=ev,bitscore=bits,qlen_supplied=qlen_supplied))
 return rows,total,malformed

def read_blast6_sqlite(path,max_evalue,min_bitscore,min_align_aa,query_meta,sqlite_path):
 query_meta=query_meta or {}; total=malformed=passed=0; qlen_supplied_any=False; connection=sqlite3.connect(str(sqlite_path)); cursor=connection.cursor()
 cursor.execute('DROP TABLE IF EXISTS hsps'); cursor.execute('CREATE TABLE hsps(query TEXT,chrom TEXT,start INT,end INT,strand TEXT,pident REAL,aln_len INT,qstart INT,qend INT,qlen INT,evalue REAL,bitscore REAL,qlen_supplied INT)')
 batch=[]
 with Path(path).open(errors='replace') as fh:
  for line in fh:
   if not line.strip() or line.startswith('#'): continue
   total+=1; f=line.rstrip('\n').split('\t')
   if len(f)<12: malformed+=1; continue
   try:
    q,chrom=f[0],f[1];pid=float(f[2]);alen=int(float(f[3]));qs=int(float(f[6]));qe=int(float(f[7]));ss=int(float(f[8]));se=int(float(f[9]));ev=float(f[10]);bits=float(f[11]);supplied=len(f)>12 and f[12] not in {'','.'};qlen=int(float(f[12])) if supplied else int(query_meta.get(q,{}).get('length') or max(qs,qe))
   except (ValueError,TypeError): malformed+=1; continue
   qlen_supplied_any=qlen_supplied_any or supplied
   if ev>max_evalue or bits<min_bitscore or alen<min_align_aa: continue
   batch.append((q,chrom,min(ss,se),max(ss,se),'+' if ss<=se else '-',pid,alen,min(qs,qe),max(qs,qe),qlen,ev,bits,int(supplied)));passed+=1
   if len(batch)>=50000: cursor.executemany('INSERT INTO hsps VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)',batch);connection.commit();batch=[]
 if batch: cursor.executemany('INSERT INTO hsps VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)',batch);connection.commit()
 cursor.execute('CREATE INDEX idx_hsp_order ON hsps(chrom,query,strand,start,end)');connection.commit()
 def iterator():
  try:
   for row in connection.execute('SELECT query,chrom,start,end,strand,pident,aln_len,qstart,qend,qlen,evalue,bitscore,qlen_supplied FROM hsps ORDER BY chrom,query,strand,start,end'):
    yield dict(query=row[0],chrom=row[1],start=row[2],end=row[3],strand=row[4],pident=row[5],aln_len=row[6],qstart=row[7],qend=row[8],qlen=row[9],evalue=row[10],bitscore=row[11],qlen_supplied=bool(row[12]))
  finally: connection.close()
 return iterator(),{'blast6_rows':total,'malformed_rows':malformed,'passed_hsps':passed,'qlen_supplied_any':qlen_supplied_any}

def nonredundant_bitscore(hsps):
 covered=[];score=0.0
 for hsp in sorted(hsps,key=lambda x:(x['bitscore'],x['aln_len']),reverse=True):
  start0=min(hsp['qstart'],hsp['qend'])-1;end0=max(hsp['qstart'],hsp['qend']);length=max(1,end0-start0)
  overlap=sum(max(0,min(end0,e)-max(start0,s)) for s,e in covered);novel=max(0,length-overlap)
  score+=hsp['bitscore']*(novel/length);covered=merge_half_open_intervals(covered+[(start0,end0)])
 return score,covered

def build_cluster(cur,key,start,end,score_mode='raw'):
 best=max(cur,key=lambda r:(r['bitscore'],r['aln_len'],-r['evalue']));raw_score=sum(r['bitscore'] for r in cur);nonred_score,qblocks=nonredundant_bitscore(cur);raw_aln=sum(r['aln_len'] for r in cur);qlen=max(r['qlen'] for r in cur);qcov=sum(e-s for s,e in qblocks)/max(1,qlen);mean=sum(r['pident']*r['aln_len'] for r in cur)/max(1,raw_aln);tblocks0=merge_half_open_intervals([(r['start']-1,r['end']) for r in cur]);tblocks=[(s+1,e) for s,e in tblocks0];evidence_len=sum(e-s for s,e in tblocks0)
 return dict(query=best['query'],chrom=key[0],strand=key[2],start=start,end=end,cluster_size=len(cur),cluster_score=raw_score if score_mode=='raw' else nonred_score,raw_bitscore_sum=raw_score,nonredundant_score=nonred_score,best_bitscore=best['bitscore'],total_aligned_aa=sum(e-s for s,e in qblocks),raw_aligned_aa=raw_aln,mean_identity=mean,query_length=qlen,query_coverage=qcov,query_blocks=[(s+1,e) for s,e in qblocks],evidence_blocks=tblocks,evidence_length=evidence_len,evalue=min(r['evalue'] for r in cur),candidate_type='blast_homology_pseudogene_candidate')


def dag_chain_hsps(hsps,strand,query_tolerance=20,max_query_gap=500):
 """Split one proximity group into deterministic collinear HSP chains by DP."""
 remaining=list(range(len(hsps)));chains=[]
 while remaining:
  ordered=sorted(remaining,key=lambda i:(hsps[i]['start'],hsps[i]['end'],hsps[i]['qstart'],hsps[i]['qend']))
  score={};previous={}
  for position,index in enumerate(ordered):
   current=hsps[index];score[index]=current['bitscore'];previous[index]=None
   for prior_index in ordered[:position]:
    prior=hsps[prior_index];target_overlap=max(0,prior['end']-current['start']+1)
    if target_overlap>max(30,int(.25*min(prior['end']-prior['start']+1,current['end']-current['start']+1))):continue
    if strand=='+':
     if current['qstart']+query_tolerance<prior['qstart']:continue
     query_gap=max(0,current['qstart']-prior['qend'])
    else:
     if current['qend']>prior['qend']+query_tolerance:continue
     query_gap=max(0,prior['qstart']-current['qend'])
    if query_gap>max_query_gap:continue
    genomic_gap=max(0,current['start']-prior['end']);candidate=score[prior_index]+current['bitscore']-0.001*genomic_gap-0.10*query_gap-0.25*target_overlap
    if candidate>score[index]:score[index]=candidate;previous[index]=prior_index
  end=max(ordered,key=lambda i:(score[i],hsps[i]['bitscore']));chain=[]
  while end is not None:chain.append(end);end=previous[end]
  chain.reverse();chains.append([hsps[i] for i in chain]);used=set(chain);remaining=[i for i in remaining if i not in used]
 return chains

def iter_clusters_sorted(rows,gap,enforce_query_order=False,query_order_tolerance_aa=20,score_mode='raw',chain_mode='proximity',dag_max_query_gap=500):
 cur=[];key=None;start=end=None
 def emit(group,group_key):
  if chain_mode=='dag':
   chains=dag_chain_hsps(group,group_key[2],query_order_tolerance_aa,dag_max_query_gap)
   for chain_index,chain in enumerate(chains,1):
    result=build_cluster(chain,group_key,min(x['start'] for x in chain),max(x['end'] for x in chain),score_mode);result['chain_mode']='dag';result['chain_index']=chain_index;yield result
  else:
   result=build_cluster(group,group_key,min(x['start'] for x in group),max(x['end'] for x in group),score_mode);result['chain_mode']='proximity';result['chain_index']=1;yield result
 for row in rows:
  k=(row['chrom'],row['query'],row['strand']);order_break=False
  if cur and k==key and enforce_query_order:
   previous=cur[-1];order_break=(row['qstart']+query_order_tolerance_aa<previous['qstart']) if row['strand']=='+' else (row['qstart']>previous['qstart']+query_order_tolerance_aa)
  if not cur or k!=key or row['start']>end+gap or order_break:
   if cur:yield from emit(cur,key)
   cur=[row];key=k;start=row['start'];end=row['end']
  else:cur.append(row);start=min(start,row['start']);end=max(end,row['end'])
 if cur:yield from emit(cur,key)

def cluster_hsps(rows,gap,enforce_query_order=False,query_order_tolerance_aa=20,score_mode='raw',chain_mode='proximity',dag_max_query_gap=500):
 rows=sorted(rows,key=lambda r:(r['chrom'],r['query'],r['strand'],r['start'],r['end']))
 return list(iter_clusters_sorted(rows,gap,enforce_query_order,query_order_tolerance_aa,score_mode,chain_mode,dag_max_query_gap))

def source_overlap_fraction(candidate,gene_loci,query_meta,use_blocks=True):
 q=candidate['query'];meta=query_meta.get(q,{});gene=str(meta.get('gene_id') or re.sub(r'\.\d+$','',q));keys=aliases(gene)+aliases(q)
 blocks=[(s-1,e) for s,e in candidate.get('evidence_blocks',[])] if use_blocks and candidate.get('evidence_blocks') else [(int(candidate['start'])-1,int(candidate['end']))];denominator=max(1,sum(e-s for s,e in blocks));best=0
 for key in dict.fromkeys(keys):
  overlap=0
  for bstart,bend in blocks:
   block_best=0
   for chrom,s,e,_strand in gene_loci.get(key,[]):
    if chrom==candidate['chrom']:block_best=max(block_best,max(0,min(bend,e)-max(bstart,s)))
   overlap+=block_best
  best=max(best,overlap)
 return min(1.0,best/denominator)

def annotate_overlap(candidate,cds_index,gene_loci,query_meta):
 start0=int(candidate['start'])-1;end0=int(candidate['end']);blocks=[(s-1,e) for s,e in candidate.get('evidence_blocks',[])] or [(start0,end0)];den=max(1,sum(e-s for s,e in blocks))
 candidate['cds_overlap_frac']=min(1.0,sum(cds_index.overlap_bp(candidate['chrom'],s,e) for s,e in blocks)/den);candidate['source_overlap_frac']=source_overlap_fraction(candidate,gene_loci,query_meta,True)
 candidate['cds_outer_overlap_frac']=cds_index.overlap_fraction(candidate['chrom'],start0,end0);candidate['source_outer_overlap_frac']=source_overlap_fraction(candidate,gene_loci,query_meta,False);candidate['evidence_length']=den;candidate['evidence_block_count']=len(blocks)
 return candidate

def weighted_median(values_weights):
 ordered=sorted((float(v),max(0.0,float(w))) for v,w in values_weights)
 total=sum(w for _v,w in ordered)
 if not ordered: return 0
 if total<=0: return int(ordered[len(ordered)//2][0])
 acc=0
 for value,weight in ordered:
  acc+=weight
  if acc>=total/2: return int(round(value))
 return int(round(ordered[-1][0]))

def merge_loci(clusters,gap=0,same_strand=False,coordinate_mode='union',cds_index=None,gene_loci=None,query_meta=None,max_merged_span=0):
 if not clusters: return []
 cds_index=cds_index or IntervalIndex(); gene_loci=gene_loci or {}; query_meta=query_meta or {}
 arr=sorted(clusters,key=lambda c:(c['chrom'],c['strand'] if same_strand else '.',c['start'],c['end'],c['query']))
 merged=[]; cur=[]; cur_chr=cur_str=None; start=end=None
 def flush():
  nonlocal cur,start,end
  if not cur: return
  by_parent=defaultdict(list)
  for candidate in cur: by_parent[candidate['query']].append(candidate)
  def parent_rank(item):
   _query,members=item
   return (max(c['query_coverage'] for c in members),max(c['best_bitscore'] for c in members),sum(c['cluster_score'] for c in members),sum(c['total_aligned_aa'] for c in members))
  parent,selected=max(by_parent.items(),key=parent_rank)
  best=max(selected,key=lambda c:(c['query_coverage'],c['best_bitscore'],c['cluster_score'],c['total_aligned_aa']))
  if coordinate_mode=='best': ns,ne=best['start'],best['end']
  elif coordinate_mode=='consensus':
   ns=weighted_median((c['start'],c['cluster_score']) for c in cur); ne=weighted_median((c['end'],c['cluster_score']) for c in cur)
   if ne<ns: ns,ne=min(c['start'] for c in cur),max(c['end'] for c in cur)
  else: ns,ne=start,end
  queries=[q for q,_members in sorted(by_parent.items(),key=parent_rank,reverse=True)]
  aln=sum(c['total_aligned_aa'] for c in selected); out=dict(best)
  support_bonus=min(100.0,10.0*max(0,len(queries)-1))
  parent_score=sum(c['cluster_score'] for c in selected)
  locus_tblocks0=merge_half_open_intervals([(s-1,e) for c in cur for s,e in c.get('evidence_blocks',[(c['start'],c['end'])])]);selected_qblocks0=merge_half_open_intervals([(s-1,e) for c in selected for s,e in c.get('query_blocks',[])]);nonred_aa=sum(e-s for s,e in selected_qblocks0) if selected_qblocks0 else sum(c['total_aligned_aa'] for c in selected);query_length=max((c.get('query_length') or round(c['total_aligned_aa']/max(1e-9,c.get('query_coverage',1))) for c in selected),default=nonred_aa)
  out.update(query=parent,start=ns,end=ne,parent_gene=parent,cluster_size=sum(c['cluster_size'] for c in selected),cluster_score=parent_score+support_bonus,raw_bitscore_sum=sum(c.get('raw_bitscore_sum',c['cluster_score']) for c in selected),nonredundant_score=sum(c.get('nonredundant_score',c['cluster_score']) for c in selected),best_bitscore=max(c['best_bitscore'] for c in selected),total_aligned_aa=nonred_aa,raw_aligned_aa=sum(c.get('raw_aligned_aa',c['total_aligned_aa']) for c in selected),mean_identity=sum(c['mean_identity']*c['total_aligned_aa'] for c in selected)/max(1,aln),query_length=query_length,query_coverage=nonred_aa/max(1,query_length),query_blocks=[(s+1,e) for s,e in selected_qblocks0],evidence_blocks=[(s+1,e) for s,e in locus_tblocks0],evidence_length=sum(e-s for s,e in locus_tblocks0),evalue=min(c['evalue'] for c in selected),merged_locus_size=len(cur),alternative_parents=','.join(queries[1:25]))
  annotate_overlap(out,cds_index,gene_loci,query_meta); merged.append(out); cur=[]
 for c in arr:
  kchr=c['chrom']; kstr=c['strand'] if same_strand else '.'
  can_join=bool(cur) and kchr==cur_chr and kstr==cur_str and c['start']<=end+gap
  if can_join and max_merged_span and max(end,c['end'])-min(start,c['start'])+1>max_merged_span: can_join=False
  if not can_join:
   flush(); cur=[c]; cur_chr=kchr; cur_str=kstr; start=c['start']; end=c['end']
  else:
   cur.append(c); start=min(start,c['start']); end=max(end,c['end'])
 flush(); return merged

def gff_escape(value):
 return str(value).replace('%','%25').replace(';','%3B').replace('=','%3D').replace(',','%2C').replace(' ','%20')

def main(argv=None):
 args=parse_args(argv); query_meta=load_protein_metadata(args.query_fasta) if args.query_fasta else {}
 cds_index,gene_loci=load_gff_indexes(args.gff); sqlite_path=None
 if args.sqlite_stream:
  args.out_tsv.parent.mkdir(parents=True,exist_ok=True);fd,name=tempfile.mkstemp(prefix='easypseudogene_hsps.',suffix='.sqlite',dir=str(args.out_tsv.parent));os.close(fd);sqlite_path=Path(name)
  hsp_iter,stats=read_blast6_sqlite(args.blast6,args.max_evalue,args.min_bitscore,args.min_align_aa,query_meta,sqlite_path)
  raw_iter=iter_clusters_sorted(hsp_iter,args.cluster_gap,args.enforce_query_order,args.query_order_tolerance_aa,args.score_mode,args.chain_mode,args.dag_max_query_gap);total=stats['blast6_rows'];malformed=stats['malformed_rows'];passed_hsps=stats['passed_hsps'];true_qlen=bool(args.query_fasta or stats['qlen_supplied_any'])
 else:
  hsps,total,malformed=read_blast6(args.blast6,args.max_evalue,args.min_bitscore,args.min_align_aa,query_meta);passed_hsps=len(hsps);true_qlen=bool(args.query_fasta or any(r.get('qlen_supplied') for r in hsps));raw_iter=iter(cluster_hsps(hsps,args.cluster_gap,args.enforce_query_order,args.query_order_tolerance_aa,args.score_mode,args.chain_mode,args.dag_max_query_gap))
 clusters=[];removed_score=removed_cds=removed_source=raw_count=0
 for c in raw_iter:
  raw_count+=1
  if c['cluster_score']<args.min_cluster_score: removed_score+=1; continue
  annotate_overlap(c,cds_index,gene_loci,query_meta)
  cds_filter={'block':c['cds_overlap_frac'],'outer':c['cds_outer_overlap_frac'],'both':max(c['cds_overlap_frac'],c['cds_outer_overlap_frac'])}[args.overlap_filter_basis]
  source_filter={'block':c['source_overlap_frac'],'outer':c['source_outer_overlap_frac'],'both':max(c['source_overlap_frac'],c['source_outer_overlap_frac'])}[args.overlap_filter_basis]
  if cds_filter>args.max_cds_overlap: removed_cds+=1; continue
  if source_filter>=args.source_overlap_cutoff: removed_source+=1; continue
  clusters.append(c)
 if args.max_candidates_per_query:
  byq=defaultdict(list)
  for c in clusters: byq[c['query']].append(c)
  clusters=[]
  for arr in byq.values():
   arr.sort(key=lambda c:(c['cluster_score'],c['best_bitscore'],c['total_aligned_aa']),reverse=True); clusters.extend(arr[:args.max_candidates_per_query])
 pre_merge=len(clusters)
 if args.merge_loci:
  clusters=merge_loci(clusters,args.merge_gap,args.merge_same_strand,args.merge_coordinate_mode,cds_index,gene_loci,query_meta,args.max_merged_span)
  kept=[]
  for c in clusters:
   if c['cds_overlap_frac']>args.max_cds_overlap: removed_cds+=1; continue
   if c['source_overlap_frac']>=args.source_overlap_cutoff: removed_source+=1; continue
   kept.append(c)
  clusters=kept
 clusters.sort(key=lambda c:(c['chrom'],c['start'],c['end'],c['query']))
 fields=['candidate_id','query_id','parent_gene','chr','start','end','strand','candidate_type','length','cluster_size','cluster_score','raw_bitscore_sum','nonredundant_score','best_bitscore','total_aligned_aa','raw_aligned_aa','mean_identity','query_length','query_coverage','cds_overlap_frac','source_overlap_frac','cds_outer_overlap_frac','source_outer_overlap_frac','evidence_block_count','evidence_length','outer_span_length','merged_locus_size','alternative_parents','evalue','confidence_score','evidence_tags']
 with atomic_text_writer(args.out_tsv) as out:
  w=csv.DictWriter(out,fieldnames=fields,delimiter='\t'); w.writeheader()
  for i,c in enumerate(clusters,1):
   conf=min(100.0,c['cluster_score']/3.0)
   w.writerow(dict(candidate_id=f'BF{i:06d}',query_id=c['query'],parent_gene=c['query'],chr=c['chrom'],start=c['start'],end=c['end'],strand=c['strand'],candidate_type=c['candidate_type'],length=c['end']-c['start']+1,cluster_size=c['cluster_size'],cluster_score=round(c['cluster_score'],3),raw_bitscore_sum=round(c.get('raw_bitscore_sum',c['cluster_score']),3),nonredundant_score=round(c.get('nonredundant_score',c['cluster_score']),3),best_bitscore=round(c['best_bitscore'],3),total_aligned_aa=c['total_aligned_aa'],raw_aligned_aa=c.get('raw_aligned_aa',c['total_aligned_aa']),mean_identity=round(c['mean_identity'],3),query_length=c.get('query_length',0),query_coverage=round(c['query_coverage'],4),cds_overlap_frac=round(c['cds_overlap_frac'],6),source_overlap_frac=round(c['source_overlap_frac'],6),cds_outer_overlap_frac=round(c.get('cds_outer_overlap_frac',0),6),source_outer_overlap_frac=round(c.get('source_outer_overlap_frac',0),6),evidence_block_count=c.get('evidence_block_count',1),evidence_length=c.get('evidence_length',c['end']-c['start']+1),outer_span_length=c['end']-c['start']+1,merged_locus_size=c.get('merged_locus_size',1),alternative_parents=c.get('alternative_parents',''),evalue=c['evalue'],confidence_score=round(conf,3),evidence_tags=f'homology_cluster;cluster_size={c["cluster_size"]};cds_overlap={c["cds_overlap_frac"]:.3f};source_overlap={c["source_overlap_frac"]:.3f}'))
 with atomic_text_writer(args.out_gff3) as out:
  out.write('##gff-version 3\n')
  for i,c in enumerate(clusters,1):
   attrs=';'.join([f'ID=BF{i:06d}',f'Name=BF{i:06d}_{gff_escape(c["query"])}',f'Derives_from={gff_escape(c["query"])}',f'candidate_type={c["candidate_type"]}',f'cluster_size={c["cluster_size"]}',f'cluster_score={c["cluster_score"]:.2f}',f'query_coverage={c["query_coverage"]:.3f}',f'cds_overlap={c["cds_overlap_frac"]:.3f}',f'source_overlap={c["source_overlap_frac"]:.3f}'])
   out.write(f'{c["chrom"]}\tEasyPseudogene\tpseudogene\t{c["start"]}\t{c["end"]}\t{c["cluster_score"]:.2f}\t{c["strand"]}\t.\t{attrs}\n')
 summary=dict(blast6_rows=total,malformed_rows=malformed,passed_hsps=passed_hsps,raw_clusters=raw_count,removed_low_score=removed_score,removed_cds_overlap=removed_cds,removed_source_locus=removed_source,pre_merge_candidates=pre_merge,output_candidates=len(clusters),merge_loci=args.merge_loci,merge_gap=args.merge_gap,merge_same_strand=args.merge_same_strand,merge_coordinate_mode=args.merge_coordinate_mode,max_merged_span=args.max_merged_span,max_evalue=args.max_evalue,min_bitscore=args.min_bitscore,min_align_aa=args.min_align_aa,cluster_gap=args.cluster_gap,enforce_query_order=args.enforce_query_order,query_order_tolerance_aa=args.query_order_tolerance_aa,min_cluster_score=args.min_cluster_score,score_mode=args.score_mode,overlap_filter_basis=args.overlap_filter_basis,max_cds_overlap=args.max_cds_overlap,source_overlap_cutoff=args.source_overlap_cutoff,true_query_coverage=true_qlen,sqlite_stream=args.sqlite_stream)
 with atomic_text_writer(args.summary_json) as out: out.write(json.dumps(summary,indent=2,sort_keys=True)+'\n')
 if sqlite_path and sqlite_path.exists() and not args.keep_sqlite: sqlite_path.unlink()
 print(json.dumps(summary,indent=2,sort_keys=True))
if __name__=='__main__': main()
