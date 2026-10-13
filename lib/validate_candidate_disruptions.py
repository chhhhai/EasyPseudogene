#!/usr/bin/env python3
"""Targeted family-aware miniprot validation and pseudogene subtype annotation."""
from __future__ import annotations
import argparse,csv,json,re,subprocess,sys
from collections import defaultdict
from pathlib import Path
HERE=Path(__file__).resolve().parent
if str(HERE) not in sys.path: sys.path.insert(0,str(HERE))
from catalog_utils import IntervalIndex,atomic_text_writer,read_fasta,write_fasta
from self_catalog_blast_fast import aliases,load_gff_indexes

def parse_tags(fields):
 tags={}
 for item in fields[12:]:
  parts=item.split(':',2)
  if len(parts)!=3: continue
  key,typ,val=parts
  if typ=='i':
   try: val=int(val)
   except ValueError: pass
  tags[key]=val
 return tags

def load_candidates(path):
 with Path(path).open() as fh: return list(csv.DictReader(fh,delimiter='\t'))

def parse_attrs(text):
 out={}
 for part in text.split(';'):
  part=part.strip()
  if not part: continue
  if '=' in part: key,value=part.split('=',1)
  else:
   m=re.match(r'([^\s]+)\s+"?(.+?)"?$',part)
   if not m: continue
   key,value=m.groups()
  out[key]=value.strip('"')
 return out

def load_parent_exon_counts(gff):
 counts=defaultdict(int)
 if not gff: return counts
 with Path(gff).open(errors='replace') as fh:
  for line in fh:
   if not line.strip() or line.startswith('#'): continue
   f=line.rstrip('\n').split('\t')
   if len(f)<9 or f[2]!='CDS': continue
   attrs=parse_attrs(f[8]); parents=[]
   for key in ('Parent','transcript_id','protein_id'):
    parents.extend(x.strip().replace('transcript:','') for x in attrs.get(key,'').split(',') if x.strip())
   for parent in set(parents): counts[parent]+=1
 return counts

def alternative_queries(row,max_alternatives):
 initial=row['query_id']; values=[initial]
 for query in (row.get('alternative_parents') or '').split(','):
  query=query.strip()
  if query and query not in values: values.append(query)
 return values[:1+max_alternatives]

def longest_run(sequence,base):
 runs=re.findall(f'{base}+',sequence.upper()); return max((len(x) for x in runs),default=0)

def prepare_inputs(candidates,genome_path,protein_path,workdir,base_flank,max_flank,max_alternatives):
 genome={name:seq for name,_desc,seq in read_fasta(Path(genome_path))}
 proteins={name:(name,desc,seq) for name,desc,seq in read_fasta(Path(protein_path))}
 windows=[]; mapping={}; missing_chr=missing_query=0; wanted={}; dynamic_flanks=[]
 for row in candidates:
  cid=row['candidate_id']; chrom=row['chr']; initial=row['query_id']
  if chrom not in genome: missing_chr+=1; continue
  expected=[q for q in alternative_queries(row,max_alternatives) if q in proteins]
  if not expected: missing_query+=1; continue
  protein_aa=max(len(proteins[q][2]) for q in expected); candidate_len=int(row['end'])-int(row['start'])+1
  flank=min(max_flank,max(base_flank,int(protein_aa*3*0.25),int(candidate_len*0.20))); dynamic_flanks.append(flank)
  start=max(1,int(row['start'])-flank); end=min(len(genome[chrom]),int(row['end'])+flank); target='EPWIN_'+cid; sequence=genome[chrom][start-1:end]
  windows.append((target,f'{target} candidate={cid} chr={chrom} offset={start}',sequence))
  mapping[target]={'candidate_id':cid,'initial_query':initial,'expected_queries':expected,'chrom':chrom,'offset1':start,'window_end':end,'candidate_start':int(row['start']),'candidate_end':int(row['end']),'sequence':sequence}
  for query in expected: wanted[query]=proteins[query]
 winfa=workdir/'candidate_windows.fa'; qfa=workdir/'candidate_family_proteins.fa'; write_fasta(windows,winfa); write_fasta(wanted.values(),qfa)
 return winfa,qfa,mapping,{'windows':len(windows),'queries':len(wanted),'missing_chrom':missing_chr,'missing_query':missing_query,'mean_dynamic_flank':sum(dynamic_flanks)/max(1,len(dynamic_flanks)),'max_dynamic_flank':max(dynamic_flanks,default=0)}

def poly_a_evidence(sequence,target_start,target_end,strand,width=100):
 if strand=='+': downstream=sequence[target_end:target_end+width].upper(); base='A'
 else: downstream=sequence[max(0,target_start-width):target_start].upper(); base='T'
 return longest_run(downstream,base),downstream.count(base)/max(1,len(downstream))

def target_site_duplication(sequence,target_start,target_end,min_len=4,max_len=20):
 left=sequence[max(0,target_start-max_len):target_start].upper();right=sequence[target_end:target_end+max_len].upper()
 for length in range(min(max_len,len(left),len(right)),min_len-1,-1):
  if left[-length:]==right[:length] and set(left[-length:])<={"A","C","G","T"}:return length,left[-length:]
 return 0,""

def classify_repeat_label(label):
 text=str(label).lower()
 if any(token in text for token in ('low_complexity','low complexity')):return 'low_complexity'
 if any(token in text for token in ('simple_repeat','simple repeat','satellite','microsatellite')):return 'simple_repeat'
 if any(token in text for token in ('line','sine','ltr','retrotrans','retropos','gypsy','copia')):return 'retroelement'
 if any(token in text for token in ('dna_transpos','helitron','hat','mutator','mudr','tc1','mariner','cacta','en-spm')):return 'dna_transposon'
 return 'unknown'

class RepeatEvidenceIndex:
 def __init__(self,category_intervals=None):
  category_intervals=category_intervals or {};self.categories={};all_by_chrom=defaultdict(list)
  for name,intervals in category_intervals.items():
   by_chrom=defaultdict(list)
   for chrom,start,end in intervals:by_chrom[chrom].append((start,end));all_by_chrom[chrom].append((start,end))
   self.categories[name]=IntervalIndex(by_chrom)
  self.total=IntervalIndex(all_by_chrom)
 def overlap_bp(self,chrom,start,end):return self.total.overlap_bp(chrom,start,end)
 def overlap_fraction(self,chrom,start,end):return self.total.overlap_fraction(chrom,start,end)
 def fractions(self,chrom,start,end):
  length=max(1,end-start);return {name:index.overlap_bp(chrom,start,end)/length for name,index in self.categories.items()}

def load_repeat_index(path):
 categories=defaultdict(list)
 if not path:return RepeatEvidenceIndex()
 with Path(path).open() as handle:
  for line in handle:
   if not line.strip() or line.startswith('#'):continue
   f=line.rstrip('\n').split('\t')
   try:
    whitespace=line.split()
    if len(f)<3 and len(whitespace)>=11 and whitespace[0].isdigit():
     chrom=whitespace[4];start=int(whitespace[5])-1;end=int(whitespace[6]);label=' '.join(whitespace[9:11])
    elif len(f)>=9 and not f[3].isdigit():
     chrom=f[0];start=int(f[3])-1;end=int(f[4]);label=' '.join([f[2],f[8]])
    else:
     chrom=f[0];start=int(f[1]);end=int(f[2]);label=' '.join(f[3:])
   except (ValueError,IndexError):continue
   if end>start:categories[classify_repeat_label(label)].append((chrom,start,end))
 return RepeatEvidenceIndex(categories)

def read_best_family_paf(path,mapping):
 best={}; rows=expected=exact_rows=0
 with Path(path).open(errors='replace') as fh:
  for line in fh:
   if not line.strip(): continue
   rows+=1; f=line.rstrip('\n').split('\t')
   if len(f)<12 or f[5] not in mapping: continue
   target=f[5]; meta=mapping[target]
   if f[0] not in meta['expected_queries']: continue
   expected+=1
   if f[0]==meta['initial_query']: exact_rows+=1
   try:
    qlen=int(f[1]); qs=int(f[2]); qe=int(f[3]); ts=int(f[7]); te=int(f[8]); nmatch=int(f[9]); alen=int(f[10]); mapq=int(f[11])
   except ValueError: continue
   tags=parse_tags(f); fs=int(tags.get('fs',0) or 0); st=int(tags.get('st',0) or 0); score=int(tags.get('AS',0) or 0); cg=str(tags.get('cg',''))
   introns=len(re.findall(r'\d+[NUV]',cg)); qcov=(qe-qs)/max(1,qlen); ident=nmatch/max(1,alen); polyrun,polyfrac=poly_a_evidence(meta['sequence'],min(ts,te),max(ts,te),f[4]);tsd_len,tsd_seq=target_site_duplication(meta['sequence'],min(ts,te),max(ts,te))
   hit={'query_id':f[0],'initial_query_id':meta['initial_query'],'parent_changed':f[0]!=meta['initial_query'],'target_id':target,'query_coverage':qcov,'identity':ident,'aligned_length':alen,'mapq':mapq,'frameshift_count':fs,'stop_count':st,'alignment_score':score,'strand':f[4],'candidate_intron_count':introns,'polyA_run':polyrun,'polyA_fraction':polyfrac,'tsd_length':tsd_len,'tsd_sequence':tsd_seq,'qstart':qs,'qend':qe,'qlen':qlen,'alignment_source':'miniprot','target_start0':ts,'target_end0':te,'validated_genomic_start':meta.get('offset1',1)+ts,'validated_genomic_end':meta.get('offset1',1)+te-1,'match_fraction':ident,'query_left_unaligned_fraction':qs/max(1,qlen),'query_right_unaligned_fraction':max(0,qlen-qe)/max(1,qlen)}
   hit['structural_tags_complete']=int('fs' in tags and 'st' in tags)
   cid=meta['candidate_id']; rank=(qcov,ident,score,alen,1 if f[0]==meta['initial_query'] else 0)
   if cid not in best or rank>best[cid][0]: best[cid]=(rank,hit)
 return {cid:value[1] for cid,value in best.items()},{'paf_rows':rows,'family_pair_rows':expected,'exact_parent_rows':exact_rows,'validated_candidates':len(best)}

def read_mmseqs_family_fallback(path,mapping,existing_hits):
 by_key=defaultdict(list)
 for target,meta in mapping.items():
  if meta['candidate_id'] in existing_hits: continue
  for query in meta['expected_queries']: by_key[(meta['chrom'],query)].append((target,meta))
 best={};rows=family_rows=0
 with Path(path).open(errors='replace') as fh:
  for line in fh:
   if not line.strip() or line.startswith('#'): continue
   rows+=1;f=line.rstrip('\n').split('\t')
   if len(f)<12: continue
   key=(f[1],f[0]); candidates=by_key.get(key)
   if not candidates: continue
   try: pid=float(f[2])/100;alen=int(float(f[3]));qs=int(float(f[6]));qe=int(float(f[7]));ss=int(float(f[8]));se=int(float(f[9]));bits=float(f[11]);qlen=int(float(f[12])) if len(f)>12 else max(qs,qe)
   except ValueError: continue
   hs,he=min(ss,se),max(ss,se)
   for target,meta in candidates:
    if he<meta['candidate_start'] or hs>meta['candidate_end']: continue
    family_rows+=1;qcov=(max(qs,qe)-min(qs,qe)+1)/max(1,qlen);query=f[0]
    hit={'query_id':query,'initial_query_id':meta['initial_query'],'parent_changed':query!=meta['initial_query'],'target_id':target,'query_coverage':qcov,'identity':pid,'aligned_length':alen,'mapq':0,'frameshift_count':0,'stop_count':0,'alignment_score':bits,'strand':'+' if ss<=se else '-','candidate_intron_count':0,'polyA_run':0,'polyA_fraction':0.0,'tsd_length':0,'tsd_sequence':'','qstart':min(qs,qe),'qend':max(qs,qe),'qlen':qlen,'alignment_source':'mmseqs_fallback','validated_genomic_start':hs,'validated_genomic_end':he,'match_fraction':pid,'query_left_unaligned_fraction':min(qs,qe)/max(1,qlen),'query_right_unaligned_fraction':max(0,qlen-max(qs,qe))/max(1,qlen)}
    cid=meta['candidate_id'];rank=(qcov,pid,bits,alen,1 if query==meta['initial_query'] else 0)
    if cid not in best or rank>best[cid][0]:best[cid]=(rank,hit)
 return {cid:value[1] for cid,value in best.items()},{'fallback_m8_rows':rows,'fallback_family_rows':family_rows,'fallback_candidates':len(best)}

def evidence_level(hit):
 if not hit: return 'unsupported','no_family_parent_alignment'
 if hit.get('alignment_source')=='mmseqs_fallback':
  return 'medium','mmseqs_family_fallback'
 disruptions=hit['frameshift_count']+hit['stop_count']; qcov=hit['query_coverage']; ident=hit['identity']; left=hit.get('query_left_unaligned_fraction',1); right=hit.get('query_right_unaligned_fraction',1)
 if disruptions>0 and qcov>=0.50 and left<=0.20 and right<=0.20: return 'high',f'core_disruption_in_near_complete_alignment;fs={hit["frameshift_count"]};stop={hit["stop_count"]}'
 if disruptions>0 and qcov>=0.20: return 'medium',f'disruption_in_partial_or_terminal_alignment;fs={hit["frameshift_count"]};stop={hit["stop_count"]}'
 if disruptions>0: return 'medium','weak_partial_alignment_with_disruption'
 if qcov<0.35: return 'medium','fragmentary_parent_alignment'
 if qcov>=0.80 and left<=0.20 and right<=0.20: return 'low','intact_orf_like_alignment'
 if ident<0.35: return 'medium','low_match_fraction_parent_alignment'
 return 'medium','partial_parent_alignment'

def classify_subtype(hit,parent_exons):
 if not hit: return 'unsupported'
 if hit.get('alignment_source')=='mmseqs_fallback': return 'fragmentary' if hit['query_coverage']<0.35 else 'homology_candidate'
 qcov=hit['query_coverage']; disruptions=hit['frameshift_count']+hit['stop_count']; candidate_exons=hit['candidate_intron_count']+1
 if qcov<0.35: return 'fragmentary'
 if qcov>=0.80 and disruptions==0: return 'intact_orf_like'
 if parent_exons>=2 and candidate_exons==1 and qcov>=0.50: return 'processed_like'
 if candidate_exons>=2 and disruptions>0: return 'duplicated_like'
 if disruptions>0: return 'disrupted_unclassified'
 return 'homology_candidate'

def source_distance(row,selected_parent,gene_loci):
 loci=[]
 for key in aliases(selected_parent):loci.extend(gene_loci.get(key,[]))
 if not loci:return '',-1,0,''
 same=[x for x in loci if x[0]==row['chr']]
 if not same:return loci[0][0],-1,0,loci[0][3]
 start0=int(row['start'])-1;end0=int(row['end'])
 def distance(locus):
  _chrom,start,end,_strand=locus;return 0 if max(start0,start)<min(end0,end) else start-end0 if start>=end0 else start0-end
 best=min(same,key=distance);return row['chr'],max(0,distance(best)),1,best[3]

def load_interval_bed(path):
 by_chrom=defaultdict(list)
 if not path:return IntervalIndex()
 with Path(path).open() as handle:
  for line in handle:
   if not line.strip() or line.startswith('#'):continue
   f=line.rstrip('\n').split('\t')
   try:by_chrom[f[0]].append((int(f[1]),int(f[2])))
   except (ValueError,IndexError):continue
 return IntervalIndex(by_chrom)

def subtype_evidence(hit,parent_exons,source_dist,same_source,source_strand,repeat_fractions,synteny_fraction):
 if not hit:return {'subtype':'unsupported','candidate_exons':0,'exon_ratio':0.0,'intron_loss_count':0,'processed_score':0.0,'processed_confidence':'none','duplicated_score':0.0,'duplicated_confidence':'none','strand_match':''}
 candidate_exons=hit.get('candidate_intron_count',0)+1;qcov=hit.get('query_coverage',0.0);disruptions=hit.get('frameshift_count',0)+hit.get('stop_count',0);exon_ratio=min(parent_exons,candidate_exons)/max(1,max(parent_exons,candidate_exons));intron_loss=max(0,parent_exons-candidate_exons);strand_match=int(bool(source_strand) and hit.get('strand')==source_strand)
 processed=0.0
 if parent_exons>=2 and candidate_exons==1:processed+=2
 if intron_loss>=1:processed+=min(2,intron_loss)*0.5
 if qcov>=0.50:processed+=1
 if hit.get('polyA_run',0)>=8:processed+=1
 if hit.get('tsd_length',0)>=4:processed+=1
 if repeat_fractions.get('retroelement',0)>=0.10:processed+=0.5
 if disruptions>0:processed+=0.5
 if max(repeat_fractions.get('simple_repeat',0),repeat_fractions.get('low_complexity',0))>=0.50:processed-=1
 duplicated=0.0
 if candidate_exons>=2:duplicated+=2
 if exon_ratio>=0.50:duplicated+=1
 if disruptions>0:duplicated+=1
 if qcov>=0.50:duplicated+=1
 if synteny_fraction>=0.20:duplicated+=2
 if same_source and source_dist>=1000:duplicated+=0.5
 subtype=classify_subtype(hit,parent_exons)
 return {'subtype':subtype,'candidate_exons':candidate_exons,'exon_ratio':exon_ratio,'intron_loss_count':intron_loss,'processed_score':max(0,processed),'processed_confidence':(('high' if processed>=5 else 'medium' if processed>=3 else 'low') if subtype=='processed_like' else 'none'),'duplicated_score':duplicated,'duplicated_confidence':(('high' if duplicated>=6 and synteny_fraction>=0.20 else 'medium' if duplicated>=4 else 'low') if subtype=='duplicated_like' else 'none'),'strand_match':strand_match}

def write_validated(candidates,hits,parent_exon_counts,out_path,repeat_index=None,gene_loci=None,synteny_index=None):
 gene_loci=gene_loci or {}
 repeat_index=repeat_index or RepeatEvidenceIndex();synteny_index=synteny_index or IntervalIndex()
 base_fields=list(candidates[0].keys()) if candidates else []
 extra=['alignment_source','validator_status','validated_genomic_start','validated_genomic_end','validated_match_fraction','query_left_unaligned_fraction','query_right_unaligned_fraction','disruption_evidence_level','pseudogene_subtype','selected_parent','initial_parent','parent_changed','parent_exon_count','candidate_exon_count','validated_query_coverage','validated_identity','validated_aligned_length','validated_mapq','validated_frameshift_count','validated_stop_count','validated_alignment_score','validated_strand','polyA_run','polyA_fraction','tsd_length','tsd_sequence','repeat_overlap_fraction','repeat_retroelement_fraction','repeat_dna_transposon_fraction','repeat_simple_fraction','repeat_low_complexity_fraction','repeat_unknown_fraction','repeat_primary_class','repeat_supportive_fraction','repeat_penalty_fraction','te_associated','organellar_contig','source_chrom','source_distance','same_chromosome_as_source','source_strand','candidate_source_strand_match','synteny_overlap_fraction','parent_candidate_exon_ratio','intron_loss_count','processed_evidence_score','processed_confidence','duplicated_evidence_score','duplicated_confidence','validator_evidence']
 counts=defaultdict(int); subtype_counts=defaultdict(int)
 extra.append('structural_tags_complete')
 with atomic_text_writer(Path(out_path)) as out:
  w=csv.DictWriter(out,fieldnames=base_fields+extra,delimiter='\t'); w.writeheader()
  for row in candidates:
   hit=hits.get(row['candidate_id']); level,evidence=evidence_level(hit); selected=hit['query_id'] if hit else ''; parent_exons=parent_exon_counts.get(selected,0);counts[level]+=1
   repeat_fraction=repeat_index.overlap_fraction(row['chr'],int(row['start'])-1,int(row['end']));repeat_fractions=repeat_index.fractions(row['chr'],int(row['start'])-1,int(row['end']));primary_repeat=max(repeat_fractions,key=repeat_fractions.get) if repeat_fractions and max(repeat_fractions.values())>0 else '';supportive=repeat_fractions.get('retroelement',0.0);penalty=max(repeat_fractions.get('simple_repeat',0.0),repeat_fractions.get('low_complexity',0.0));organellar=int(str(row['chr']).lower() in {'mt','pt','chrm','chrc','mitochondrion','chloroplast'} or 'mitochond' in str(row['chr']).lower() or 'chloroplast' in str(row['chr']).lower());source_chrom,source_dist,same_source,source_strand=source_distance(row,selected,gene_loci);synteny_fraction=synteny_index.overlap_fraction(row['chr'],int(row['start'])-1,int(row['end']));subev=subtype_evidence(hit,parent_exons,source_dist,same_source,source_strand,repeat_fractions,synteny_fraction);subtype=subev['subtype'];subtype_counts[subtype]+=1
   values={key:'' for key in extra}; values.update(repeat_overlap_fraction=f'{repeat_fraction:.4f}',repeat_retroelement_fraction=f"{repeat_fractions.get('retroelement',0.0):.4f}",repeat_dna_transposon_fraction=f"{repeat_fractions.get('dna_transposon',0.0):.4f}",repeat_simple_fraction=f"{repeat_fractions.get('simple_repeat',0.0):.4f}",repeat_low_complexity_fraction=f"{repeat_fractions.get('low_complexity',0.0):.4f}",repeat_unknown_fraction=f"{repeat_fractions.get('unknown',0.0):.4f}",repeat_primary_class=primary_repeat,repeat_supportive_fraction=f'{supportive:.4f}',repeat_penalty_fraction=f'{penalty:.4f}',te_associated=int(repeat_fraction>=0.5),organellar_contig=organellar,source_chrom=source_chrom,source_distance=source_dist,same_chromosome_as_source=same_source,source_strand=source_strand,candidate_source_strand_match=subev['strand_match'],synteny_overlap_fraction=f'{synteny_fraction:.4f}',parent_candidate_exon_ratio=f"{subev['exon_ratio']:.4f}",intron_loss_count=subev['intron_loss_count'],processed_evidence_score=f"{subev['processed_score']:.2f}",processed_confidence=subev['processed_confidence'],duplicated_evidence_score=f"{subev['duplicated_score']:.2f}",duplicated_confidence=subev['duplicated_confidence'],alignment_source=hit.get('alignment_source','') if hit else '',validator_status=('validated_miniprot' if hit and hit.get('alignment_source')=='miniprot' else 'supported_mmseqs' if hit else 'no_family_parent_alignment'),disruption_evidence_level=level,pseudogene_subtype=subtype,selected_parent=selected,initial_parent=row['query_id'],parent_changed=int(bool(hit and hit['parent_changed'])),parent_exon_count=parent_exons,validator_evidence=evidence)
   if hit: values.update(candidate_exon_count=subev['candidate_exons'],validated_query_coverage=f'{hit["query_coverage"]:.4f}',validated_identity=f'{hit["identity"]:.4f}',validated_genomic_start=hit.get('validated_genomic_start',''),validated_genomic_end=hit.get('validated_genomic_end',''),validated_match_fraction=f'{hit.get("match_fraction",hit["identity"]):.4f}',query_left_unaligned_fraction=f'{hit.get("query_left_unaligned_fraction",0):.4f}',query_right_unaligned_fraction=f'{hit.get("query_right_unaligned_fraction",0):.4f}',validated_aligned_length=hit['aligned_length'],validated_mapq=hit['mapq'],validated_frameshift_count=hit['frameshift_count'],validated_stop_count=hit['stop_count'],validated_alignment_score=hit['alignment_score'],validated_strand=hit['strand'],polyA_run=hit['polyA_run'],polyA_fraction=f'{hit["polyA_fraction"]:.4f}',tsd_length=hit.get('tsd_length',0),tsd_sequence=hit.get('tsd_sequence',''))
   values['structural_tags_complete']=hit.get('structural_tags_complete',0) if hit else ''
   w.writerow({**row,**values})
 return dict(counts),dict(subtype_counts)

def main(argv=None):
 ap=argparse.ArgumentParser(description=__doc__)
 ap.add_argument('--candidates',required=True,type=Path); ap.add_argument('--repeat-bed',type=Path); ap.add_argument('--synteny-bed',type=Path); ap.add_argument('--genome',required=True,type=Path); ap.add_argument('--proteins',required=True,type=Path); ap.add_argument('--gff',type=Path); ap.add_argument('--workdir',required=True,type=Path); ap.add_argument('--out-tsv',required=True,type=Path); ap.add_argument('--summary-json',required=True,type=Path)
 ap.add_argument('--fallback-m8',type=Path); ap.add_argument('--catalog-summary',type=Path); ap.add_argument('--threads',type=int,default=16); ap.add_argument('--flank',type=int,default=500); ap.add_argument('--max-flank',type=int,default=5000); ap.add_argument('--max-alternative-parents',type=int,default=10); ap.add_argument('--miniprot-outn',type=int,default=500); ap.add_argument('--miniprot-outs',type=float,default=0.1); ap.add_argument('--reuse-paf',type=Path)
 args=ap.parse_args(argv); args.workdir.mkdir(parents=True,exist_ok=True); candidates=load_candidates(args.candidates); parent_exons=load_parent_exon_counts(args.gff); _cds,gene_loci=load_gff_indexes(args.gff)
 winfa,qfa,mapping,prep=prepare_inputs(candidates,args.genome,args.proteins,args.workdir,args.flank,args.max_flank,args.max_alternative_parents); paf=args.reuse_paf or args.workdir/'candidate_family_validation.paf'
 if not args.reuse_paf and not mapping:
  paf.write_text('')
 elif not args.reuse_paf:
  tmp=paf.with_name(paf.name+f'.tmp.{__import__("os").getpid()}'); cmd=['miniprot','--outn',str(args.miniprot_outn),'--outs',str(args.miniprot_outs),'-t',str(args.threads),str(winfa),str(qfa)]
  with tmp.open('w') as out,(args.workdir/'miniprot.stderr.log').open('w') as err: process=subprocess.run(cmd,stdout=out,stderr=err,text=True)
  if process.returncode!=0: raise SystemExit(f'miniprot failed: exit {process.returncode}')
  tmp.replace(paf)
 hits,paf_summary=read_best_family_paf(paf,mapping); fallback_summary={}
 fallback_m8=args.fallback_m8
 if not fallback_m8 and args.catalog_summary:
  try: fallback_m8=Path(json.loads(args.catalog_summary.read_text()).get('mmseqs_fast',{}).get('m8',''))
  except Exception: fallback_m8=None
 if fallback_m8 and fallback_m8.is_file():
  fallback,fallback_summary=read_mmseqs_family_fallback(fallback_m8,mapping,hits);hits.update(fallback)
 levels,subtypes=write_validated(candidates,hits,parent_exons,args.out_tsv,load_repeat_index(args.repeat_bed),gene_loci,load_interval_bed(args.synteny_bed))
 summary={'input_candidates':len(candidates),'prep':prep,'paf':str(paf),**paf_summary,**fallback_summary,'total_supported_candidates':len(hits),'evidence_levels':levels,'pseudogene_subtypes':subtypes,'parameters':{'threads':args.threads,'base_flank':args.flank,'max_flank':args.max_flank,'max_alternative_parents':args.max_alternative_parents,'miniprot_outn':args.miniprot_outn,'miniprot_outs':args.miniprot_outs}}
 with atomic_text_writer(args.summary_json) as out: out.write(json.dumps(summary,indent=2,sort_keys=True)+'\n')
 print(json.dumps(summary,indent=2,sort_keys=True))
if __name__=='__main__': main()
