#!/usr/bin/env python3
"""Create deterministic evidence-stratified manual-review candidates with IGV links."""
import argparse,csv,random
from pathlib import Path

def read_pairs(path):
 if not path:return {}
 with Path(path).open() as handle:return {row['prediction_id']:row for row in csv.DictReader(handle,delimiter='\t')}
def choose(rows,n,rng):
 rows=sorted(rows,key=lambda r:(r.get('chr',''),int(r.get('start',0)),int(r.get('end',0)),r.get('candidate_id','')))
 return rows if len(rows)<=n else sorted(rng.sample(rows,n),key=lambda r:(r.get('chr',''),int(r.get('start',0))))
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--validated-tsv',required=True,type=Path);ap.add_argument('--matched-1bp',required=True,type=Path);ap.add_argument('--matched-50',required=True,type=Path);ap.add_argument('--out-tsv',required=True,type=Path);ap.add_argument('--per-stratum',type=int,default=20);ap.add_argument('--seed',type=int,default=20260715);ap.add_argument('--igv-template',default='http://localhost:60151/goto?locus={chrom}:{start}-{end}');args=ap.parse_args();rng=random.Random(args.seed)
 with args.validated_tsv.open() as handle:
  reader=csv.DictReader(handle,delimiter='\t');rows=list(reader);base_fields=reader.fieldnames or []
 match1=read_pairs(args.matched_1bp);match50=read_pairs(args.matched_50)
 def matched(row):return row.get('candidate_id') in match1
 strata={
  'matched_strict_core':[r for r in rows if matched(r) and r.get('disruption_evidence_level')=='high'],
  'unmatched_strict_core':[r for r in rows if not matched(r) and r.get('disruption_evidence_level')=='high'],
  'matched_fallback_only':[r for r in rows if matched(r) and r.get('alignment_source')=='mmseqs_fallback'],
  'unmatched_fallback_only':[r for r in rows if not matched(r) and r.get('alignment_source')=='mmseqs_fallback'],
  'boundary_mismatch_1bp_not_50':[r for r in rows if r.get('candidate_id') in match1 and r.get('candidate_id') not in match50],
  'processed_like':[r for r in rows if r.get('pseudogene_subtype')=='processed_like'],
  'duplicated_like':[r for r in rows if r.get('pseudogene_subtype')=='duplicated_like'],
 }
 selected=[]
 for stratum,candidates in strata.items():
  for row in choose(candidates,args.per_stratum,rng):selected.append((stratum,row))
 extra=['review_stratum','matched_1bp','matched_50','truth_id_1bp','truth_start_1bp','truth_end_1bp','overlap_bp_1bp','prediction_overlap_fraction_1bp','truth_overlap_fraction_1bp','igv_url'];fields=extra+base_fields;args.out_tsv.parent.mkdir(parents=True,exist_ok=True)
 with args.out_tsv.open('w',newline='') as handle:
  writer=csv.DictWriter(handle,fieldnames=fields,delimiter='\t',lineterminator='\n');writer.writeheader()
  for stratum,row in selected:
   candidate=row.get('candidate_id');pair=match1.get(candidate,{});start=row.get('start','');end=row.get('end','');chrom=row.get('chr') or row.get('chrom','');extra_row={'review_stratum':stratum,'matched_1bp':int(candidate in match1),'matched_50':int(candidate in match50),'truth_id_1bp':pair.get('truth_id',''),'truth_start_1bp':pair.get('truth_start',''),'truth_end_1bp':pair.get('truth_end',''),'overlap_bp_1bp':pair.get('overlap_bp',''),'prediction_overlap_fraction_1bp':pair.get('prediction_overlap_fraction',''),'truth_overlap_fraction_1bp':pair.get('truth_overlap_fraction',''),'igv_url':args.igv_template.format(chrom=chrom,start=start,end=end,candidate_id=candidate)};writer.writerow({**extra_row,**row})
 print(f'wrote {len(selected)} review rows across {len(strata)} strata to {args.out_tsv}')
if __name__=='__main__':main()
