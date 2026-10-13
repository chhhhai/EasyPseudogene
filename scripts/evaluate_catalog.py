#!/usr/bin/env python3
"""Evaluate candidate loci with auditable one-to-one prediction-truth matching."""
import argparse,csv,json,sys
from collections import defaultdict,deque
from pathlib import Path
LIB=Path(__file__).resolve().parents[1]/'lib'
sys.path.insert(0,str(LIB))
from catalog_utils import atomic_text_writer
MATCHING_POLICY='one_to_one_maximum_cardinality'
def load_truth(path):
 truth=[];by=defaultdict(list)
 with Path(path).open() as fh:
  reader=csv.DictReader(fh,delimiter='\t')
  for i,row in enumerate(reader):
   chrom=row.get('chr') or row.get('chrom');start=int(row['start'])-1;end=int(row['end']);identifier=row.get('id') or row.get('name') or row.get('gene_id') or f'TRUTH_{i+1}';item=(chrom,start,end,identifier);truth.append(item);by[chrom].append(item)
 for chrom in by:by[chrom].sort(key=lambda x:(x[1],x[2],str(x[3])))
 return truth,by
def load_prediction_records(path):
 records=[]
 with Path(path).open() as fh:
  for i,row in enumerate(csv.DictReader(fh,delimiter='\t')):
   records.append({'prediction_id':row.get('candidate_id') or row.get('id') or row.get('name') or f'PRED_{i+1}','chrom':row.get('chr') or row.get('chrom'),'start0':int(row['start'])-1,'end':int(row['end']),'row':row})
 return records
def load_predictions(path):return [(r['chrom'],r['start0'],r['end']) for r in load_prediction_records(path)]
def maximum_cardinality_matching(adjacency,n_truth):
 n_pred=len(adjacency);pair_pred=[-1]*n_pred;pair_truth=[-1]*n_truth;distance=[0]*n_pred
 def bfs():
  queue=deque();found=False
  for pred in range(n_pred):
   if pair_pred[pred]<0:distance[pred]=0;queue.append(pred)
   else:distance[pred]=-1
  while queue:
   pred=queue.popleft()
   for truth in adjacency[pred]:
    mate=pair_truth[truth]
    if mate<0:found=True
    elif distance[mate]<0:distance[mate]=distance[pred]+1;queue.append(mate)
  return found
 def dfs(pred):
  for truth in adjacency[pred]:
   mate=pair_truth[truth]
   if mate<0 or (distance[mate]==distance[pred]+1 and dfs(mate)):
    pair_pred[pred]=truth;pair_truth[truth]=pred;return True
  distance[pred]=-1;return False
 while bfs():
  for pred in range(n_pred):
   if pair_pred[pred]<0:dfs(pred)
 return [(pred,truth) for pred,truth in enumerate(pair_pred) if truth>=0]
def match_loci(predictions,truth,truth_by_chrom,reciprocal_cutoff=None,prediction_ids=None):
 truth_index={item:i for i,item in enumerate(truth)};adjacency=[];edges={}
 for pred_index,(chrom,start,end) in enumerate(predictions):
  eligible=[]
  for truth_item in truth_by_chrom.get(chrom,[]):
   if len(truth_item)==4:_truth_chrom,truth_start,truth_end,truth_id=truth_item
   else:truth_start,truth_end,truth_id=truth_item;_truth_chrom=chrom
   overlap=max(0,min(end,truth_end)-max(start,truth_start))
   if overlap<=0:continue
   pred_fraction=overlap/max(1,end-start);truth_fraction=overlap/max(1,truth_end-truth_start)
   if reciprocal_cutoff is not None and (pred_fraction<reciprocal_cutoff or truth_fraction<reciprocal_cutoff):continue
   index=truth_index[(chrom,truth_start,truth_end,truth_id)];quality=min(pred_fraction,truth_fraction);eligible.append((quality,overlap,index));edges[(pred_index,index)]={'prediction_index':pred_index,'prediction_id':prediction_ids[pred_index] if prediction_ids else f'PRED_{pred_index+1}','prediction_chrom':chrom,'prediction_start':start+1,'prediction_end':end,'truth_index':index,'truth_id':truth_id,'truth_chrom':chrom,'truth_start':truth_start+1,'truth_end':truth_end,'overlap_bp':overlap,'prediction_overlap_fraction':pred_fraction,'truth_overlap_fraction':truth_fraction,'reciprocal_min_fraction':quality}
  eligible.sort(key=lambda x:(-x[0],-x[1],x[2]));adjacency.append([index for _quality,_overlap,index in eligible])
 matching=maximum_cardinality_matching(adjacency,len(truth));pairs=[edges[pair] for pair in matching];pairs.sort(key=lambda x:(x['prediction_chrom'],x['prediction_start'],x['prediction_end'],str(x['prediction_id'])));tp=len(pairs);precision=tp/max(1,len(predictions));recall=tp/max(1,len(truth));f1=2*precision*recall/max(1e-12,precision+recall);metrics={'predictions':len(predictions),'truth_total':len(truth),'matched_pairs':tp,'true_positives':tp,'precision':precision,'recall':recall,'f1':f1,'matching_policy':MATCHING_POLICY};return metrics,pairs
def evaluate(predictions,truth,truth_by_chrom,reciprocal_cutoff=None):return match_loci(predictions,truth,truth_by_chrom,reciprocal_cutoff)[0]
def write_pairs(path,pairs,threshold):
 fields=['threshold','matching_policy','prediction_index','prediction_id','prediction_chrom','prediction_start','prediction_end','truth_index','truth_id','truth_chrom','truth_start','truth_end','overlap_bp','prediction_overlap_fraction','truth_overlap_fraction','reciprocal_min_fraction']
 with atomic_text_writer(path) as out:
  writer=csv.DictWriter(out,fieldnames=fields,delimiter='\t',lineterminator='\n');writer.writeheader()
  for pair in pairs:writer.writerow({'threshold':threshold,'matching_policy':MATCHING_POLICY,**pair})
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--predictions',required=True,type=Path);ap.add_argument('--truth',required=True,type=Path);ap.add_argument('--out-tsv',required=True,type=Path);ap.add_argument('--out-json',type=Path);ap.add_argument('--matched-pairs-prefix',type=Path);args=ap.parse_args();truth,by=load_truth(args.truth);records=load_prediction_records(args.predictions);pred=[(r['chrom'],r['start0'],r['end']) for r in records];ids=[r['prediction_id'] for r in records];results={}
 for label,cutoff,suffix in [('>=1bp',None,'1bp'),('>=50% reciprocal',.5,'reciprocal50'),('>=80% reciprocal',.8,'reciprocal80')]:
  metrics,pairs=match_loci(pred,truth,by,cutoff,ids);results[label]=metrics
  if args.matched_pairs_prefix:write_pairs(Path(f'{args.matched_pairs_prefix}.{suffix}.tsv'),pairs,label)
 fields=['threshold','predictions','truth_total','matched_pairs','true_positives','precision','recall','f1','matching_policy'];args.out_tsv.parent.mkdir(parents=True,exist_ok=True)
 with atomic_text_writer(args.out_tsv) as out:
  writer=csv.DictWriter(out,fieldnames=fields,delimiter='\t',lineterminator='\n');writer.writeheader()
  for label,metrics in results.items():writer.writerow({'threshold':label,**metrics})
 if args.out_json:
  with atomic_text_writer(args.out_json) as out:out.write(json.dumps(results,indent=2)+'\n')
 print(args.out_tsv.read_text())
if __name__=='__main__':main()
