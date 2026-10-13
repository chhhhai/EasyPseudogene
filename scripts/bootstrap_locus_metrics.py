#!/usr/bin/env python3
"""Chromosome-bootstrap confidence intervals for one-to-one locus metrics."""
import argparse,importlib.util,json,random,statistics
from pathlib import Path
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('evaluate_catalog',HERE/'evaluate_catalog.py');ev=importlib.util.module_from_spec(spec);spec.loader.exec_module(ev)
def load_predictions(path,fmt):
 if fmt=='tsv1':return ev.load_predictions(path)
 rows=[]
 with Path(path).open() as handle:
  for line in handle:
   if not line.strip() or line.startswith('#'):continue
   f=line.rstrip('\n').split('\t')
   try:rows.append((f[0],int(f[1]),int(f[2])))
   except (ValueError,IndexError):continue
 return rows
def quantile(values,q):
 values=sorted(values);position=(len(values)-1)*q;lo=int(position);hi=min(len(values)-1,lo+1);fraction=position-lo;return values[lo]*(1-fraction)+values[hi]*fraction
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--predictions',required=True,type=Path);ap.add_argument('--truth',required=True,type=Path);ap.add_argument('--pred-format',choices=['tsv1','bed0'],default='tsv1');ap.add_argument('--chromosomes',default='1,2,3,4,5');ap.add_argument('--replicates',type=int,default=2000);ap.add_argument('--seed',type=int,default=20260715);ap.add_argument('--out-json',required=True,type=Path);args=ap.parse_args();pred=load_predictions(args.predictions,args.pred_format);truth,truth_by=ev.load_truth(args.truth);chroms=args.chromosomes.split(',');criteria={'1bp':None,'50':.5,'80':.8};per={}
 for chrom in chroms:
  cp=[x for x in pred if x[0]==chrom];ct=[x for x in truth if x[0]==chrom];cb={chrom:ct};per[chrom]={}
  for label,cutoff in criteria.items():per[chrom][label]=ev.evaluate(cp,ct,cb,cutoff)
 rng=random.Random(args.seed);samples={label:{'precision':[],'recall':[],'f1':[]} for label in criteria}
 for _ in range(args.replicates):
  chosen=[rng.choice(chroms) for _ in chroms]
  for label in criteria:
   np=sum(per[c][label]['predictions'] for c in chosen);nt=sum(per[c][label]['truth_total'] for c in chosen);tp=sum(per[c][label]['true_positives'] for c in chosen);precision=tp/max(1,np);recall=tp/max(1,nt);f1=2*precision*recall/max(1e-12,precision+recall)
   for key,value in [('precision',precision),('recall',recall),('f1',f1)]:samples[label][key].append(value)
 result={'replicates':args.replicates,'seed':args.seed,'chromosomes':chroms,'matching':'one_to_one_maximum_cardinality','per_chromosome':per,'bootstrap_ci95':{label:{key:{'median':statistics.median(values),'lower':quantile(values,.025),'upper':quantile(values,.975)} for key,values in metrics.items()} for label,metrics in samples.items()}}
 args.out_json.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result['bootstrap_ci95'],indent=2))
if __name__=='__main__':main()
