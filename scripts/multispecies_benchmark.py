#!/usr/bin/env python3
"""Manifest validation, family-aware folds, and leave-one-species-out logistic evaluation."""
import argparse,csv,hashlib,json,math
from pathlib import Path

def stable_fold(species,family,folds):return int(hashlib.sha256(str(family).encode()).hexdigest()[:16],16)%folds
def assign_folds(input_tsv,out_tsv,species,family_column,folds):
 with Path(input_tsv).open() as handle:reader=csv.DictReader(handle,delimiter='\t');rows=list(reader);fields=reader.fieldnames or []
 with Path(out_tsv).open('w',newline='') as out:
  writer=csv.DictWriter(out,fieldnames=list(dict.fromkeys(fields+['benchmark_species','family_id','family_fold'])),delimiter='\t',lineterminator='\n');writer.writeheader()
  for row in rows:
   family=row.get(family_column) or row.get('selected_parent') or row.get('query_id') or row.get('candidate_id');writer.writerow({**row,'benchmark_species':species,'family_id':family,'family_fold':stable_fold(species,family,folds)})
 return len(rows)
def sigmoid(x):return 1/(1+math.exp(-max(-35,min(35,x))))
def fit_logistic(x,y,epochs=500,rate=.05,l2=.001):
 if not x:return []
 try:
  import numpy as np
  matrix=np.asarray(x,dtype=float);target=np.asarray(y,dtype=float);design=np.column_stack([np.ones(len(matrix)),matrix]);weights=np.zeros(design.shape[1])
  for _ in range(epochs):
   probability=1/(1+np.exp(-np.clip(design@weights,-35,35)));gradient=design.T@(probability-target)/len(target);gradient[1:]+=l2*weights[1:];weights-=rate*gradient
  return weights.tolist()
 except ImportError:
  weights=[0.0]*(len(x[0])+1);n=len(x)
  for _ in range(epochs):
   grad=[0.0]*len(weights)
   for row,target in zip(x,y):
    error=sigmoid(weights[0]+sum(w*v for w,v in zip(weights[1:],row)))-target;grad[0]+=error
    for j,value in enumerate(row,1):grad[j]+=error*value
   for j in range(len(weights)):weights[j]-=rate*(grad[j]/n+(l2*weights[j] if j else 0))
  return weights
def auc(y,score):
 pairs=sorted(zip(score,y),key=lambda x:x[0]);positives=sum(y);negatives=len(y)-positives
 if not positives or not negatives:return float('nan')
 rank_sum=0.0;i=0
 while i<len(pairs):
  j=i+1
  while j<len(pairs) and pairs[j][0]==pairs[i][0]:j+=1
  average_rank=(i+1+j)/2;rank_sum+=average_rank*sum(label for _score,label in pairs[i:j]);i=j
 return (rank_sum-positives*(positives+1)/2)/(positives*negatives)
def average_precision(y,score):
 ordered=sorted(zip(score,y),reverse=True);total=sum(y)
 if not total:return float('nan')
 tp=0;value=0.0
 for rank,(_score,label) in enumerate(ordered,1):
  if label:tp+=1;value+=tp/rank
 return value/total
def calibration_bins(y,score,bins=10):
 output=[]
 for index in range(bins):
  lo=index/bins;hi=(index+1)/bins;values=[(label,pred) for label,pred in zip(y,score) if lo<=pred<(hi if index<bins-1 else hi+1e-12)]
  output.append({'lower':lo,'upper':hi,'count':len(values),'mean_probability':sum(v[1] for v in values)/len(values) if values else None,'observed_fraction':sum(v[0] for v in values)/len(values) if values else None})
 return output
def metrics(y,score,threshold=.5):
 pred=[s>=threshold for s in score];tp=sum(a and b for a,b in zip(y,pred));fp=sum((not a) and b for a,b in zip(y,pred));fn=sum(a and not b for a,b in zip(y,pred));tn=sum((not a) and (not b) for a,b in zip(y,pred));p=tp/max(1,tp+fp);r=tp/max(1,tp+fn);return {'tp':tp,'fp':fp,'tn':tn,'fn':fn,'prevalence':sum(y)/max(1,len(y)),'precision':p,'recall':r,'specificity':tn/max(1,tn+fp),'f1':2*p*r/max(1e-12,p+r),'roc_auc':auc(y,score),'pr_auc':average_precision(y,score),'brier_score':sum((a-b)**2 for a,b in zip(y,score))/max(1,len(y)),'calibration':calibration_bins(y,score)}
def best_f1_threshold(y,score):
 ordered=sorted(zip(score,y),reverse=True);total_positive=sum(y);tp=fp=0;best=(0.0,.5)
 for index,(threshold,label) in enumerate(ordered):
  if label:tp+=1
  else:fp+=1
  if index+1<len(ordered) and ordered[index+1][0]==threshold:continue
  fn=total_positive-tp;precision=tp/max(1,tp+fp);recall=tp/max(1,tp+fn);f1=2*precision*recall/max(1e-12,precision+recall)
  if f1>best[0]:best=(f1,threshold)
 return best
def load_labeled(path,features,label_column):
 rows=[]
 with Path(path).open() as handle:
  for row in csv.DictReader(handle,delimiter='\t'):
   try:rows.append(([float(row.get(f) or 0) for f in features],int(float(row[label_column])>0),row))
   except (ValueError,KeyError):continue
 return rows
def loso(manifest_path,features,label_column,out_json,predictions_dir=None):
 manifest=json.loads(Path(manifest_path).read_text());datasets={item['species']:load_labeled(item['labeled_tsv'],features,label_column) for item in manifest['datasets']};results={}
 # Input fold columns may come from legacy species-dependent hashing. Recompute.
 for species,records in datasets.items():
  for _,_,row in records:
   family=row.get('family_id') or row.get('selected_parent') or row.get('query_id')
   if not family:raise ValueError('Missing family_id/query_id for family-aware evaluation')
   row['family_id']=family;row['family_fold']=str(stable_fold(species,family,5))
 def fit_records(records):
  means=[sum(r[0][j] for r in records)/max(1,len(records)) for j in range(len(features))];sds=[math.sqrt(sum((r[0][j]-means[j])**2 for r in records)/max(1,len(records))) or 1 for j in range(len(features))];transform=lambda row:[(row[j]-means[j])/sds[j] for j in range(len(features))];x=[transform(r[0]) for r in records];y=[r[1] for r in records];return means,sds,transform,fit_logistic(x,y)
 def predict(records,transform,weights):return [sigmoid(weights[0]+sum(w*v for w,v in zip(weights[1:],transform(r[0])))) for r in records]
 for held,test in datasets.items():
  train=[record for species,records in datasets.items() if species!=held for record in records];fold_values=sorted({int(record[2].get('family_fold') or stable_fold(record[2].get('benchmark_species',''),record[2].get('family_id',''),5)) for record in train});oof_y=[];oof_score=[]
  for fold in fold_values:
   fit=[r for r in train if int(r[2].get('family_fold') or stable_fold(r[2].get('benchmark_species',''),r[2].get('family_id',''),5))!=fold];valid=[r for r in train if int(r[2].get('family_fold') or stable_fold(r[2].get('benchmark_species',''),r[2].get('family_id',''),5))==fold]
   if not fit or not valid:continue
   _mean,_sd,transform,weights=fit_records(fit);oof_y.extend(r[1] for r in valid);oof_score.extend(predict(valid,transform,weights))
  train_best_f1,threshold=best_f1_threshold(oof_y,oof_score) if oof_y else (0.0,.5);means,sds,transform,weights=fit_records(train);scores=predict(test,transform,weights);
  if predictions_dir:
   directory=Path(predictions_dir);directory.mkdir(parents=True,exist_ok=True);path=directory/f'{held}.held_out_predictions.tsv'
   with path.open('w',newline='') as handle:
    writer=csv.writer(handle,delimiter='\t',lineterminator='\n');writer.writerow(['candidate_id','family_id','family_fold','label','score','locked_threshold','predicted_positive'])
    for record,score in zip(test,scores):writer.writerow([record[2].get('candidate_id',''),record[2].get('family_id',''),record[2].get('family_fold',''),record[1],score,threshold,int(score>=threshold)])
  results[held]={'train_species':[species for species in datasets if species!=held],'train_rows':len(train),'test_rows':len(test),'family_folds':fold_values,'threshold_selection':'family_fold_out_of_fold_training_species_only','locked_threshold':threshold,'train_oof_best_f1':train_best_f1,'feature_means':dict(zip(features,means)),'feature_sds':dict(zip(features,sds)),'weights':{'intercept':weights[0],**dict(zip(features,weights[1:]))},**metrics([r[1] for r in test],scores,threshold)}
 output={'mode':'leave_one_species_out','features':features,'label_column':label_column,'fold_policy':'sha256_family_only_recomputed_5fold','results':results};Path(out_json).write_text(json.dumps(output,indent=2)+'\n');return output

def validate_manifest(path):
 data=json.loads(Path(path).read_text());issues=[];seen=set()
 for item in data.get('datasets',[]):
  species=item.get('species');
  if not species or species in seen:issues.append(f'invalid or duplicate species: {species}')
  seen.add(species)
  for key in ('genome','proteins','gff','truth'):
   if key in item and not Path(item[key]).is_file():issues.append(f'{species}: missing {key}: {item[key]}')
 return {'datasets':len(data.get('datasets',[])),'issues':issues,'valid':not issues}
def main():
 ap=argparse.ArgumentParser(description=__doc__);sub=ap.add_subparsers(dest='command',required=True)
 a=sub.add_parser('validate-manifest');a.add_argument('--manifest',required=True,type=Path)
 a=sub.add_parser('assign-folds');a.add_argument('--input',required=True,type=Path);a.add_argument('--output',required=True,type=Path);a.add_argument('--species',required=True);a.add_argument('--family-column',default='selected_parent');a.add_argument('--folds',type=int,default=5)
 a=sub.add_parser('loso');a.add_argument('--manifest',required=True,type=Path);a.add_argument('--features',required=True);a.add_argument('--label-column',default='truth_overlap');a.add_argument('--out-json',required=True,type=Path);a.add_argument('--predictions-dir',type=Path)
 args=ap.parse_args()
 if args.command=='validate-manifest':print(json.dumps(validate_manifest(args.manifest),indent=2))
 elif args.command=='assign-folds':print(assign_folds(args.input,args.output,args.species,args.family_column,args.folds))
 else:print(json.dumps(loso(args.manifest,args.features.split(','),args.label_column,args.out_json,args.predictions_dir),indent=2))
if __name__=='__main__':main()
