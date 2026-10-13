#!/usr/bin/env python3
"""Unblind two independent computational reviews and calculate agreement."""
import argparse,csv,json
from collections import Counter,defaultdict
from pathlib import Path

def load(path):
 with path.open() as h:return {r['review_id']:r for r in csv.DictReader(h,delimiter='\t')}
def kappa(a,b):
 labels=sorted(set(a)|set(b));n=len(a);observed=sum(x==y for x,y in zip(a,b))/max(1,n);ca=Counter(a);cb=Counter(b);expected=sum(ca[x]/max(1,n)*cb[x]/max(1,n) for x in labels);return (observed-expected)/max(1e-12,1-expected),observed

def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--reviewer-a',required=True,type=Path);ap.add_argument('--reviewer-b',required=True,type=Path);ap.add_argument('--key',required=True,type=Path);ap.add_argument('--out-tsv',required=True,type=Path);ap.add_argument('--out-json',required=True,type=Path);args=ap.parse_args();a=load(args.reviewer_a);b=load(args.reviewer_b);key=load(args.key);ids=sorted(set(a)&set(b)&set(key));fields=['review_id','candidate_id','review_stratum','matched_truth_id','a_candidate','b_candidate','a_subtype','b_subtype','candidate_agree','subtype_agree','a_reason','b_reason']
 with args.out_tsv.open('w',newline='') as h:
  w=csv.DictWriter(h,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader()
  for i in ids:w.writerow({'review_id':i,'candidate_id':key[i].get('candidate_id',''),'review_stratum':key[i].get('review_stratum',''),'matched_truth_id':key[i].get('matched_truth_id',''),'a_candidate':a[i]['candidate_support'],'b_candidate':b[i]['candidate_support'],'a_subtype':a[i]['subtype_support'],'b_subtype':b[i]['subtype_support'],'candidate_agree':int(a[i]['candidate_support']==b[i]['candidate_support']),'subtype_agree':int(a[i]['subtype_support']==b[i]['subtype_support']),'a_reason':a[i].get('reason',''),'b_reason':b[i].get('reason','')})
 ck,ca=kappa([a[i]['candidate_support'] for i in ids],[b[i]['candidate_support'] for i in ids]);sk,sa=kappa([a[i]['subtype_support'] for i in ids],[b[i]['subtype_support'] for i in ids]);strata=defaultdict(lambda:{'n':0,'a_supported':0,'b_supported':0,'agreement':0})
 for i in ids:
  x=strata[key[i].get('review_stratum','')];x['n']+=1;x['a_supported']+=a[i]['candidate_support']=='supported';x['b_supported']+=b[i]['candidate_support']=='supported';x['agreement']+=a[i]['candidate_support']==b[i]['candidate_support']
 for x in strata.values():
  n=x['n'];x['a_supported_fraction']=x.pop('a_supported')/n;x['b_supported_fraction']=x.pop('b_supported')/n;x['agreement_fraction']=x.pop('agreement')/n
 result={'review_type':'blinded independent computational evidence audit; not human or experimental validation','records':len(ids),'candidate_support':{'cohen_kappa':ck,'agreement_fraction':ca},'subtype_support':{'cohen_kappa':sk,'agreement_fraction':sa},'strata':dict(strata)};args.out_json.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
