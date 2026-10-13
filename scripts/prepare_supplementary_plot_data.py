#!/usr/bin/env python3
"""Prepare auditable Phase-3 supplementary plotting tables."""
import argparse,csv,json
from pathlib import Path

def pr_points(path,max_points=250):
 rows=[]
 with path.open() as h:
  for r in csv.DictReader(h,delimiter='\t'):rows.append((float(r['score']),int(r['label'])))
 rows.sort(reverse=True);positives=sum(y for _,y in rows);tp=fp=0;points=[(0.0,1.0,1.0)]
 for score,y in rows:
  tp+=y;fp+=1-y;points.append((tp/max(1,positives),tp/max(1,tp+fp),score))
 step=max(1,len(points)//max_points);sample=points[::step]
 if sample[-1]!=points[-1]:sample.append(points[-1])
 return sample
def calibration(path,bins=10):
 values=[]
 with path.open() as h:
  for r in csv.DictReader(h,delimiter='\t'):values.append((float(r['score']),int(r['label'])))
 out=[]
 for i in range(bins):
  lo=i/bins;hi=(i+1)/bins;v=[x for x in values if lo<=x[0]<(hi if i<bins-1 else hi+1e-12)]
  out.append((lo,hi,len(v),sum(x for x,_ in v)/len(v) if v else '',sum(y for _,y in v)/len(v) if v else ''))
 return out
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--root',type=Path,required=True);args=ap.parse_args();root=args.root;out=root/'phase3_publication_figures/data';out.mkdir(parents=True,exist_ok=True);pred=root/'multispecies_model_data/held_out_predictions';species_map={'Arabidopsis_thaliana_TAIR10':'A. thaliana','Oryza_sativa_IRGSP1':'O. sativa','Arabidopsis_lyrata_GCF0000042552':'A. lyrata'}
 with (out/'pr_curves.tsv').open('w',newline='') as h:
  w=csv.writer(h,delimiter='\t',lineterminator='\n');w.writerow(['species','recall','precision','threshold'])
  for stem,label in species_map.items():
   for rec,prec,thr in pr_points(pred/f'{stem}.held_out_predictions.tsv'):w.writerow([label,rec,prec,thr])
 with (out/'calibration.tsv').open('w',newline='') as h:
  w=csv.writer(h,delimiter='\t',lineterminator='\n');w.writerow(['species','lower','upper','count','mean_probability','observed_fraction'])
  for stem,label in species_map.items():
   for row in calibration(pred/f'{stem}.held_out_predictions.tsv'):w.writerow([label,*row])
 boot=json.loads((root/'phase3_uncertainty/phase3_bootstrap_summary.json').read_text())
 with (out/'bootstrap_forest.tsv').open('w',newline='') as h:
  w=csv.writer(h,delimiter='\t',lineterminator='\n');w.writerow(['species','tier','median','lower','upper'])
  for species,tiers in boot.items():
   label={'TAIR10':'A. thaliana','Rice':'O. sativa','Alyrata':'A. lyrata'}[species]
   for tier,data in tiers.items():
    ci=data['1bp']['f1'];w.writerow([label,{'sensitive':'Sensitive','strict_core_disruption':'Strict core'}[tier],ci['median'],ci['lower'],ci['upper']])
 summaries={'A. thaliana':root/'benchmark_TAIR10_EasyPseudogene_vs_PGRe/runs/EasyPseudogene_phase3_release_full_20260715/output/TAIR10_phase3_release/disruption_validation/tier_summary.json','O. sativa':root/'external_benchmark_Oryza_sativa_IRGSP1_release59/runs/mmseqs_balanced_phase3_release/output/Oryza_sativa_IRGSP1_phase3_release/disruption_validation/tier_summary.json','A. lyrata':root/'external_benchmark_Arabidopsis_lyrata_GCF0000042552/runs/mmseqs_balanced_phase3/output/Arabidopsis_lyrata_GCF0000042552_phase3/disruption_validation/tier_summary.json'}
 with (out/'attrition.tsv').open('w',newline='') as h:
  w=csv.writer(h,delimiter='\t',lineterminator='\n');w.writerow(['species','stage','count','fraction_of_sensitive'])
  for species,path in summaries.items():
   c=json.loads(path.read_text())['tier_counts'];total=c['sensitive']
   for stage,key in [('Sensitive','sensitive'),('Miniprot-supported','miniprot_supported'),('Strict core','strict_core_disruption'),('Subtype candidates','subtype_candidates')]:w.writerow([species,stage,c[key],c[key]/total])
 print(out)
if __name__=='__main__':main()
