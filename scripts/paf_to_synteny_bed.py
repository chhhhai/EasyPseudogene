#!/usr/bin/env python3
"""Convert long, unique whole-genome PAF alignments into merged target synteny intervals."""
import argparse,json
from pathlib import Path

def merge(intervals,gap):
 out=[]
 for chrom,start,end in sorted(intervals):
  if out and out[-1][0]==chrom and start<=out[-1][2]+gap:out[-1]=(chrom,out[-1][1],max(out[-1][2],end))
  else:out.append((chrom,start,end))
 return out
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--paf',required=True,type=Path);ap.add_argument('--out-bed',required=True,type=Path);ap.add_argument('--summary-json',required=True,type=Path);ap.add_argument('--min-alignment',type=int,default=5000);ap.add_argument('--min-identity',type=float,default=.70);ap.add_argument('--min-mapq',type=int,default=20);ap.add_argument('--merge-gap',type=int,default=5000);args=ap.parse_args();kept=[];total=0
 with args.paf.open() as h:
  for line in h:
   if not line.strip():continue
   f=line.rstrip('\n').split('\t');total+=1
   try:target=f[5];start=int(f[7]);end=int(f[8]);matches=int(f[9]);block=int(f[10]);mapq=int(f[11]);identity=matches/max(1,block)
   except (ValueError,IndexError):continue
   if block>=args.min_alignment and identity>=args.min_identity and mapq>=args.min_mapq:kept.append((target,min(start,end),max(start,end)))
 merged=merge(kept,args.merge_gap);args.out_bed.parent.mkdir(parents=True,exist_ok=True);args.out_bed.write_text(''.join(f'{c}\t{s}\t{e}\n' for c,s,e in merged));summary={'paf_records':total,'kept_alignments':len(kept),'merged_intervals':len(merged),'syntenic_bp':sum(e-s for _,s,e in merged),'thresholds':{'min_alignment':args.min_alignment,'min_identity':args.min_identity,'min_mapq':args.min_mapq,'merge_gap':args.merge_gap}};args.summary_json.write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
