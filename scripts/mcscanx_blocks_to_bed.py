#!/usr/bin/env python3
"""Convert MCScanX collinear blocks to merged BED intervals on one species prefix."""
import argparse,json
from pathlib import Path

def merge(rows,gap=0):
 out=[]
 for chrom,start,end in sorted(rows):
  if out and out[-1][0]==chrom and start<=out[-1][2]+gap:out[-1]=(chrom,out[-1][1],max(out[-1][2],end))
  else:out.append((chrom,start,end))
 return out
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--gff',required=True,type=Path);ap.add_argument('--collinearity',required=True,type=Path);ap.add_argument('--species-prefix',default='A|');ap.add_argument('--out-bed',required=True,type=Path);ap.add_argument('--summary-json',required=True,type=Path);ap.add_argument('--padding',type=int,default=10000);args=ap.parse_args();genes={}
 with args.gff.open() as h:
  for line in h:
   chrom,gene,start,end=line.rstrip('\n').split('\t');genes[gene]=(chrom[2:] if chrom.startswith(args.species_prefix[0]+'_') else chrom,int(start)-1,int(end))
 blocks=[];current=[];headers=0
 def flush():
  if not current:return
  loci=[genes[g] for g in current if g in genes]
  if not loci:return
  by={}
  for chrom,start,end in loci:by.setdefault(chrom,[]).append((start,end))
  for chrom,v in by.items():blocks.append((chrom,max(0,min(x[0] for x in v)-args.padding),max(x[1] for x in v)+args.padding))
 with args.collinearity.open() as h:
  for line in h:
   if line.startswith('## Alignment'):
    flush();current=[];headers+=1
   elif line.startswith(' ') and '\t' in line:
    f=line.rstrip('\n').split('\t')
    if len(f)>=3:
     for gene in f[1:3]:
      if gene.startswith(args.species_prefix):current.append(gene)
 flush();merged=merge(blocks);args.out_bed.write_text(''.join(f'{c}\t{s}\t{e}\n' for c,s,e in merged));summary={'mcscanx_blocks':headers,'species_block_intervals':len(blocks),'merged_intervals':len(merged),'syntenic_bp':sum(e-s for _,s,e in merged),'padding':args.padding};args.summary_json.write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
