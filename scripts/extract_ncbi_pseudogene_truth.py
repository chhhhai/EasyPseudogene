#!/usr/bin/env python3
"""Extract top-level protein-related pseudogene loci from an NCBI GFF3."""
import argparse
from pathlib import Path

def attrs(text):
 out={}
 for part in text.split(';'):
  if '=' in part:k,v=part.split('=',1);out[k]=v
 return out
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--gff',required=True,type=Path);ap.add_argument('--out',required=True,type=Path);ap.add_argument('--include-rna-pseudogenes',action='store_true',help='include tRNA/rRNA/ncRNA pseudogenes; excluded by default');args=ap.parse_args();rows=[]
 with args.gff.open(errors='replace') as handle:
  for line in handle:
   if not line.strip() or line.startswith('#'):continue
   f=line.rstrip('\n').split('\t')
   if len(f)<9 or f[2]!='pseudogene':continue
   a=attrs(f[8]);biotype=(a.get('gene_biotype') or a.get('gbkey') or '').lower()
   if not args.include_rna_pseudogenes and any(x in biotype for x in ('trna','rrna','snrna','snorna','ncrna')):continue
   rows.append((f[0],int(f[3]),int(f[4]),a.get('gene') or a.get('Name') or a.get('ID') or f'pseudogene_{len(rows)+1}',biotype))
 args.out.parent.mkdir(parents=True,exist_ok=True)
 with args.out.open('w') as out:
  out.write('chr\tstart\tend\tid\tbiotype\n')
  for row in rows:out.write('\t'.join(map(str,row))+'\n')
 print(len(rows))
if __name__=='__main__':main()
