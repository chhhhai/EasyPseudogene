#!/usr/bin/env python3
"""Prepare MCScanX GFF/BLAST anchors from cross-species orthogroups and genome annotations."""
import argparse,csv,sys
from collections import defaultdict
from pathlib import Path
LIB=Path(__file__).resolve().parents[1]/'lib';sys.path.insert(0,str(LIB))
from self_catalog_blast_fast import load_gff_indexes

def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--orthogroups',required=True,type=Path);ap.add_argument('--species-a',required=True);ap.add_argument('--gff-a',required=True,type=Path);ap.add_argument('--species-b',required=True);ap.add_argument('--gff-b',required=True,type=Path);ap.add_argument('--prefix',required=True,type=Path);args=ap.parse_args();_,a_loci=load_gff_indexes(args.gff_a);_,b_loci=load_gff_indexes(args.gff_b);groups=defaultdict(lambda:defaultdict(list))
 with args.orthogroups.open() as h:
  for r in csv.DictReader(h,delimiter='\t'):
   if r['species'] in {args.species_a,args.species_b}:groups[r['orthogroup_id']][r['species']].append(r['protein_id'])
 selected={};records=[]
 for species,loci_index,label in [(args.species_a,a_loci,'A'),(args.species_b,b_loci,'B')]:
  used=set()
  for group,members in groups.items():
   for protein in members.get(species,[]):
    loci=loci_index.get(protein,[])
    if not loci:continue
    chrom,start,end,_strand=loci[0];key=(chrom,start,end)
    if key in used:continue
    used.add(key);gene=f'{label}|{protein}';selected[(species,protein)]=gene;records.append((f'{label}_{chrom}',gene,start+1,end));break
 args.prefix.parent.mkdir(parents=True,exist_ok=True)
 with Path(str(args.prefix)+'.gff').open('w') as h:
  for row in sorted(records):h.write('\t'.join(map(str,row))+'\n')
 pairs=0
 with Path(str(args.prefix)+'.blast').open('w') as h:
  for group,members in groups.items():
   aa=[selected[(args.species_a,p)] for p in members.get(args.species_a,[]) if (args.species_a,p) in selected];bb=[selected[(args.species_b,p)] for p in members.get(args.species_b,[]) if (args.species_b,p) in selected]
   for a in aa:
    for b in bb:h.write(f'{a}\t{b}\t1e-50\t200\n{b}\t{a}\t1e-50\t200\n');pairs+=1
 print(f'{len(records)} positioned genes; {pairs} cross-species orthogroup anchors')
if __name__=='__main__':main()
