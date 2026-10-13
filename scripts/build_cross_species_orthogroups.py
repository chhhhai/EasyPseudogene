#!/usr/bin/env python3
"""Build cross-species protein families for leakage-safe benchmark folds with MMseqs2."""
import argparse,csv,subprocess,sys
from pathlib import Path
LIB=Path(__file__).resolve().parents[1]/'lib';sys.path.insert(0,str(LIB))
from catalog_utils import read_fasta,write_fasta

def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--proteome',action='append',required=True,help='SPECIES=FASTA');ap.add_argument('--workdir',required=True,type=Path);ap.add_argument('--out-map',required=True,type=Path);ap.add_argument('--threads',type=int,default=16);ap.add_argument('--min-seq-id',type=float,default=.3);ap.add_argument('--coverage',type=float,default=.5);args=ap.parse_args();args.workdir.mkdir(parents=True,exist_ok=True);combined=args.workdir/'combined_proteins.fa';records=[];metadata={}
 for item in args.proteome:
  species,path=item.split('=',1)
  for protein,desc,seq in read_fasta(Path(path)):
   prefixed=f'{species}|{protein}';records.append((prefixed,prefixed,seq));metadata[prefixed]=(species,protein)
 write_fasta(records,combined);prefix=args.workdir/'orthogroups';tmp=args.workdir/'mmseqs_tmp';cmd=['mmseqs','easy-cluster',str(combined),str(prefix),str(tmp),'--threads',str(args.threads),'--min-seq-id',str(args.min_seq_id),'-c',str(args.coverage),'--cov-mode','1'];process=subprocess.run(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
 if process.returncode:raise SystemExit(process.stderr[-3000:])
 cluster_tsv=Path(str(prefix)+'_cluster.tsv');member_to_group={}
 with cluster_tsv.open() as handle:
  for row in csv.reader(handle,delimiter='\t'):
   if len(row)>=2:member_to_group[row[1]]=row[0]
 args.out_map.parent.mkdir(parents=True,exist_ok=True)
 with args.out_map.open('w',newline='') as out:
  writer=csv.writer(out,delimiter='\t',lineterminator='\n');writer.writerow(['species','protein_id','orthogroup_id'])
  for member,(species,protein) in sorted(metadata.items()):writer.writerow([species,protein,member_to_group.get(member,member)])
 print(f'{len(records)} proteins; {len(set(member_to_group.values()))} orthogroups')
if __name__=='__main__':main()
