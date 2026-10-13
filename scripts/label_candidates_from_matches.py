#!/usr/bin/env python3
"""Append one-to-one reference-concordance labels and family-aware folds to candidate TSV."""
import argparse,csv,hashlib
from pathlib import Path

def stable_fold(species,family,folds):
    if folds < 2: raise ValueError('At least two folds are required')
    key=str(family).encode()
    return int(hashlib.sha256(key).hexdigest()[:16],16)%folds

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--candidates',required=True,type=Path)
    ap.add_argument('--matched-pairs',required=True,type=Path)
    ap.add_argument('--species',required=True)
    ap.add_argument('--out',required=True,type=Path)
    ap.add_argument('--folds',type=int,default=5);ap.add_argument('--family-map',type=Path,help='TSV from build_cross_species_orthogroups.py')
    args=ap.parse_args()
    family_map={}
    if args.family_map:
        with args.family_map.open() as handle:
            for item in csv.DictReader(handle,delimiter='\t'):
                if item['species']==args.species:family_map[item['protein_id']]=item['orthogroup_id']
    with args.matched_pairs.open() as handle:
        matched={row['prediction_id'] for row in csv.DictReader(handle,delimiter='\t')}
    with args.candidates.open() as handle:
        reader=csv.DictReader(handle,delimiter='\t');rows=list(reader);fields=reader.fieldnames or []
    args.out.parent.mkdir(parents=True,exist_ok=True)
    with args.out.open('w',newline='') as out:
        writer=csv.DictWriter(out,fieldnames=list(dict.fromkeys(fields+['benchmark_species','family_id','family_fold','truth_overlap'])),delimiter='\t',lineterminator='\n');writer.writeheader()
        for row in rows:
            candidate=row.get('candidate_id') or row.get('id');parent=row.get('selected_parent') or row.get('query_id') or candidate;family=family_map.get(parent,parent)
            writer.writerow({**row,'benchmark_species':args.species,'family_id':family,'family_fold':stable_fold(args.species,family,args.folds),'truth_overlap':int(candidate in matched)})
    print(f'{len(rows)} rows; {sum((r.get("candidate_id") or r.get("id")) in matched for r in rows)} positives')
if __name__=='__main__':main()
