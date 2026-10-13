#!/usr/bin/env python3
"""Write explicit evidence-tier catalogs from an EasyPseudogene validated TSV."""
import argparse,csv,json
from collections import Counter
from pathlib import Path
from catalog_utils import atomic_text_writer

TIER_RULES={
 'sensitive':lambda row:True,
 'miniprot_supported':lambda row:row.get('alignment_source')=='miniprot',
 'strict_core_disruption':lambda row:row.get('disruption_evidence_level')=='high' and row.get('alignment_source')=='miniprot',
 'subtype_candidates':lambda row:row.get('pseudogene_subtype') in {'processed_like','duplicated_like'},
}

def write_tiers(validated_tsv,prefix,summary_json):
 with Path(validated_tsv).open() as handle:
  reader=csv.DictReader(handle,delimiter='\t');rows=list(reader);fields=reader.fieldnames or []
 if not fields:raise ValueError('validated TSV has no header')
 outputs={};counts={};subtypes=Counter(row.get('pseudogene_subtype','') for row in rows)
 for tier,rule in TIER_RULES.items():
  selected=[row for row in rows if rule(row)];path=Path(f'{prefix}.{tier}.tsv');path.parent.mkdir(parents=True,exist_ok=True)
  with atomic_text_writer(path) as handle:
   writer=csv.DictWriter(handle,fieldnames=fields,delimiter='\t',lineterminator='\n');writer.writeheader();writer.writerows(selected)
  outputs[tier]=str(path);counts[tier]=len(selected)
 summary={'input_validated_tsv':str(validated_tsv),'input_candidates':len(rows),'tier_counts':counts,'subtype_counts':dict(sorted(subtypes.items())),'outputs':outputs,'semantics':{'sensitive':'all candidates; high-recall screening output','miniprot_supported':'structural protein-genome alignment support; not necessarily disrupted','strict_core_disruption':'miniprot-supported near-complete alignment with internal disruption evidence','subtype_candidates':'provisional processed-like or duplicated-like evidence classes'}}
 with atomic_text_writer(summary_json) as handle:handle.write(json.dumps(summary,indent=2,sort_keys=True)+'\n')
 return summary

def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--validated-tsv',required=True,type=Path);ap.add_argument('--prefix',required=True,type=Path);ap.add_argument('--summary-json',required=True,type=Path);args=ap.parse_args();print(json.dumps(write_tiers(args.validated_tsv,args.prefix,args.summary_json),indent=2,sort_keys=True))
if __name__=='__main__':main()
