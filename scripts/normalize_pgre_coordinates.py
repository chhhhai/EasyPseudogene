#!/usr/bin/env python3
"""Convert P-GRe's headered 1-based inclusive export to BED6; never guess from .bed extension."""
import argparse,csv
from pathlib import Path

def convert(source,output):
    rows=[]
    with Path(source).open() as f:
        reader=csv.DictReader(f,delimiter='\t')
        if not reader.fieldnames or not {'chr','start','end'}.issubset(reader.fieldnames):
            raise ValueError('Expected headered P-GRe export with chr/start/end (1-based inclusive)')
        for i,r in enumerate(reader,1):
            start,end=int(r['start']),int(r['end'])
            if start<1 or end<start:raise ValueError(f'Invalid coordinates in row {i}')
            strand=r.get('strand','.')
            if strand not in ['+','-','.']:raise ValueError(f'Invalid strand in row {i}')
            rows.append([r['chr'],start-1,end,r.get('candidate_id') or r.get('id') or r.get('name') or f'PGRe_{i}',0,strand])
    with Path(output).open('x') as f:
        csv.writer(f,delimiter='\t',lineterminator='\n').writerows(rows)
    return len(rows)
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--input',required=True);p.add_argument('--output',required=True);a=p.parse_args();print(convert(a.input,a.output))
