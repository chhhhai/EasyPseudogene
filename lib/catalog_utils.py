#!/usr/bin/env python3
"""Shared interval, FASTA, metadata, and atomic-I/O helpers."""
from __future__ import annotations
import bisect, contextlib, os, re, tempfile
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Mapping, TextIO, Tuple
HalfOpenInterval = Tuple[int, int]

def merge_half_open_intervals(intervals: Iterable[HalfOpenInterval]) -> List[HalfOpenInterval]:
    arr = sorted((int(s), int(e)) for s, e in intervals if int(e) > int(s))
    if not arr: return []
    merged=[arr[0]]
    for start,end in arr[1:]:
        ls,le=merged[-1]
        if start <= le: merged[-1]=(ls,max(le,end))
        else: merged.append((start,end))
    return merged

class IntervalIndex:
    def __init__(self, intervals_by_chrom: Mapping[str, Iterable[HalfOpenInterval]] | None = None):
        self.intervals={}; self.starts={}
        for chrom,intervals in (intervals_by_chrom or {}).items():
            merged=merge_half_open_intervals(intervals)
            self.intervals[str(chrom)]=merged
            self.starts[str(chrom)]=[s for s,_ in merged]
    def overlap_bp(self, chrom: str, start0: int, end0: int) -> int:
        if end0 <= start0 or chrom not in self.intervals: return 0
        intervals=self.intervals[chrom]; starts=self.starts[chrom]
        i=bisect.bisect_left(starts,end0)
        j=max(0,bisect.bisect_right(starts,start0)-1)
        total=0
        while j < i:
            start,end=intervals[j]
            if end > start0: total += max(0,min(end0,end)-max(start0,start))
            j += 1
        return total
    def overlap_fraction(self, chrom: str, start0: int, end0: int) -> float:
        return self.overlap_bp(chrom,start0,end0)/max(1,end0-start0)

def union_coverage_fraction(intervals: Iterable[Tuple[int,int]], query_length: int, *, one_based_inclusive: bool=False) -> float:
    half=[]
    for start,end in intervals:
        lower,upper=sorted((int(start),int(end)))
        start_i=lower-1 if one_based_inclusive else lower
        end_i=upper
        if end_i > start_i: half.append((start_i,end_i))
    covered=sum(e-s for s,e in merge_half_open_intervals(half))
    return min(1.0,covered/max(1,int(query_length)))

def infer_gene_id(sequence_id: str, description: str) -> str:
    text=description or ''
    patterns=[
        r'(?:^|[\s\[])gene[=:]([^\s;\]]+)',
        r'(?:^|[\s\[])gene_id[=:]([^\s;\]]+)',
        r'(?:^|[\s\[])locus_tag[=:]([^\s;\]]+)',
        r'(?:^|[\s\[])gene_symbol=([^\s;\]]+)',
        r'(?:^|\s)GN=([^\s;]+)',
    ]
    for pattern in patterns:
        m=re.search(pattern,text,flags=re.IGNORECASE)
        if m: return m.group(1).strip('"[]')
    return re.sub(r'\.\d+$','',sequence_id)

def read_fasta(path: Path) -> Iterator[Tuple[str,str,str]]:
    name=None; desc=None; seq=[]
    with Path(path).open(errors='replace') as handle:
        for line in handle:
            line=line.rstrip('\n')
            if line.startswith('>'):
                if name is not None: yield name,desc or name,''.join(seq)
                desc=line[1:]; name=desc.split()[0]; seq=[]
            else: seq.append(line.strip())
        if name is not None: yield name,desc or name,''.join(seq)

def write_fasta(records, path: Path, width: int=80) -> None:
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    with atomic_text_writer(path) as out:
        for _name,desc,seq in records:
            out.write('>'+desc+'\n')
            for i in range(0,len(seq),width): out.write(seq[i:i+width]+'\n')

def load_protein_metadata(path: Path) -> Dict[str,Dict[str,object]]:
    return {name:{'length':len(seq),'gene_id':infer_gene_id(name,desc),'description':desc} for name,desc,seq in read_fasta(Path(path))}

def file_sha256(path,chunk_size=1024*1024):
 import hashlib
 digest=hashlib.sha256()
 with Path(path).open('rb') as handle:
  for chunk in iter(lambda:handle.read(chunk_size),b''):digest.update(chunk)
 return digest.hexdigest()

@contextlib.contextmanager
def atomic_text_writer(path: Path) -> Iterator[TextIO]:
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    fd,tmp_name=tempfile.mkstemp(prefix=path.name+'.tmp.',dir=str(path.parent),text=True)
    tmp=Path(tmp_name)
    try:
        with os.fdopen(fd,'w') as handle:
            yield handle; handle.flush(); os.fsync(handle.fileno())
        os.replace(tmp,path)
    except Exception:
        try: os.close(fd)
        except OSError: pass
        try: tmp.unlink()
        except FileNotFoundError: pass
        raise
