"""Create collision-free per-locus jobs and the query-to-locus provenance table."""
import argparse
import csv
import hashlib
import re
from collections import Counter, defaultdict
from pathlib import Path
from prefilter_proteins import read_fasta


def prepare(proteins, hits, flanked, regions, jobs, paf):
    rows=[s.split('\t') for s in hits.read_text().splitlines() if s.strip()]
    flanks=[s.split('\t') for s in flanked.read_text().splitlines() if s.strip()]
    if len(rows)!=len(flanks):raise ValueError('Unflanked and flanked hit counts differ')
    counts=Counter(r[0] for r in rows)
    seqs={}
    with proteins.open() as f:
        for header,sequence in read_fasta(f):
            query=header[1:].split()[0]
            if query in seqs:raise ValueError('Duplicate FASTA sequence ID: '+query)
            seqs[query]=(header,sequence)
    aliases=defaultdict(list);used=set();mapping=[]
    with jobs.open('w') as out, (regions/'locus_proteins.fa').open('w') as fa:
        for r,f in zip(rows,flanks):
            query,strand,chrom,start,end=r[:5]
            if r[:3]!=f[:3]:raise ValueError('Mismatched flanked hit order')
            if not (0<=int(f[3])<=int(start)<int(end)<=int(f[4])):
                raise ValueError('Flanked interval does not contain its original hit')
            key='\t'.join(r[:5]);digest=hashlib.sha256(key.encode()).hexdigest()[:20]
            # Preserve ordinary single-hit query IDs; reserve the EPlocus_ namespace.
            locus=query if counts[query]==1 and re.fullmatch(r'[A-Za-z0-9_.-]+',query) and not query.startswith('EPlocus_') else 'EPlocus_'+digest
            if locus in used:continue # identical locus, not a second independent prediction
            used.add(locus);aliases[query].append(locus)
            h,seq=seqs[query];description=h[1:].partition(' ')[2]
            if not description:description='gene='+query
            # Query identity in GeneWise output becomes the locus ID.
            peptide=f'>{locus} {description}\n{seq}\n'
            (regions/(locus+'.pep.fa')).write_text(peptide);fa.write(peptide)
            out.write('\t'.join([query,strand,chrom,str(int(f[3])+1),f[4],locus])+'\n')
            mapping.append([locus,query,strand,chrom,start,end])
    with (regions/'locus_map.tsv').open('w') as out:
        w=csv.writer(out,delimiter='\t',lineterminator='\n');w.writerow(['locus_id','query_id','strand','chrom','start0','end0']);w.writerows(mapping)
    (regions/'id_map.tsv').write_text(''.join(f'{r[0]}\t{r[1]}\n' for r in mapping))
    with paf.open() as inp,(regions/'locus_alignments.paf').open('w') as out:
        for line in inp:
            q,sep,rest=line.partition('\t')
            for locus in aliases.get(q,[]):out.write(locus+sep+rest)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['proteins','hits','flanked','regions','jobs','paf']:p.add_argument('--'+name,required=True,type=Path)
    a=p.parse_args();prepare(**vars(a))
