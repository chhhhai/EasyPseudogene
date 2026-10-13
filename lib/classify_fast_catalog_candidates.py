#!/usr/bin/env python3
"""Classify EasyPseudogene fast-catalog candidates from miniprot/self-mode outputs.

Inputs are intentionally simple so this can run immediately after screening or
self-mode masking without GeneWise. It adds annotation/parent context, polyA scan,
optional repeat overlap, a coarse candidate type, and a confidence score.
"""
from __future__ import annotations
import argparse, re
from collections import defaultdict
from pathlib import Path
from typing import Dict


def read_fasta(path: Path) -> Dict[str, str]:
    seqs={}; name=None; chunks=[]
    with path.open(errors='replace') as fh:
        for line in fh:
            line=line.strip()
            if not line: continue
            if line.startswith('>'):
                if name is not None: seqs[name]=''.join(chunks).upper()
                name=line[1:].split()[0]; chunks=[]
            else: chunks.append(line)
    if name is not None: seqs[name]=''.join(chunks).upper()
    return seqs

def parse_attrs(s):
    d={}
    for part in s.split(';'):
        if not part: continue
        if '=' in part:
            k,v=part.split('=',1); d[k]=v
        elif ' ' in part:
            k,v=part.split(' ',1); d[k]=v.strip().strip('"')
    return d

def strip_ver(x): return re.sub(r'\.\d+$','',x)

def ids(attrs):
    vals=[]
    for k in ['ID','Name','Parent','protein_id','transcript_id','gene_id','gene','locus_tag','Dbxref']:
        v=attrs.get(k,'')
        for t in re.split(r'[,|]',v):
            t=t.strip()
            if not t: continue
            vals += [t, strip_ver(t)]
            if ':' in t: vals += [t.split(':')[-1], strip_ver(t.split(':')[-1])]
    seen=set(); out=[]
    for x in vals:
        if x and x not in seen: seen.add(x); out.append(x)
    return out

def ov(a,b,c,d): return max(0, min(b,d)-max(a,c)+1)
def frac_overlap(start,end, intervals):
    q=max(1,end-start+1); o=0
    for s,e in intervals:
        if e < start:
            continue
        if s > end:
            break
        o += ov(start,end,s,e)
    return o, o/q

def load_gff(gff: Path):
    exons=defaultdict(list); genes=defaultdict(list); id_to_gene={}; feature_by_chr=defaultdict(list)
    tx_parent={}
    if not gff or not gff.exists(): return exons, genes, id_to_gene, feature_by_chr
    for line in gff.open(errors='replace'):
        if not line.strip() or line.startswith('#'): continue
        f=line.rstrip('\n').split('\t')
        if len(f)<9: continue
        chrom,src,ftype,s,e,score,strand,phase,attr=f
        try: s=int(s); e=int(e)
        except: continue
        at=parse_attrs(attr); idlist=ids(at)
        if ftype.lower()=='gene':
            gid=at.get('ID') or at.get('Name') or (idlist[0] if idlist else f'{chrom}:{s}-{e}')
            for x in idlist: id_to_gene[x]=gid
            genes[gid].append((chrom,s,e,strand)); feature_by_chr[chrom].append((s,e,ftype,gid))
        elif ftype in {'mRNA','transcript'}:
            parent=at.get('Parent') or at.get('gene_id') or at.get('gene')
            for x in idlist:
                if parent: tx_parent[x]=parent; id_to_gene[x]=parent
            feature_by_chr[chrom].append((s,e,ftype,parent or 'NA'))
        elif ftype in {'CDS','exon'}:
            parent=at.get('Parent') or at.get('gene_id') or at.get('gene')
            gid=parent
            if parent in tx_parent: gid=tx_parent[parent]
            for x in idlist:
                if gid: id_to_gene[x]=gid
            if gid: exons[gid].append((chrom,s,e,strand))
            feature_by_chr[chrom].append((s,e,ftype,gid or 'NA'))
    return exons, genes, id_to_gene, feature_by_chr

def load_bed(path: Path):
    d=defaultdict(list)
    if not path or not path.exists(): return d
    for line in path.open(errors='replace'):
        if not line.strip() or line.startswith('#'): continue
        f=line.split('\t')
        if len(f)<3: continue
        try: d[f[0]].append((int(f[1])+1,int(f[2])))
        except: pass
    return d

def polyA(seq, chrom, start, end, strand, window=80, min_run=8):
    s=seq.get(chrom,'')
    if not s: return 0, 'none'
    if strand=='-':
        region=s[max(0,start-1-window):max(0,start-1)]
        chars='T'; side='upstream_T_on_minus'
    else:
        region=s[end:min(len(s),end+window)]
        chars='A'; side='downstream_A_on_plus'
    best=0; cur=0
    for c in region:
        if c==chars: cur+=1; best=max(best,cur)
        else: cur=0
    return (best if best>=min_run else 0), side


def revcomp(seq: str) -> str:
    tbl = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
    return seq.translate(tbl)[::-1].upper()

CODON_TABLE = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L','TCT':'S','TCC':'S','TCA':'S','TCG':'S',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','TGT':'C','TGC':'C','TGA':'*','TGG':'W',
    'CTT':'L','CTC':'L','CTA':'L','CTG':'L','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M','ACT':'T','ACC':'T','ACA':'T','ACG':'T',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K','AGT':'S','AGC':'S','AGA':'R','AGG':'R',
    'GTT':'V','GTC':'V','GTA':'V','GTG':'V','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'GAT':'D','GAC':'D','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G'
}

def translate_dna(dna: str, frame: int) -> str:
    aa=[]
    for i in range(frame, len(dna)-2, 3):
        aa.append(CODON_TABLE.get(dna[i:i+3].upper(), 'X'))
    return ''.join(aa)

def orf_lite(seq_by_chr, chrom, start, end, strand, qlen=0, max_nt=60000):
    """Lightweight ORF disruption proxy on the candidate genomic interval.

    This is intentionally conservative: it is most informative for intronless/processed
    candidates. For intron-containing duplicated candidates raw genomic DNA can contain
    intronic stops, so callers should treat this as supporting evidence, not a hard rule.
    """
    chrom_seq = seq_by_chr.get(chrom, '')
    if not chrom_seq:
        return dict(longest_orf_aa=0, orf_coverage=0.0, internal_stop_count=0, best_frame='NA', orf_disrupted=0)
    s=max(1,start); e=min(len(chrom_seq),end)
    if e < s:
        return dict(longest_orf_aa=0, orf_coverage=0.0, internal_stop_count=0, best_frame='NA', orf_disrupted=0)
    dna=chrom_seq[s-1:e].upper()
    if len(dna) > max_nt:
        # Keep the ends and center roughly; giant raw intervals are usually not ORF-lite-friendly.
        flank=max_nt//3
        mid_start=max(0, len(dna)//2 - flank//2)
        dna=dna[:flank] + dna[mid_start:mid_start+flank] + dna[-flank:]
    if strand == '-':
        dna=revcomp(dna)
    best=(0, 10**9, 0, 0)  # longest, stops, frame, aa_len
    for frame in range(3):
        aa=translate_dna(dna, frame)
        stops=aa.count('*')
        longest=max((len(x) for x in aa.split('*')), default=0)
        cand=(longest, -stops, frame, len(aa))
        if cand > (best[0], -best[1], best[2], best[3]):
            best=(longest, stops, frame, len(aa))
    longest, stops, frame, aa_len = best
    cov=longest/max(1, qlen) if qlen else 0.0
    # internal stops in best frame are disruption evidence; ORF much shorter than query is also suspicious.
    disrupted = 1 if (stops > 0 or (qlen and cov < 0.7)) else 0
    return dict(longest_orf_aa=longest, orf_coverage=round(cov,4), internal_stop_count=stops, best_frame=frame, orf_disrupted=disrupted)


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--best-hits', required=True, type=Path)
    ap.add_argument('--genome', required=True, type=Path)
    ap.add_argument('--gff', type=Path)
    ap.add_argument('--repeat-bed', type=Path)
    ap.add_argument('--out-tsv', required=True, type=Path)
    ap.add_argument('--out-gff3', required=True, type=Path)
    ap.add_argument('--polyA-window', type=int, default=80)
    ap.add_argument('--orf-lite-max-nt', type=int, default=60000, help='Maximum candidate nucleotides used for raw ORF-lite scan')
    args=ap.parse_args()
    seq=read_fasta(args.genome)
    exons, genes, id_to_gene, feats=load_gff(args.gff) if args.gff else ({},{},{},{})
    repeats=load_bed(args.repeat_bed) if args.repeat_bed else {}
    rows=[]
    for idx,line in enumerate(args.best_hits.open(errors='replace'),1):
        if not line.strip(): continue
        f=line.rstrip('\n').split('\t')
        if len(f)<8: continue
        q,strand,chrom,start0,end0,tlen,alen,mapq=f[:8]
        # Optional evidence columns emitted by newer run_screening.sh:
        # qlen,qstart,qend,nmatch,[unused alignment tags],fs,st,[unused CIGAR tags]
        qlen=qstart=qend=nmatch=0; fs=st=0
        if len(f) >= 19:
            try:
                qlen=int(float(f[8])); qstart=int(float(f[9])); qend=int(float(f[10])); nmatch=int(float(f[11]))
            except Exception:
                pass
            try: fs=int(float(f[15]))
            except Exception: fs=0
            try: st=int(float(f[16]))
            except Exception: st=0
        start=min(int(start0),int(end0))+1; end=max(int(start0),int(end0))
        alen=int(float(alen)); mapq=int(float(mapq)); length=end-start+1
        parent=id_to_gene.get(q) or id_to_gene.get(strip_ver(q)) or q
        parent_exons=exons.get(parent,[])
        parent_exon_n=len(parent_exons)
        # miniprot hit rows are collapsed candidate intervals; use rough structural proxies
        candidate_single_exon = length <= max(alen*4, 2000)
        intron_loss = 1 if parent_exon_n>=2 and candidate_single_exon else 0
        rep_ov, rep_frac=frac_overlap(start,end,repeats.get(chrom,[]))
        pa_len, pa_side=polyA(seq,chrom,start,end,strand,args.polyA_window)
        orf=orf_lite(seq, chrom, start, end, strand, qlen, args.orf_lite_max_nt)
        query_coverage=(qend-qstart)/max(1,qlen) if qlen else 0.0
        match_identity=nmatch/max(1,alen) if nmatch else 0.0
        disruption_count=fs+st
        if disruption_count>0 or (orf['orf_disrupted'] and query_coverage >= 0.3 and candidate_single_exon):
            ctype='disrupted_pseudogene_candidate'
        elif length < 500 or alen < 300:
            ctype='fragmentary_pseudogene_candidate'
        elif intron_loss and pa_len>=8:
            ctype='processed_pseudogene_candidate'
        elif intron_loss:
            ctype='processed_like_pseudogene_candidate'
        else:
            ctype='duplicated_pseudogene_candidate'
        score=0.0
        score += min(30.0, alen/50.0)
        score += min(15.0, mapq/3.0)
        score += min(20.0, query_coverage*20.0)
        score += min(20.0, disruption_count*10.0)
        score += 12.0 if orf['orf_disrupted'] and candidate_single_exon else 0.0
        score += min(8.0, max(0.0, (1.0 - float(orf['orf_coverage'])) * 8.0)) if qlen and candidate_single_exon else 0.0
        score += 10.0 if intron_loss else 0.0
        score += 10.0 if pa_len>=8 else 0.0
        score -= min(30.0, rep_frac*30.0)
        score=max(0.0, min(100.0, score))
        evid=[]
        if fs: evid.append(f'frameshift_{fs}')
        if st: evid.append(f'stop_{st}')
        if orf['orf_disrupted']: evid.append(f"orf_disrupted_longestAA_{orf['longest_orf_aa']}_cov_{orf['orf_coverage']}")
        if orf['internal_stop_count']: evid.append(f"orf_internal_stops_{orf['internal_stop_count']}")
        if query_coverage: evid.append(f'query_cov_{query_coverage:.3f}')
        if match_identity: evid.append(f'match_id_{match_identity:.3f}')
        if intron_loss: evid.append('intron_loss_like')
        if pa_len>=8: evid.append(f'polyA_{pa_len}')
        if rep_frac>0: evid.append(f'repeat_overlap_{rep_frac:.3f}')
        if not evid: evid.append('protein_to_genome_alignment')
        rows.append(dict(candidate_id=f'{q}.fastcat_{idx}',query_id=q,parent_gene=parent,chr=chrom,start=start,end=end,strand=strand,candidate_type=ctype,aln_len=alen,mapq=mapq,query_length=qlen,query_coverage=round(query_coverage,4),match_identity=round(match_identity,4),frameshift_count=fs,stop_count=st,longest_orf_aa=orf['longest_orf_aa'],orf_coverage=orf['orf_coverage'],orf_internal_stop_count=orf['internal_stop_count'],orf_disrupted=orf['orf_disrupted'],length=length,parent_exon_count=parent_exon_n,intron_loss=intron_loss,polyA_length=pa_len,polyA_side=pa_side,repeat_overlap_frac=round(rep_frac,4),confidence_score=round(score,3),evidence_tags=';'.join(evid)))
    fields=['candidate_id','query_id','parent_gene','chr','start','end','strand','candidate_type','aln_len','mapq','query_length','query_coverage','match_identity','frameshift_count','stop_count','longest_orf_aa','orf_coverage','orf_internal_stop_count','orf_disrupted','length','parent_exon_count','intron_loss','polyA_length','polyA_side','repeat_overlap_frac','confidence_score','evidence_tags']
    args.out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with args.out_tsv.open('w') as out:
        out.write('\t'.join(fields)+'\n')
        for r in rows: out.write('\t'.join(str(r[k]) for k in fields)+'\n')
    with args.out_gff3.open('w') as out:
        out.write('##gff-version 3\n')
        for r in rows:
            attrs=';'.join([f'ID={r["candidate_id"]}',f'Name={r["query_id"]}',f'Parent={r["parent_gene"]}',f'candidate_type={r["candidate_type"]}',f'confidence_score={r["confidence_score"]}',f'query_coverage={r["query_coverage"]}',f'frameshift_count={r["frameshift_count"]}',f'stop_count={r["stop_count"]}',f'orf_coverage={r["orf_coverage"]}',f'orf_disrupted={r["orf_disrupted"]}',f'evidence={r["evidence_tags"]}'])
            out.write('\t'.join([r['chr'],'EasyPseudogene-fast',r['candidate_type'],str(r['start']),str(r['end']),str(r['confidence_score']),r['strand'],'.',attrs])+'\n')
    print(f'Fast catalog candidates: {len(rows)}')
    print(f'  TSV: {args.out_tsv}')
    print(f'  GFF3: {args.out_gff3}')

if __name__=='__main__': main()
