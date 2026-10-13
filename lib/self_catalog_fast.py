#!/usr/bin/env python3
"""Integrated same-species pseudogene fast-catalog engine for EasyPseudogene.

This is a dedicated fast path for self-species cataloging. It avoids the old
GeneWise-oriented multi-file pipeline: parse miniprot PAF directly, retain top-N
hits per query, remove only the inferred source functional locus, compute CDS
annotation overlap and miniprot disruption evidence, optionally apply a score
cutoff, and write TSV/GFF3/summary outputs.
"""
from __future__ import annotations
import argparse, json, re
from collections import defaultdict, Counter
from pathlib import Path
from typing import Dict, List, Tuple, Any


def strip_ver(x: str) -> str:
    return re.sub(r"\.\d+$", "", x or "")


def parse_attrs(s: str) -> Dict[str, str]:
    d: Dict[str, str] = {}
    for part in s.split(';'):
        if not part:
            continue
        if '=' in part:
            k, v = part.split('=', 1)
            d[k] = v
        elif ' ' in part:
            k, v = part.split(' ', 1)
            d[k] = v.strip().strip('"')
    return d


def attr_ids(attrs: Dict[str, str]) -> List[str]:
    vals: List[str] = []
    keys = ['ID', 'Name', 'Parent', 'protein_id', 'transcript_id', 'gene_id', 'gene', 'locus_tag', 'Dbxref']
    for k in keys:
        v = attrs.get(k, '')
        for t in re.split(r'[,|]', v):
            t = t.strip()
            if not t:
                continue
            vals.append(t)
            vals.append(strip_ver(t))
            if ':' in t:
                vals.append(t.split(':')[-1])
                vals.append(strip_ver(t.split(':')[-1]))
    out: List[str] = []
    seen = set()
    for x in vals:
        if x and x not in seen:
            seen.add(x); out.append(x)
    return out


def merge_intervals(v: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    if not v:
        return []
    v = sorted(v)
    out = [v[0]]
    for s, e in v[1:]:
        ls, le = out[-1]
        if s <= le + 1:
            out[-1] = (ls, max(le, e))
        else:
            out.append((s, e))
    return out


def overlap_bp(start: int, end: int, intervals: List[Tuple[int, int]]) -> int:
    total = 0
    for s, e in intervals:
        if e < start:
            continue
        if s > end:
            break
        total += max(0, min(end, e) - max(start, s) + 1)
    return total


def max_overlap_frac(start: int, end: int, intervals: List[Tuple[int, int]]) -> Tuple[int, float]:
    q = max(1, end - start + 1)
    best = 0
    for s, e in intervals:
        if e < start:
            continue
        if s > end:
            break
        ov = max(0, min(end, e) - max(start, s) + 1)
        if ov > best:
            best = ov
    return best, best / q


def load_gff(gff: Path):
    cds_by_chr: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    gene_loci_by_id: Dict[str, List[Tuple[str, int, int, str]]] = defaultdict(list)
    exon_count_by_gene: Counter[str] = Counter()
    transcript_parent: Dict[str, str] = {}
    pending_exons: List[Tuple[List[str], str, int, int, str]] = []

    for line in gff.open(errors='replace'):
        if not line.strip() or line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9:
            continue
        chrom, src, ftype, s, e, score, strand, phase, attr_s = f
        try:
            s_i, e_i = int(s), int(e)
        except Exception:
            continue
        attrs = parse_attrs(attr_s)
        ids = attr_ids(attrs)
        ftl = ftype.lower()
        if ftl == 'gene':
            locus = (chrom, s_i, e_i, strand)
            for x in ids:
                gene_loci_by_id[x].append(locus)
            # common Arabidopsis shortcut: AT1G01010 from gene:AT1G01010
            for x in list(ids):
                if ':' in x:
                    gene_loci_by_id[x.split(':')[-1]].append(locus)
                    gene_loci_by_id[strip_ver(x.split(':')[-1])].append(locus)
        elif ftl in {'mrna', 'transcript', 'rrna', 'trna'}:
            parent = attrs.get('Parent') or attrs.get('gene_id') or attrs.get('gene') or ''
            for x in ids:
                if parent:
                    transcript_parent[x] = parent
                    transcript_parent[strip_ver(x)] = parent
        elif ftl in {'cds', 'exon'}:
            parent = attrs.get('Parent') or attrs.get('gene_id') or attrs.get('gene') or ''
            pending_exons.append((ids + ([parent] if parent else []), chrom, s_i, e_i, strand))
            if ftl == 'cds':
                cds_by_chr[chrom].append((s_i, e_i))

    for ids, chrom, s_i, e_i, strand in pending_exons:
        gid = None
        for x in ids:
            gid = transcript_parent.get(x) or transcript_parent.get(strip_ver(x))
            if gid:
                break
        if gid:
            exon_count_by_gene[gid] += 1
            exon_count_by_gene[strip_ver(gid)] += 1
            if ':' in gid:
                exon_count_by_gene[gid.split(':')[-1]] += 1

    cds_by_chr = {c: merge_intervals(v) for c, v in cds_by_chr.items()}
    for k, loci in list(gene_loci_by_id.items()):
        # Deduplicate while preserving list shape.
        seen = set(); uniq = []
        for loc in loci:
            if loc not in seen:
                seen.add(loc); uniq.append(loc)
        gene_loci_by_id[k] = uniq
    return cds_by_chr, gene_loci_by_id, exon_count_by_gene



def load_bed_intervals(path: Path | None) -> Dict[str, List[Tuple[int, int]]]:
    if not path:
        return {}
    out: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    with path.open(errors='replace') as fh:
        for line in fh:
            if not line.strip() or line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 3:
                continue
            try:
                chrom = f[0]
                start = int(float(f[1])) + 1
                end = int(float(f[2]))
            except Exception:
                continue
            if end >= start:
                out[chrom].append((start, end))
    return {c: merge_intervals(v) for c, v in out.items()}


def parse_paf_tags(fields: List[str]) -> Dict[str, Any]:
    tags: Dict[str, Any] = {'AS': None, 'ms': None, 'np': None, 'fs': 0, 'st': 0, 'cg': '', 'cs': ''}
    for x in fields[12:]:
        if len(x) < 5 or ':' not in x:
            continue
        parts = x.split(':', 2)
        if len(parts) != 3:
            continue
        key, typ, val = parts
        if key in {'AS', 'ms', 'np', 'fs', 'st'}:
            try: tags[key] = int(val)
            except Exception: tags[key] = val
        elif key in {'cg', 'cs'}:
            tags[key] = val
    return tags


def read_paf_topn(paf: Path, max_hits: int, min_aln_len: int, min_qcov: float, min_mapq: int):
    by_q: Dict[str, List[Dict[str, Any]]] = defaultdict(list)
    total = 0; passed_basic = 0
    for line in paf.open(errors='replace'):
        if not line.strip():
            continue
        total += 1
        f = line.rstrip('\n').split('\t')
        if len(f) < 12:
            continue
        try:
            q, qlen, qs, qe, strand, chrom, tlen, ts, te, nmatch, alen, mapq = (
                f[0], int(f[1]), int(f[2]), int(f[3]), f[4], f[5], int(f[6]), int(f[7]), int(f[8]), int(f[9]), int(f[10]), int(f[11])
            )
        except Exception:
            continue
        qcov = (qe - qs) / max(1, qlen)
        if alen < min_aln_len or qcov < min_qcov or mapq < min_mapq:
            continue
        passed_basic += 1
        tags = parse_paf_tags(f)
        hit = dict(query=q, qlen=qlen, qstart=qs, qend=qe, strand=strand, chrom=chrom, tlen=tlen,
                   start0=min(ts, te), end0=max(ts, te), nmatch=nmatch, alen=alen, mapq=mapq,
                   query_coverage=qcov, match_identity=nmatch / max(1, alen), **tags)
        # Prefer longer/covered/disrupted hits, stable by input order implicitly.
        score = (alen, qcov, int(tags.get('fs') or 0) + int(tags.get('st') or 0), nmatch)
        hit['_rank'] = score
        by_q[q].append(hit)
    hits = []
    for q, arr in by_q.items():
        arr.sort(key=lambda h: h['_rank'], reverse=True)
        hits.extend(arr[:max_hits])
    return hits, {'paf_rows': total, 'passed_basic_filters': passed_basic, 'topn_hits': len(hits), 'queries_with_hits': len(by_q)}


def source_overlap(hit, gene_loci_by_id) -> Tuple[int, float]:
    q = hit['query']
    keys = [q, strip_ver(q), q.split(':')[-1] if ':' in q else q, strip_ver(q.split(':')[-1] if ':' in q else q)]
    start = hit['start0'] + 1; end = hit['end0']; chrom = hit['chrom']
    best_ov = 0; best_frac = 0.0
    for k in keys:
        for loc_chrom, s, e, strand in gene_loci_by_id.get(k, []):
            if loc_chrom != chrom:
                continue
            ov, frac = max_overlap_frac(start, end, [(s, e)])
            if frac > best_frac:
                best_ov, best_frac = ov, frac
    return best_ov, best_frac


def classify_and_score(hit, cds_frac: float, src_frac: float, exon_count: int, repeat_frac: float = 0.0) -> Tuple[str, float, List[str]]:
    qcov = hit['query_coverage']; ident = hit['match_identity']; fs = int(hit.get('fs') or 0); st = int(hit.get('st') or 0)
    alen = hit['alen']; length = hit['end0'] - hit['start0']
    disruption = fs + st
    evidence = [f'query_cov_{qcov:.3f}', f'match_id_{ident:.3f}', f'cds_overlap_{cds_frac:.3f}', f'repeat_overlap_{repeat_frac:.3f}']
    if fs: evidence.append(f'frameshift_{fs}')
    if st: evidence.append(f'stop_{st}')
    if src_frac: evidence.append(f'source_overlap_{src_frac:.3f}')
    # A rough processed-like proxy: genomic footprint is compact relative to protein alignment.
    compact = length <= max(alen * 4, 2000)
    intron_loss_like = 1 if exon_count >= 2 and compact else 0
    if intron_loss_like: evidence.append('intron_loss_like')

    if disruption > 0:
        ctype = 'disrupted_pseudogene_candidate'
    elif alen < 300 or qcov < 0.35:
        ctype = 'fragmentary_pseudogene_candidate'
    elif intron_loss_like:
        ctype = 'processed_like_pseudogene_candidate'
    else:
        ctype = 'duplicated_pseudogene_candidate'

    score = 0.0
    score += min(25.0, qcov * 25.0)
    score += min(20.0, ident * 20.0)
    score += min(25.0, disruption * 12.5)
    score += 10.0 if intron_loss_like else 0.0
    score += 10.0 if cds_frac < 0.2 else 0.0
    score -= min(30.0, cds_frac * 30.0) if disruption == 0 else min(10.0, cds_frac * 10.0)
    score -= min(25.0, repeat_frac * 25.0)
    score = max(0.0, min(100.0, score))
    return ctype, score, evidence



def cluster_locus_rows(rows: List[Dict[str, Any]], max_gap: int = 0, by_strand: bool = False, coordinate_mode: str = "best") -> Tuple[List[Dict[str, Any]], Counter]:
    """Merge overlapping/nearby protein hits into genomic pseudogene loci.

    Parent selection is performed at the cluster level: the representative parent/query is
    the row with the strongest score, then disruption evidence, query coverage, and
    alignment length. Coordinates are the cluster union, so the output represents a locus
    rather than a single protein-to-genome fragment.
    """
    counts: Counter = Counter()
    if not rows:
        return rows, counts
    grouped: Dict[Tuple[str, str], List[Dict[str, Any]]] = defaultdict(list)
    for r in rows:
        key = (r['chr'], r['strand'] if by_strand else '.')
        grouped[key].append(r)
    clustered: List[Dict[str, Any]] = []
    cluster_id = 0
    for key, arr in grouped.items():
        arr.sort(key=lambda r: (int(r['start']), int(r['end'])))
        cur: List[Dict[str, Any]] = []
        cur_end = -1
        for r in arr:
            s_i, e_i = int(r['start']), int(r['end'])
            if not cur or s_i <= cur_end + max_gap:
                cur.append(r)
                cur_end = max(cur_end, e_i)
            else:
                cluster_id += 1
                clustered.append(collapse_cluster(cur, cluster_id, coordinate_mode))
                cur = [r]
                cur_end = e_i
        if cur:
            cluster_id += 1
            clustered.append(collapse_cluster(cur, cluster_id, coordinate_mode))
    counts['clusters'] = len(clustered)
    counts['clustered_input_rows'] = len(rows)
    counts['multi_hit_clusters'] = sum(1 for r in clustered if int(r.get('cluster_size', 1)) > 1)
    return clustered, counts


def collapse_cluster(cur: List[Dict[str, Any]], cluster_id: int, coordinate_mode: str = "best") -> Dict[str, Any]:
    def rank(r: Dict[str, Any]):
        disruption = int(r.get('frameshift_count', 0)) + int(r.get('stop_count', 0))
        return (float(r.get('confidence_score', 0)), disruption, float(r.get('query_coverage', 0)), int(r.get('aln_len', 0)), -float(r.get('cds_overlap_frac', 0)))
    best = max(cur, key=rank)
    out = dict(best)
    if coordinate_mode == 'union':
        out['start'] = min(int(r['start']) for r in cur)
        out['end'] = max(int(r['end']) for r in cur)
    else:
        # Use the best parent/query alignment as the representative interval; this
        # removes duplicate locus calls without inflating intervals.
        out['start'] = int(best['start'])
        out['end'] = int(best['end'])
    out['candidate_id'] = f"locus_cluster_{cluster_id:06d}"
    out['cluster_size'] = len(cur)
    out['cluster_best_query'] = best['query_id']
    parents = []
    seen = set()
    for r in sorted(cur, key=rank, reverse=True):
        p = str(r.get('parent_gene', ''))
        if p and p not in seen:
            seen.add(p); parents.append(p)
        if len(parents) >= 5:
            break
    out['alternative_parents'] = ','.join(parents)
    out['cluster_member_ids'] = ','.join(str(r['candidate_id']) for r in cur[:50])
    # Aggregate evidence in a compact form.
    out['frameshift_count'] = max(int(r.get('frameshift_count', 0)) for r in cur)
    out['stop_count'] = max(int(r.get('stop_count', 0)) for r in cur)
    out['aln_len'] = max(int(r.get('aln_len', 0)) for r in cur)
    out['query_coverage'] = max(float(r.get('query_coverage', 0)) for r in cur)
    out['match_identity'] = max(float(r.get('match_identity', 0)) for r in cur)
    out['cds_overlap_frac'] = round(min(float(r.get('cds_overlap_frac', 0)) for r in cur), 4)
    out['repeat_overlap_frac'] = round(min(float(r.get('repeat_overlap_frac', 0)) for r in cur), 4)
    out['confidence_score'] = max(float(r.get('confidence_score', 0)) for r in cur)
    type_priority = {
        'disrupted_pseudogene_candidate': 4,
        'processed_like_pseudogene_candidate': 3,
        'fragmentary_pseudogene_candidate': 2,
        'duplicated_pseudogene_candidate': 1,
    }
    out['candidate_type'] = max((r['candidate_type'] for r in cur), key=lambda t: type_priority.get(t, 0))
    tags = []
    for r in cur[:20]:
        tags.extend(str(r.get('evidence_tags', '')).split(';'))
    compact_tags = []
    seen_tags = set()
    for t in tags:
        if t and t not in seen_tags:
            seen_tags.add(t); compact_tags.append(t)
    compact_tags.append(f'cluster_size_{len(cur)}')
    compact_tags.append(f'best_parent_{best["parent_gene"]}')
    out['evidence_tags'] = ';'.join(compact_tags[:40])
    return out


def main():
    ap = argparse.ArgumentParser(description='Integrated EasyPseudogene self-species fast catalog engine')
    ap.add_argument('--paf', required=True, type=Path)
    ap.add_argument('--gff', required=True, type=Path)
    ap.add_argument('--out-tsv', required=True, type=Path)
    ap.add_argument('--out-gff3', required=True, type=Path)
    ap.add_argument('--summary-json', required=True, type=Path)
    ap.add_argument('--max-hits-per-protein', type=int, default=5)
    ap.add_argument('--min-aln-len', type=int, default=30)
    ap.add_argument('--min-query-coverage', type=float, default=0.05)
    ap.add_argument('--min-mapq', type=int, default=0)
    ap.add_argument('--source-overlap-cutoff', type=float, default=0.8)
    ap.add_argument('--score-cutoff', type=float, default=0.0)
    ap.add_argument('--require-evidence', action='store_true', help='Require disruption evidence or low CDS overlap/fragmentary evidence')
    ap.add_argument('--max-cds-overlap', type=float, default=1.01, help='Discard candidates with CDS overlap fraction above this value [disabled]')
    ap.add_argument('--repeat-bed', type=Path, default=None, help='Optional RepeatMasker/repeat BED for repeat-overlap penalty')
    ap.add_argument('--cluster-loci', action='store_true', help='Cluster overlapping retained hits into locus-level candidates and select one parent/query per locus')
    ap.add_argument('--cluster-max-gap', type=int, default=0, help='Maximum bp gap for locus clustering [0]')
    ap.add_argument('--cluster-by-strand', action='store_true', help='Cluster hits separately by strand')
    ap.add_argument('--cluster-coordinate-mode', choices=['best','union'], default='best', help='Use best representative interval or union span for clustered loci [best]')
    args = ap.parse_args()

    cds_by_chr, gene_loci_by_id, exon_count_by_gene = load_gff(args.gff)
    repeat_by_chr = load_bed_intervals(args.repeat_bed)
    hits, summary = read_paf_topn(args.paf, args.max_hits_per_protein, args.min_aln_len, args.min_query_coverage, args.min_mapq)
    rows = []
    counts = Counter()
    for idx, h in enumerate(hits, 1):
        start = h['start0'] + 1; end = h['end0']; length = end - start + 1
        src_ov, src_frac = source_overlap(h, gene_loci_by_id)
        if src_frac >= args.source_overlap_cutoff:
            counts['removed_source_locus'] += 1
            continue
        cds_ov = overlap_bp(start, end, cds_by_chr.get(h['chrom'], []))
        cds_frac = cds_ov / max(1, length)
        if cds_frac > args.max_cds_overlap:
            counts['removed_cds_overlap_cutoff'] += 1
            continue
        repeat_ov = overlap_bp(start, end, repeat_by_chr.get(h['chrom'], [])) if repeat_by_chr else 0
        repeat_frac = repeat_ov / max(1, length)
        parent = strip_ver(h['query'])
        exon_count = max(exon_count_by_gene.get(h['query'], 0), exon_count_by_gene.get(parent, 0))
        ctype, score, evidence = classify_and_score(h, cds_frac, src_frac, exon_count, repeat_frac)
        if score < args.score_cutoff:
            counts['removed_score_cutoff'] += 1
            continue
        if args.require_evidence:
            has_evidence = (int(h.get('fs') or 0) + int(h.get('st') or 0) > 0) or cds_frac < 0.2 or ctype == 'fragmentary_pseudogene_candidate'
            if not has_evidence:
                counts['removed_no_evidence'] += 1
                continue
        counts['kept'] += 1
        counts[f'type:{ctype}'] += 1
        row = dict(
            candidate_id=f'{h["query"]}.selfcat_{idx}', query_id=h['query'], parent_gene=parent,
            chr=h['chrom'], start=start, end=end, strand=h['strand'], candidate_type=ctype,
            query_length=h['qlen'], query_coverage=round(h['query_coverage'], 4), match_identity=round(h['match_identity'], 4),
            aln_len=h['alen'], mapq=h['mapq'], frameshift_count=int(h.get('fs') or 0), stop_count=int(h.get('st') or 0),
            cds_overlap_frac=round(cds_frac, 4), repeat_overlap_frac=round(repeat_frac, 4), source_overlap_frac=round(src_frac, 4), parent_exon_count=exon_count,
            confidence_score=round(score, 3), evidence_tags=';'.join(evidence)
        )
        rows.append(row)

    if args.cluster_loci:
        rows, cluster_counts = cluster_locus_rows(rows, args.cluster_max_gap, args.cluster_by_strand, args.cluster_coordinate_mode)
        counts.update(cluster_counts)

    fields = ['candidate_id','query_id','parent_gene','chr','start','end','strand','candidate_type','query_length','query_coverage','match_identity','aln_len','mapq','frameshift_count','stop_count','cds_overlap_frac','repeat_overlap_frac','source_overlap_frac','parent_exon_count','confidence_score','cluster_size','cluster_best_query','alternative_parents','cluster_member_ids','evidence_tags']
    for r in rows:
        r.setdefault('cluster_size', 1)
        r.setdefault('cluster_best_query', r.get('query_id', ''))
        r.setdefault('alternative_parents', r.get('parent_gene', ''))
        r.setdefault('cluster_member_ids', r.get('candidate_id', ''))
    args.out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with args.out_tsv.open('w') as out:
        out.write('\t'.join(fields)+'\n')
        for r in rows:
            out.write('\t'.join(str(r[k]) for k in fields)+'\n')
    with args.out_gff3.open('w') as out:
        out.write('##gff-version 3\n')
        for r in rows:
            attrs = ';'.join([f'ID={r["candidate_id"]}', f'Name={r["query_id"]}', f'Parent={r["parent_gene"]}', f'candidate_type={r["candidate_type"]}', f'confidence_score={r["confidence_score"]}', f'query_coverage={r["query_coverage"]}', f'match_identity={r["match_identity"]}', f'frameshift_count={r["frameshift_count"]}', f'stop_count={r["stop_count"]}', f'cds_overlap_frac={r["cds_overlap_frac"]}', f'repeat_overlap_frac={r["repeat_overlap_frac"]}', f'cluster_size={r["cluster_size"]}', f'best_query={r["cluster_best_query"]}', f'alternative_parents={r["alternative_parents"]}', f'evidence={r["evidence_tags"]}'])
            out.write('\t'.join([r['chr'], 'EasyPseudogene-selfcat-fast', r['candidate_type'], str(r['start']), str(r['end']), str(r['confidence_score']), r['strand'], '.', attrs])+'\n')
    summary.update({k:int(v) for k,v in counts.items()})
    summary['output_candidates'] = len(rows)
    summary['out_tsv'] = str(args.out_tsv)
    summary['out_gff3'] = str(args.out_gff3)
    args.summary_json.write_text(json.dumps(summary, indent=2, sort_keys=True) + '\n')
    print(json.dumps(summary, indent=2, sort_keys=True))

if __name__ == '__main__':
    main()
