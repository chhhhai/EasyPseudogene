#!/usr/bin/env python3
"""Filter EasyPseudogene broad candidates for same-species/self mode.

This module removes protein-to-genome candidates that overlap annotated functional
loci in the target GFF, including the query protein's own gene locus when it can
be inferred from gene/transcript/protein IDs. It rewrites the miniprot best-hit
files so downstream GeneWise/filtering steps operate only on off-annotation
homologous candidate regions.
"""
from __future__ import annotations

import argparse
import re
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

from annotation_io import aliases, open_text, parse_attributes, stable_id


def parse_attrs(s: str) -> Dict[str, str]:
    return parse_attributes(s)


def strip_version(x: str) -> str:
    return stable_id(x)


def ids_from_attrs(attrs: Dict[str, str]) -> List[str]:
    keys = ['ID','Name','Parent','protein_id','transcript_id','gene_id','gene','locus_tag','pacid','Dbxref']
    ids: List[str] = []
    for k in keys:
        v = attrs.get(k)
        if not v:
            continue
        for token in re.split(r'[,|]', v):
            token = token.strip()
            if not token:
                continue
            vals = [token, strip_version(token)]
            if ':' in token:
                vals += [token.split(':')[-1], strip_version(token.split(':')[-1])]
            ids.extend(vals)
    seen = set(); ret = []
    for x in ids:
        if x and x not in seen:
            seen.add(x); ret.append(x)
    return ret


def is_pseudogene_feature(ftype: str, attrs: Dict[str, str]) -> bool:
    text = ' '.join([ftype] + [str(v) for v in attrs.values()]).lower()
    pseudo_terms = [
        'pseudogene', 'pseudogenic', 'transcribed_pseudogene',
        'unitary_pseudogene', 'processed_pseudogene', 'duplicated_pseudogene'
    ]
    return any(term in text for term in pseudo_terms)


def merge_intervals(iv: Iterable[Tuple[int, int]]) -> List[Tuple[int, int]]:
    arr = sorted((int(s), int(e)) for s, e in iv if int(e) >= int(s))
    if not arr:
        return []
    merged = [arr[0]]
    for s, e in arr[1:]:
        ls, le = merged[-1]
        if s <= le + 1:
            merged[-1] = (ls, max(le, e))
        else:
            merged.append((s, e))
    return merged


def overlap_len(a_start: int, a_end: int, b_start: int, b_end: int) -> int:
    return max(0, min(a_end, b_end) - max(a_start, b_start) + 1)


def max_overlap_fraction(qs: int, qe: int, intervals: List[Tuple[int, int]]) -> Tuple[int, float]:
    # intervals are sorted/merged by start. Stop scanning once interval start is beyond query end.
    qlen = max(1, qe - qs + 1)
    ov = 0
    for s, e in intervals:
        if e < qs:
            continue
        if s > qe:
            break
        ov += overlap_len(qs, qe, s, e)
    return ov, ov / qlen


def load_gff(gff: Path, mask_features: set[str], keep_pseudogene_annotations: bool = False):
    mask_by_chr: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    pseudo_by_chr: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    locus_by_id: Dict[str, List[Tuple[str, int, int]]] = defaultdict(list)
    transcript_parent: Dict[str, str] = {}
    gene_loci: Dict[str, Tuple[str, int, int]] = {}
    pending_feature_ids = []

    with open_text(gff) as fh:
        for line in fh:
            if line.startswith('##FASTA') or line.startswith('>'):
                break
            if not line.strip() or line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9:
                continue
            chrom, _src, ftype, start, end, _score, _strand, _phase, attr_s = parts
            try:
                start_i, end_i = int(start), int(end)
            except ValueError:
                continue
            attrs = parse_attrs(attr_s)
            ids = ids_from_attrs(attrs)
            is_pseudo = is_pseudogene_feature(ftype, attrs)
            if is_pseudo:
                pseudo_by_chr[chrom].append((start_i, end_i))
            if ftype in mask_features and not (keep_pseudogene_annotations and is_pseudo):
                mask_by_chr[chrom].append((start_i, end_i))
            if ftype.lower() == 'gene':
                for x in ids:
                    gene_loci[x] = (chrom, start_i, end_i)
                    locus_by_id[x].append((chrom, start_i, end_i))
            elif ftype in {'mRNA', 'transcript'}:
                parent = attrs.get('Parent') or attrs.get('gene_id') or attrs.get('gene')
                for x in ids:
                    if parent:
                        transcript_parent[x] = parent
                pending_feature_ids.append((ids, chrom, start_i, end_i))
            elif ftype in {'CDS', 'exon'}:
                pending_feature_ids.append((ids, chrom, start_i, end_i))

    for ids, chrom, start_i, end_i in pending_feature_ids:
        candidate_loci = []
        for x in ids:
            parent = transcript_parent.get(x) or transcript_parent.get(strip_version(x))
            if parent:
                for key in {parent, strip_version(parent), parent.split(':')[-1]}:
                    if key in gene_loci:
                        candidate_loci.append(gene_loci[key])
            for key in {x, strip_version(x), x.split(':')[-1] if ':' in x else x}:
                if key in gene_loci:
                    candidate_loci.append(gene_loci[key])
        if not candidate_loci:
            candidate_loci = [(chrom, start_i, end_i)]
        for x in ids:
            for loc in candidate_loci:
                locus_by_id[x].append(loc)

    mask_by_chr = {c: merge_intervals(v) for c, v in mask_by_chr.items()}
    pseudo_by_chr = {c: merge_intervals(v) for c, v in pseudo_by_chr.items()}
    return mask_by_chr, pseudo_by_chr, locus_by_id


def main() -> None:
    ap = argparse.ArgumentParser(description='Filter broad miniprot candidates for EasyPseudogene self mode')
    ap.add_argument('--best-hits', required=True, type=Path)
    ap.add_argument('--best-hits-flanked', required=True, type=Path)
    ap.add_argument('--gff', required=True, type=Path)
    ap.add_argument('--out-best-hits', required=True, type=Path)
    ap.add_argument('--out-best-hits-flanked', required=True, type=Path)
    ap.add_argument('--report', required=True, type=Path)
    ap.add_argument('--bed', required=True, type=Path)
    ap.add_argument('--gff3-out', type=Path, help='Optional GFF3 output for kept broad self candidates')
    ap.add_argument('--removed', required=True, type=Path)
    ap.add_argument('--mask-features', default='gene,CDS')
    ap.add_argument('--annotated-overlap-cutoff', type=float, default=0.8)
    ap.add_argument('--self-overlap-cutoff', type=float, default=0.8)
    ap.add_argument('--keep-pseudogene-annotations', action='store_true', help='Do not mask GFF features whose type/attributes indicate pseudogene/pseudogenic annotation; retain candidates overlapping them.')
    ap.add_argument('--pseudogene-overlap-cutoff', type=float, default=0.1, help='Candidate-overlap fraction with pseudogene annotation required to force retain when --keep-pseudogene-annotations is set.')
    args = ap.parse_args()

    mask_features = {x.strip() for x in args.mask_features.split(',') if x.strip()}
    mask_by_chr, pseudo_by_chr, locus_by_id = load_gff(args.gff, mask_features, args.keep_pseudogene_annotations)

    # Screening emits matching rows in order. A query can have several loci,
    # so joining by query ID would substitute its last flank for every locus.
    hits = [line for line in args.best_hits.read_text().splitlines(keepends=True) if line.strip()]
    flanks = [line for line in args.best_hits_flanked.read_text().splitlines(keepends=True) if line.strip()]
    if len(hits) != len(flanks):
        raise ValueError('Unflanked and flanked hit counts differ')
    for line, flank in zip(hits, flanks):
        h, f = line.rstrip('\n').split('\t'), flank.rstrip('\n').split('\t')
        if len(h) < 8 or len(f) < 5 or h[:3] != f[:3]:
            raise ValueError('Mismatched or invalid flanked hit order')
        if not (0 <= int(f[3]) <= int(h[3]) < int(h[4]) <= int(f[4])):
            raise ValueError('Flanked interval does not contain its original hit')

    kept = []
    kept_flanked = []
    report_rows = []
    removed_rows = []
    bed_rows = []

    for line, flank in zip(hits, flanks):
        f = line.rstrip('\n').split('\t')
        if len(f) < 8:
            continue
        q, strand, chrom, start0, end0, tlen, alen, mapq = f[:8]
        s1 = min(int(start0), int(end0)) + 1
        e1 = max(int(start0), int(end0))
        q_keys = sorted(aliases(q))
        q_keys.append(q.split(':')[-1] if ':' in q else q)
        q_keys.append(strip_version(q_keys[-1]))
        mask_ov, mask_frac = max_overlap_fraction(s1, e1, mask_by_chr.get(chrom, []))
        pseudo_ov, pseudo_frac = max_overlap_fraction(s1, e1, pseudo_by_chr.get(chrom, []))
        self_frac = 0.0
        for key in q_keys:
            for loc_chrom, ls, le in locus_by_id.get(key, []):
                if loc_chrom != chrom:
                    continue
                ov, frac = max_overlap_fraction(s1, e1, [(ls, le)])
                if frac > self_frac:
                    self_frac = frac
        reason = 'kept'
        force_keep_pseudo = args.keep_pseudogene_annotations and pseudo_frac >= args.pseudogene_overlap_cutoff
        if force_keep_pseudo:
            reason = 'kept_pseudogene_annotation'
        elif self_frac >= args.self_overlap_cutoff:
            reason = 'removed_self_locus'
        elif mask_frac >= args.annotated_overlap_cutoff:
            reason = 'removed_annotated_locus'
        report_rows.append([q, chrom, str(s1), str(e1), strand, alen, mapq, f'{mask_frac:.4f}', f'{self_frac:.4f}', f'{pseudo_frac:.4f}', reason])
        if reason.startswith('kept'):
            kept.append(line)
            kept_flanked.append(flank)
            bed_rows.append([chrom, str(s1 - 1), str(e1), q, '0', strand, alen, mapq])
        else:
            removed_rows.append([q, chrom, str(s1), str(e1), strand, f'{mask_frac:.4f}', f'{self_frac:.4f}', reason])

    args.out_best_hits.write_text(''.join(kept))
    args.out_best_hits_flanked.write_text(''.join(x for x in kept_flanked if x))
    with args.report.open('w') as out:
        out.write('query\tchr\tstart\tend\tstrand\taln_len\tmapq\tannotated_overlap_frac\tself_overlap_frac\tpseudogene_overlap_frac\tdecision\n')
        for r in report_rows:
            out.write('\t'.join(r) + '\n')
    with args.removed.open('w') as out:
        out.write('query\tchr\tstart\tend\tstrand\tannotated_overlap_frac\tself_overlap_frac\treason\n')
        for r in removed_rows:
            out.write('\t'.join(r) + '\n')
    with args.bed.open('w') as out:
        for r in bed_rows:
            out.write('\t'.join(r) + '\n')
    if args.gff3_out:
        with args.gff3_out.open('w') as out:
            out.write('##gff-version 3\n')
            for idx, r in enumerate(bed_rows, 1):
                chrom, start0, end, q, score, strand, alen, mapq = r
                start1 = str(int(start0) + 1)
                attrs = f'ID={q}.broad_self_candidate_{idx};Name={q};Parent={q};candidate_type=broad_self_pseudogene_candidate;aln_len={alen};mapq={mapq}'
                out.write('\t'.join([chrom, 'EasyPseudogene', 'pseudogene_candidate', start1, end, '.', strand, '.', attrs]) + '\n')

    print(f'Self-mode filtering: kept {len(kept)} / {len(report_rows)} candidates; removed {len(removed_rows)}')
    print(f'  report: {args.report}')
    print(f'  broad BED: {args.bed}')
    if args.gff3_out:
        print(f'  broad GFF3: {args.gff3_out}')

if __name__ == '__main__':
    main()
