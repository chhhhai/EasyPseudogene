#!/usr/bin/env python3
"""
Diagnose why reference pseudogenes are missing from the results.
"""
import argparse
import re
from pathlib import Path
from collections import defaultdict


def parse_header(header):
    """Extract transcript_id and gene_symbol from FASTA header."""
    tid = header.split()[0].lstrip(">")
    m = re.search(r"gene_symbol=(\S+)", header)
    gene_symbol = m.group(1) if m else None
    return tid, gene_symbol


def read_fasta(fh):
    """Yield (header, sequence) tuples."""
    header = None
    seq_parts = []
    for line in fh:
        line = line.rstrip("\n")
        if line.startswith(">"):
            if header is not None:
                yield header, "".join(seq_parts)
            header = line
            seq_parts = []
        else:
            seq_parts.append(line)
    if header is not None:
        yield header, "".join(seq_parts)


def check_genewise_disruptions(gw_path, tids):
    """Check GeneWise output for stop codons and frameshifts."""
    results = {}  # tid -> {'stop': bool, 'shift': bool, 'score': float}
    if not gw_path.exists():
        return results
    
    with open(gw_path) as fh:
        block = []
        for line in fh:
            block.append(line)
            if line.strip() == '//':
                text = ''.join(block)
                for tid in tids:
                    if any(l.startswith('Query protein:') and tid in l for l in block):
                        has_stop = 'X' in text or '*:X[' in text or ':X[' in text
                        has_shift = '!' in text
                        score_match = re.search(r'Score\s+([0-9.]+)\s+bits', text)
                        score = float(score_match.group(1)) if score_match else None
                        results[tid] = {
                            'stop': has_stop,
                            'shift': has_shift,
                            'score': score
                        }
                block = []
    return results


def main():
    ap = argparse.ArgumentParser(description="Diagnose missing reference pseudogenes")
    ap.add_argument("--ref-list", required=True, help="Reference pseudogene list (one transcript_id per line)")
    ap.add_argument("--proteins", required=True, help="Input protein FASTA")
    ap.add_argument("--base-dir", default="02.pseudogene_screen", help="Base directory for results")
    ap.add_argument("--species", nargs="+", default=["dolphin", "sperm_whale"], help="Species to check")
    args = ap.parse_args()
    
    base = Path(args.base_dir)
    
    # Load reference list
    ref_tids = set()
    with open(args.ref_list) as fh:
        for line in fh:
            tid = line.strip()
            if tid:
                ref_tids.add(tid)
    
    # Load putative results
    putative = {}  # species -> set of tids
    putative_raw = {}  # species -> set of tids
    for sp in args.species:
        p_path = base / sp / "post_genewise_filter" / "putative_pseudogene_ids.txt"
        p_raw_path = base / sp / "post_genewise_filter" / "putative_pseudogene_ids.raw.txt"
        putative[sp] = set()
        putative_raw[sp] = set()
        if p_path.exists():
            with open(p_path) as fh:
                for line in fh:
                    tid = line.strip()
                    if tid:
                        putative[sp].add(tid)
        if p_raw_path.exists():
            with open(p_raw_path) as fh:
                for line in fh:
                    tid = line.strip()
                    if tid:
                        putative_raw[sp].add(tid)
    
    # Load prefilter valid transcripts
    prefilter_valid = set()
    prefilter_path = base.parent / "prefilter" / "valid_transcript_ids.txt"
    if prefilter_path.exists():
        with open(prefilter_path) as fh:
            for line in fh:
                tid = line.strip()
                if tid:
                    prefilter_valid.add(tid)
    
    # Load protein info
    transcripts = {}  # tid -> (gene_symbol, length, header)
    symbol_to_tids = defaultdict(list)  # gene_symbol -> list of (tid, length)
    
    with open(args.proteins) as fh:
        for header, seq in read_fasta(fh):
            tid, gene_symbol = parse_header(header)
            length = len(seq)
            transcripts[tid] = (gene_symbol, length, header)
            if gene_symbol:
                symbol_to_tids[gene_symbol].append((tid, length))
    
    # Determine which transcript was kept by prefilter for each gene_symbol
    symbol_to_kept = {}  # gene_symbol -> (kept_tid, kept_length)
    for gene_symbol, tids_with_lengths in symbol_to_tids.items():
        tids_with_lengths.sort(key=lambda x: (-x[1], x[0]))  # Sort by length desc, then tid
        kept_tid, kept_length = tids_with_lengths[0]
        symbol_to_kept[gene_symbol] = (kept_tid, kept_length)
    
    # Check miniprot hits
    miniprot_hits = {}  # species -> set of tids
    for sp in args.species:
        paf_path = base / sp / "miniprot" / f"{sp}.paf"
        miniprot_hits[sp] = set()
        if paf_path.exists():
            with open(paf_path) as fh:
                for line in fh:
                    if line.strip():
                        parts = line.rstrip('\n').split('\t')
                        if len(parts) >= 1:
                            miniprot_hits[sp].add(parts[0])
    
    # Check GeneWise disruptions
    genewise_results = {}  # species -> {tid: {stop, shift, score}}
    for sp in args.species:
        gw_path = base / sp / "genewise" / "pseudogene-genewise.txt"
        genewise_results[sp] = check_genewise_disruptions(gw_path, ref_tids)
    
    # Find missing refs
    union_putative = set()
    union_putative_raw = set()
    for sp in args.species:
        union_putative |= putative[sp]
        union_putative_raw |= putative_raw[sp]
    
    missing = ref_tids - union_putative
    missing_raw = ref_tids - union_putative_raw
    
    print("=" * 80)
    print("MISSING REFERENCE PSEUDOGENES DIAGNOSIS")
    print("=" * 80)
    print(f"\nTotal reference pseudogenes: {len(ref_tids)}")
    print(f"Missing from putative (after prefilter): {len(missing)}")
    print(f"Missing from raw putative (before prefilter): {len(missing_raw)}")
    print(f"Recall (putative): {len(ref_tids & union_putative) / len(ref_tids):.3f}")
    print(f"Recall (raw putative): {len(ref_tids & union_putative_raw) / len(ref_tids):.3f}")
    
    print("\n" + "=" * 80)
    print("DETAILED DIAGNOSIS FOR EACH MISSING REF")
    print("=" * 80)
    
    for tid in sorted(missing):
        gene_symbol, length, header = transcripts.get(tid, (None, 0, None))
        
        print(f"\n{tid}")
        print(f"  Gene symbol: {gene_symbol}")
        print(f"  Length: {length} aa")
        
        # Prefilter status
        in_prefilter = tid in prefilter_valid
        if gene_symbol:
            kept_tid, kept_length = symbol_to_kept.get(gene_symbol, (None, 0))
            if kept_tid == tid:
                prefilter_reason = "KEPT (longest transcript for this gene_symbol)"
            else:
                prefilter_reason = f"SKIPPED (prefilter kept {kept_tid} with length={kept_length})"
        else:
            prefilter_reason = "SKIPPED (no gene_symbol)"
        
        print(f"  Prefilter: {prefilter_reason}")
        
        # Raw putative status
        in_raw = {}
        for sp in args.species:
            in_raw[sp] = tid in putative_raw[sp]
        print(f"  Raw putative: {', '.join([f'{sp}={in_raw[sp]}' for sp in args.species])}")
        
        # Miniprot hits
        miniprot = {}
        for sp in args.species:
            miniprot[sp] = tid in miniprot_hits[sp]
        print(f"  Miniprot hits: {', '.join([f'{sp}={miniprot[sp]}' for sp in args.species])}")
        
        # GeneWise disruptions
        print(f"  GeneWise disruptions:")
        for sp in args.species:
            gw = genewise_results[sp].get(tid, {})
            stop = gw.get('stop', False)
            shift = gw.get('shift', False)
            score = gw.get('score', None)
            print(f"    {sp}: stop={stop} shift={shift} score={score}")
        
        # Summary reason
        reasons = []
        if not in_prefilter:
            reasons.append("Prefilter removed (not longest transcript or no gene_symbol)")
        if not any(in_raw.values()):
            reasons.append("GeneWise did not detect stop/shift in any species")
        elif not all(in_raw.values()):
            missing_sp = [sp for sp in args.species if not in_raw[sp]]
            reasons.append(f"GeneWise did not detect stop/shift in: {', '.join(missing_sp)}")
        
        print(f"  Summary: {'; '.join(reasons) if reasons else 'Unknown'}")
    
    print("\n" + "=" * 80)
    print("CATEGORIZATION")
    print("=" * 80)
    
    category_prefilter = []
    category_no_disruption = []
    category_partial = []
    
    for tid in sorted(missing):
        gene_symbol, length, header = transcripts.get(tid, (None, 0, None))
        in_prefilter = tid in prefilter_valid
        in_raw = {sp: tid in putative_raw[sp] for sp in args.species}
        
        if not in_prefilter:
            category_prefilter.append(tid)
        elif not any(in_raw.values()):
            category_no_disruption.append(tid)
        else:
            category_partial.append(tid)
    
    print(f"\n1. Removed by prefilter (not longest transcript): {len(category_prefilter)}")
    for tid in category_prefilter:
        print(f"   {tid}")
    
    print(f"\n2. No stop/shift detected in any species: {len(category_no_disruption)}")
    for tid in category_no_disruption:
        print(f"   {tid}")
    
    print(f"\n3. Detected in some but not all species: {len(category_partial)}")
    for tid in category_partial:
        in_raw = {sp: tid in putative_raw[sp] for sp in args.species}
        print(f"   {tid} ({', '.join([sp for sp, v in in_raw.items() if v])})")


if __name__ == "__main__":
    main()
