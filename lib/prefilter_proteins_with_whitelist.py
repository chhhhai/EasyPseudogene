#!/usr/bin/env python3
"""
Prefilter protein FASTA with whitelist support:
1. Always keep transcripts in whitelist
2. For other transcripts: only keep transcripts with a non-empty gene_symbol
3. If multiple transcripts share the same gene_symbol, keep:
   - All whitelisted transcripts for that gene_symbol
   - The longest non-whitelisted transcript (if any)
"""
import argparse
import re
import sys
from pathlib import Path


def parse_header(header):
    """Extract transcript_id and gene_symbol from FASTA header."""
    tid = header.split()[0].lstrip(">")
    m = re.search(r"gene_symbol=(\S+)", header)
    gene_symbol = m.group(1) if m else ""
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


def load_whitelist(whitelist_path):
    """Load whitelist of transcript IDs."""
    if not whitelist_path or not Path(whitelist_path).exists():
        return set()
    whitelist = set()
    with open(whitelist_path) as fh:
        for line in fh:
            tid = line.strip()
            if tid:
                whitelist.add(tid)
    return whitelist


def main():
    ap = argparse.ArgumentParser(description="Prefilter proteins with whitelist support")
    ap.add_argument("--input", "-i", required=True, help="Input protein FASTA")
    ap.add_argument("--output", "-o", required=True, help="Output filtered FASTA")
    ap.add_argument("--report", "-r", help="Optional TSV report of kept/skipped")
    ap.add_argument("--whitelist", "-w", help="Optional whitelist file (one transcript_id per line)")
    args = ap.parse_args()

    whitelist = load_whitelist(args.whitelist)

    # First pass: collect all transcripts grouped by gene_symbol
    by_symbol = {}  # gene_symbol -> list of (tid, header, seq, length, is_whitelisted)
    no_symbol = []  # transcripts without gene_symbol
    whitelisted_no_symbol = []  # whitelisted transcripts without gene_symbol

    with open(args.input, "r") as fh:
        for header, seq in read_fasta(fh):
            tid, gene_symbol = parse_header(header)
            length = len(seq)
            is_whitelisted = tid in whitelist
            
            if not gene_symbol or gene_symbol == "":
                if is_whitelisted:
                    whitelisted_no_symbol.append((tid, header, seq, length))
                else:
                    no_symbol.append((tid, header, seq, length))
            else:
                by_symbol.setdefault(gene_symbol, []).append((tid, header, seq, length, is_whitelisted))

    # For each gene_symbol, keep:
    # - All whitelisted transcripts
    # - The longest non-whitelisted transcript (if any)
    kept = []
    skipped_dup = []
    
    for gene_symbol, transcripts in sorted(by_symbol.items()):
        # Separate whitelisted and non-whitelisted
        whitelisted = [t for t in transcripts if t[4]]
        non_whitelisted = [t for t in transcripts if not t[4]]
        
        # Always keep all whitelisted transcripts
        for tid, header, seq, length, _ in whitelisted:
            kept.append((tid, header, seq, gene_symbol, "whitelisted"))
        
        # Keep the longest non-whitelisted transcript
        if non_whitelisted:
            non_whitelisted.sort(key=lambda x: (-x[3], x[0]))  # Sort by length desc, then tid
            best = non_whitelisted[0]
            kept.append((best[0], best[1], best[2], gene_symbol, "longest"))
            for t in non_whitelisted[1:]:
                skipped_dup.append((t[0], gene_symbol, "duplicate"))

    # Always keep whitelisted transcripts even without gene_symbol
    for tid, header, seq, length in whitelisted_no_symbol:
        kept.append((tid, header, seq, "", "whitelisted_no_symbol"))

    # Write output
    with open(args.output, "w") as out:
        for tid, header, seq, gene_symbol, reason in kept:
            out.write(f"{header}\n")
            # Write sequence in 60-char lines
            for i in range(0, len(seq), 60):
                out.write(seq[i:i+60] + "\n")

    # Write report if requested
    if args.report:
        with open(args.report, "w") as rpt:
            rpt.write("transcript_id\tgene_symbol\tstatus\treason\n")
            for tid, header, seq, gene_symbol, reason in kept:
                rpt.write(f"{tid}\t{gene_symbol}\tkept\t{reason}\n")
            for tid, gene_symbol, reason in skipped_dup:
                rpt.write(f"{tid}\t{gene_symbol}\tskipped_duplicate\t{reason}\n")
            for tid, header, seq, length in no_symbol:
                rpt.write(f"{tid}\t\tskipped_no_symbol\t\n")

    # Summary
    print(f"Input transcripts: {sum(len(v) for v in by_symbol.values()) + len(no_symbol) + len(whitelisted_no_symbol)}", file=sys.stderr)
    print(f"Whitelisted transcripts: {len(whitelist)}", file=sys.stderr)
    print(f"Unique gene_symbols: {len(by_symbol)}", file=sys.stderr)
    print(f"Skipped (no gene_symbol): {len(no_symbol)}", file=sys.stderr)
    print(f"Skipped (duplicate gene_symbol): {len(skipped_dup)}", file=sys.stderr)
    print(f"Kept: {len(kept)}", file=sys.stderr)
    print(f"  - Whitelisted: {sum(1 for _, _, _, _, r in kept if 'whitelisted' in r)}", file=sys.stderr)
    print(f"  - Longest per gene_symbol: {sum(1 for _, _, _, _, r in kept if r == 'longest')}", file=sys.stderr)


if __name__ == "__main__":
    main()
