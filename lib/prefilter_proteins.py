#!/usr/bin/env python3
"""
Prefilter protein FASTA:
1. Only keep transcripts with a non-empty gene_symbol
2. If multiple transcripts share the same gene_symbol, keep only the longest one
"""
import argparse
import re
import sys


def parse_header(header):
    """Extract transcript_id and gene_symbol from FASTA header."""
    # >ENST00000622028.1 gene=ENSG00000277282.1 ensp=ENSP00000481738.1 gene_symbol=IGHV1OR21-1 cds_len=353 pep_len=117
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


def main():
    ap = argparse.ArgumentParser(description="Prefilter proteins: keep only transcripts with gene_symbol, one per gene_symbol (longest).")
    ap.add_argument("--input", "-i", required=True, help="Input protein FASTA")
    ap.add_argument("--output", "-o", required=True, help="Output filtered FASTA")
    ap.add_argument("--report", "-r", help="Optional TSV report of kept/skipped")
    args = ap.parse_args()

    # First pass: collect all transcripts grouped by gene_symbol
    by_symbol = {}  # gene_symbol -> list of (tid, header, seq, length)
    no_symbol = []  # transcripts without gene_symbol

    with open(args.input, "r") as fh:
        for header, seq in read_fasta(fh):
            tid, gene_symbol = parse_header(header)
            length = len(seq)
            if not gene_symbol or gene_symbol == "":
                no_symbol.append((tid, header, seq, length))
            else:
                by_symbol.setdefault(gene_symbol, []).append((tid, header, seq, length))

    # For each gene_symbol, keep the longest transcript
    kept = []
    skipped_dup = []
    for gene_symbol, transcripts in sorted(by_symbol.items()):
        # Sort by length descending, then by tid for reproducibility
        transcripts.sort(key=lambda x: (-x[3], x[0]))
        best = transcripts[0]
        kept.append((best[0], best[1], best[2], gene_symbol))
        for t in transcripts[1:]:
            skipped_dup.append((t[0], gene_symbol, "duplicate"))

    # Write output
    with open(args.output, "w") as out:
        for tid, header, seq, gene_symbol in kept:
            out.write(f"{header}\n")
            # Write sequence in 60-char lines
            for i in range(0, len(seq), 60):
                out.write(seq[i:i+60] + "\n")

    # Write report if requested
    if args.report:
        with open(args.report, "w") as rpt:
            rpt.write("transcript_id\tgene_symbol\tstatus\n")
            for tid, header, seq, gene_symbol in kept:
                rpt.write(f"{tid}\t{gene_symbol}\tkept\n")
            for tid, gene_symbol, reason in skipped_dup:
                rpt.write(f"{tid}\t{gene_symbol}\tskipped_duplicate\n")
            for tid, header, seq, length in no_symbol:
                rpt.write(f"{tid}\t\tskipped_no_symbol\n")

    # Summary
    print(f"Input transcripts: {sum(len(v) for v in by_symbol.values()) + len(no_symbol)}", file=sys.stderr)
    print(f"Unique gene_symbols: {len(by_symbol)}", file=sys.stderr)
    print(f"Skipped (no gene_symbol): {len(no_symbol)}", file=sys.stderr)
    print(f"Skipped (duplicate gene_symbol): {len(skipped_dup)}", file=sys.stderr)
    print(f"Kept: {len(kept)}", file=sys.stderr)


if __name__ == "__main__":
    main()
