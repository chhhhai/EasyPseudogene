#!/usr/bin/env python3
"""
Count disruptive mutations per GeneWise block.
Usage: python3 count_genewise_mutations.py pseudogene-genewise.txt > mutation_counts.tsv
"""
import argparse
import re
import sys


def main():
    ap = argparse.ArgumentParser(description="Count stop codons and frameshifts from GeneWise output")
    ap.add_argument("genewise_file", help="Input GeneWise output file")
    args = ap.parse_args()
    
    try:
        with open(args.genewise_file, 'r') as fh:
            content = fh.read()
    except IOError as e:
        print(f"ERROR: Cannot read {args.genewise_file}: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Split by "//" delimiter
    blocks = content.split("//")
    
    for block in blocks:
        # Extract transcript ID
        id_match = re.search(r'Query protein:\s+(\S+)', block)
        if not id_match:
            continue
        
        transcript_id = id_match.group(1)
        
        # Count frameshifts (!)
        shift_count = len(re.findall(r'!', block))
        
        # Count stop codons (* and :X[)
        stop_count = len(re.findall(r'\*', block))
        stopx_count = len(re.findall(r':X\[', block))
        
        # If no * but has :X[, use :X[ count
        if stop_count == 0 and stopx_count > 0:
            stop_count = stopx_count
        
        print(f"{transcript_id}\t{stop_count}\t{shift_count}")


if __name__ == "__main__":
    main()
