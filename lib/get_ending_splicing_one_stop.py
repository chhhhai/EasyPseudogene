#!/usr/bin/env python3
"""
Identify single stop mutations at splicing/ending regions.
Outputs:
  - splicing-stop.txt: stop mutations at splicing sites
  - ending-stop.txt: stop mutations at ending regions (< 50 aa remaining)
"""
import argparse
import re
import sys
from genewise_alignment import segments


def main():
    ap = argparse.ArgumentParser(description="Identify single stop mutations at splicing/ending regions")
    ap.add_argument("--one-mutation-list", required=True, help="File with transcript IDs having exactly one mutation")
    ap.add_argument("--genewise-file", required=True, help="GeneWise output file")
    ap.add_argument("--out-splicing", default="splicing-stop.txt", help="Output file for splicing stop mutations")
    ap.add_argument("--out-ending", default="ending-stop.txt", help="Output file for ending stop mutations")
    args = ap.parse_args()
    
    # Load one-mutation list
    one_mutation_set = set()
    try:
        with open(args.one_mutation_list, 'r') as fh:
            for line in fh:
                line = line.strip()
                if line:
                    # Extract first field (transcript ID)
                    tid = line.split()[0] if line.split() else line
                    one_mutation_set.add(tid)
    except IOError as e:
        print(f"ERROR: Cannot read {args.one_mutation_list}: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Open output files
    try:
        out_splicing = open(args.out_splicing, 'w')
        out_ending = open(args.out_ending, 'w')
    except IOError as e:
        print(f"ERROR: Cannot write output files: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Process GeneWise output
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
        
        # Skip if not in one-mutation list
        if transcript_id not in one_mutation_set:
            continue
        
        parsed = list(segments(block))
        if parsed:
            query = "".join(q for q, _ in parsed)
            target = "".join(t for _, t in parsed)
            if ':X[' in target:
                out_splicing.write(f"{transcript_id}\n")
            stops = [i for i, aa in enumerate(target) if aa == 'X']
            if any(sum(c.isalpha() for c in query[i + 1:]) < 50 for i in stops):
                out_ending.write(f"{transcript_id}\n")
        elif re.search(r':X\[', block):
            # Legacy compact alignment representation.
            out_splicing.write(f"{transcript_id}\n")

    out_splicing.close()
    out_ending.close()


if __name__ == "__main__":
    main()
