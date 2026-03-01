#!/usr/bin/env python3
"""
Identify single frameshift mutations at splicing/ending regions.
Outputs:
  - splicing-shift.txt: frameshift mutations at splicing sites
  - ending-shift.txt: frameshift mutations at ending regions (< 50 aa remaining)
"""
import argparse
import re
import sys


def main():
    ap = argparse.ArgumentParser(description="Identify single frameshift mutations at splicing/ending regions")
    ap.add_argument("--one-mutation-list", required=True, help="File with transcript IDs having exactly one mutation")
    ap.add_argument("--genewise-file", required=True, help="GeneWise output file")
    ap.add_argument("--out-splicing", default="splicing-shift.txt", help="Output file for splicing frameshift mutations")
    ap.add_argument("--out-ending", default="ending-shift.txt", help="Output file for ending frameshift mutations")
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
    
    # Patterns for frameshift at splicing site
    # Patterns: ![A-Z]   ,   [A-Z]!   ,   !    ,    !
    splicing_patterns = [
        r'![A-Z]\s{3}',      # !A   (with 3 spaces)
        r'\s{3}[A-Z]!',      #   A! (with 3 spaces before)
        r'!\s{3}',           # !   (with 3 spaces)
        r'\s{3}!',           #    ! (with 3 spaces before)
    ]
    
    for block in blocks:
        # Extract transcript ID
        id_match = re.search(r'Query protein:\s+(\S+)', block)
        if not id_match:
            continue
        
        transcript_id = id_match.group(1)
        
        # Skip if not in one-mutation list
        if transcript_id not in one_mutation_set:
            continue
        
        # Check for frameshift at splicing site
        has_splicing_shift = False
        for pattern in splicing_patterns:
            if re.search(pattern, block):
                has_splicing_shift = True
                break
        
        if has_splicing_shift:
            print(transcript_id)
            out_splicing.write(f"{transcript_id}\n")
        
        # Check for frameshift at ending region (pattern: .sp.tr with ! and < 50 aa remaining)
        if re.search(r'\.sp\.tr', block):
            # Extract ID from .sp.tr line
            sp_tr_match = re.search(r'>(\S+?)(?:\.sp\.tr)?', block)
            if sp_tr_match:
                block_id = sp_tr_match.group(1)
                if block_id == transcript_id:
                    # Remove whitespace, newlines, carriage returns
                    block_clean = re.sub(r'\s+', '', block)
                    block_clean = re.sub(r'\n', '', block_clean)
                    block_clean = re.sub(r'\r', '', block_clean)
                    
                    # Extract sequence after "tr"
                    tr_match = re.search(r'(.*)tr(\S+)$', block_clean)
                    if tr_match:
                        remaining_seq = tr_match.group(2)
                        remaining_len = len(remaining_seq)
                        
                        # If remaining length < 50, it's an ending frameshift
                        if remaining_len < 50:
                            out_ending.write(f"{transcript_id}\n")
    
    out_splicing.close()
    out_ending.close()


if __name__ == "__main__":
    main()
