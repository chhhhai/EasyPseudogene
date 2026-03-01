#!/usr/bin/env python3
import argparse
import re


def clean_target_seq(line):
    # Keep only amino acid letters plus X and ! (frameshift marker)
    return "".join([c for c in line if c.isalpha() or c in ("X", "!")])


def main():
    ap = argparse.ArgumentParser(description="Extract protein alignment from GeneWise text output.")
    ap.add_argument("--input", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    with open(args.input, "r") as fh, open(args.output, "w") as out:
        block = []
        for line in fh:
            block.append(line.rstrip("\n"))
            if line.strip() == "//":
                process_block(block, out)
                block = []
        if block:
            process_block(block, out)


def process_block(lines, out):
    qid = None
    for line in lines:
        if line.startswith("Query protein:"):
            qid = line.split("Query protein:")[1].strip()
            break
    if not qid:
        return
    qid_short = qid.split(".")[0]

    qseq_parts = []
    tseq_parts = []

    # Iterate with index for lookahead
    for i in range(len(lines)):
        line = lines[i]
        if line.startswith(qid_short):
            parts = line.split()
            if len(parts) < 3:
                continue
            qseg = parts[2]
            # target protein line is typically two lines below
            if i + 2 < len(lines):
                tline = lines[i + 2]
                tseg = clean_target_seq(tline)
                qseq_parts.append(qseg)
                tseq_parts.append(tseg)

    if qseq_parts and tseq_parts:
        qseq = "".join(qseq_parts)
        tseq = "".join(tseq_parts)
        out.write(f">{qid}\n{qseq}\n{tseq}\n")


if __name__ == "__main__":
    main()
