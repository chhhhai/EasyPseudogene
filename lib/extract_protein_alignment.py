#!/usr/bin/env python3
import argparse
from genewise_alignment import segments


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
    parsed = list(segments("\n".join(lines)))
    if parsed:
        qseq = "".join(q for q, _ in parsed)
        tseq = "".join(t for _, t in parsed)
        out.write(f">{qid}\n{qseq}\n{tseq}\n")


if __name__ == "__main__":
    main()
