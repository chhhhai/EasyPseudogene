#!/usr/bin/env python3
import argparse
from alignment_metrics import mutation_query_fraction


def parse_alignment(path):
    records = {}
    with open(path, "r") as fh:
        cur_id = None
        seqs = []
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if cur_id and len(seqs) >= 2:
                    records[cur_id] = (seqs[0], seqs[1])
                cur_id = line[1:].split()[0]
                seqs = []
            else:
                seqs.append(line)
        if cur_id and len(seqs) >= 2:
            records[cur_id] = (seqs[0], seqs[1])
    return records


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--protein-alignment", required=True)
    ap.add_argument("--out-list", required=True)
    ap.add_argument("--out-report", required=True)
    ap.add_argument("--min-frac", type=float, default=0.9)
    args = ap.parse_args()

    records = parse_alignment(args.protein_alignment)
    with open(args.out_list, "w") as outp, open(args.out_report, "w") as rep:
        rep.write("query\tmut_pos\tquery_len\tfraction\n")
        for qid, (qseq, tseq) in records.items():
            pos_x = tseq.find("X")
            pos_bang = tseq.find("!")
            positions = [p for p in (pos_x, pos_bang) if p >= 0]
            if not positions:
                continue
            pos = min(positions)
            frac, covered, qlen = mutation_query_fraction(qseq, tseq, pos)
            if not qlen:
                continue
            if frac >= args.min_frac:
                outp.write(f"{qid}\n")
                rep.write(f"{qid}\t{covered}\t{qlen}\t{frac:.4f}\n")


if __name__ == "__main__":
    main()
