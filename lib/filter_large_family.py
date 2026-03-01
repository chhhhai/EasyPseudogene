#!/usr/bin/env python3
import argparse
import re
import sys


def parse_desc(path):
    if not path:
        return {}
    desc = {}
    try:
        with open(path, "r") as fh:
            for line in fh:
                if not line.strip():
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 2:
                    continue
                key = parts[0]
                desc[key] = "\t".join(parts[1:])
    except FileNotFoundError:
        return {}
    return desc


def parse_patterns(patterns):
    if not patterns:
        return []
    if "|" in patterns:
        parts = [p for p in patterns.split("|") if p]
    else:
        parts = [p for p in patterns.split(",") if p]
    return [re.compile(p, re.IGNORECASE) for p in parts]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--proteins", required=True)
    ap.add_argument("--out-list", required=True)
    ap.add_argument("--out-report", required=True)
    ap.add_argument("--patterns", required=True)
    ap.add_argument("--desc", default="")
    args = ap.parse_args()

    patterns = parse_patterns(args.patterns)
    desc_map = parse_desc(args.desc)

    def match_any(text):
        if not text:
            return None
        for pat in patterns:
            if pat.search(text):
                return pat.pattern
        return None

    out_ids = set()
    with open(args.out_report, "w") as rep:
        rep.write("transcript_id\tgene_id\tgene_symbol\tmatched_pattern\tsource\n")
        with open(args.proteins, "r") as fh:
            for line in fh:
                if not line.startswith(">"):
                    continue
                header = line[1:].strip()
                parts = header.split()
                tid = parts[0]
                gene_id = ""
                gene_symbol = ""
                for tok in parts[1:]:
                    if tok.startswith("gene="):
                        gene_id = tok.split("=", 1)[1]
                    elif tok.startswith("gene_symbol="):
                        gene_symbol = tok.split("=", 1)[1]

                matched = None
                source = ""
                matched = match_any(gene_symbol)
                if matched:
                    source = "gene_symbol"
                else:
                    matched = match_any(gene_id)
                    if matched:
                        source = "gene_id"
                    else:
                        desc = ""
                        if gene_id in desc_map:
                            desc = desc_map.get(gene_id, "")
                        elif gene_symbol in desc_map:
                            desc = desc_map.get(gene_symbol, "")
                        matched = match_any(desc)
                        if matched:
                            source = "description"

                if matched:
                    out_ids.add(tid)
                    rep.write(f"{tid}\t{gene_id}\t{gene_symbol}\t{matched}\t{source}\n")

    with open(args.out_list, "w") as outp:
        for tid in sorted(out_ids):
            outp.write(f"{tid}\n")


if __name__ == "__main__":
    sys.exit(main())
