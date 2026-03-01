#!/usr/bin/env python3
import argparse
import re


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


def match_any(text, patterns):
    if not text:
        return None
    for pat in patterns:
        if pat.search(text):
            return pat.pattern
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--proteins", required=True)
    ap.add_argument("--out-list", required=True)
    ap.add_argument("--out-report", required=True)
    ap.add_argument("--patterns", required=True)
    ap.add_argument("--desc", default="")
    ap.add_argument("--biotype", default="")
    ap.add_argument("--exclude-biotypes", default="")
    ap.add_argument("--gtf", default="")
    args = ap.parse_args()

    patterns = parse_patterns(args.patterns)
    desc_map = parse_desc(args.desc)

    exclude_biotypes = set()
    if args.exclude_biotypes:
        for b in args.exclude_biotypes.split(","):
            b = b.strip().lower()
            if b:
                exclude_biotypes.add(b)

    biotype_hits = {}
    if args.biotype:
        try:
            with open(args.biotype, "r") as fh:
                for line in fh:
                    if not line.strip():
                        continue
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) < 4:
                        continue
                    gene_id, transcript_id, gene_bio, tx_bio = parts[:4]
                    gene_bio_l = gene_bio.lower()
                    tx_bio_l = tx_bio.lower()
                    matched = None
                    if gene_bio_l in exclude_biotypes:
                        matched = f"gene_biotype:{gene_bio}"
                    elif tx_bio_l in exclude_biotypes:
                        matched = f"transcript_biotype:{tx_bio}"
                    if matched:
                        biotype_hits[transcript_id] = matched
        except FileNotFoundError:
            pass

    intronless = set()
    if args.gtf:
        try:
            exon_counts = {}
            with open(args.gtf, "r") as fh:
                for line in fh:
                    if not line.strip() or line.startswith("#"):
                        continue
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) < 9:
                        continue
                    feature = parts[2]
                    if feature != "exon":
                        continue
                    attr = parts[8]
                    tid = None
                    if "transcript_id" in attr:
                        import re
                        m = re.search(r'transcript_id \"([^\"]+)\"', attr)
                        if m:
                            tid = m.group(1)
                    if tid:
                        exon_counts[tid] = exon_counts.get(tid, 0) + 1
            for tid, n in exon_counts.items():
                if n <= 1:
                    intronless.add(tid)
        except FileNotFoundError:
            pass

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

                matched = match_any(gene_symbol, patterns)
                source = "gene_symbol" if matched else ""
                if not matched:
                    matched = match_any(gene_id, patterns)
                    if matched:
                        source = "gene_id"
                if not matched:
                    desc = ""
                    if gene_id in desc_map:
                        desc = desc_map.get(gene_id, "")
                    elif gene_symbol in desc_map:
                        desc = desc_map.get(gene_symbol, "")
                    matched = match_any(desc, patterns)
                    if matched:
                        source = "description"

                if matched:
                    out_ids.add(tid)
                    rep.write(f"{tid}\t{gene_id}\t{gene_symbol}\t{matched}\t{source}\n")

                if tid in biotype_hits:
                    out_ids.add(tid)
                    rep.write(f"{tid}\t{gene_id}\t{gene_symbol}\t{biotype_hits[tid]}\tbiotype\n")

                if tid in intronless:
                    out_ids.add(tid)
                    rep.write(f"{tid}\t{gene_id}\t{gene_symbol}\tintronless\tgtf\n")

    with open(args.out_list, "w") as outp:
        for tid in sorted(out_ids):
            outp.write(f"{tid}\n")


if __name__ == "__main__":
    main()
