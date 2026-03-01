#!/usr/bin/env python3
import argparse
import sys
from collections import defaultdict


def non_overlapping_hits(hits, overlap_frac):
    # hits: list of (tname, tstart, tend)
    if not hits:
        return 0
    hits.sort(key=lambda x: (x[0], x[1], x[2]))
    count = 0
    last = None
    for tname, s, e in hits:
        if last is None or tname != last[0]:
            count += 1
            last = [tname, s, e]
            continue
        last_s, last_e = last[1], last[2]
        overlap = max(0, min(e, last_e) - max(s, last_s))
        min_len = max(1, min(e - s, last_e - last_s))
        if overlap / min_len < overlap_frac:
            count += 1
            last = [tname, s, e]
        else:
            if e > last_e:
                last[2] = e
    return count


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--paf", required=True)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--min-aln-len", type=int, default=30)
    ap.add_argument("--min-cov", type=float, default=0.3)
    ap.add_argument("--min-id", type=float, default=0.3)
    ap.add_argument("--max-hits", type=int, default=1)
    ap.add_argument("--overlap-frac", type=float, default=0.5)
    args = ap.parse_args()

    best = {}
    hits = defaultdict(list)

    with open(args.paf, "r") as fh:
        for line in fh:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 12:
                continue
            qname = parts[0]
            qlen = int(parts[1])
            qstart = int(parts[2])
            qend = int(parts[3])
            tname = parts[5]
            tstart = int(parts[7])
            tend = int(parts[8])
            nmatch = int(parts[9])
            alen = int(parts[10])
            if alen <= 0 or qlen <= 0:
                continue
            qcov = (qend - qstart) / qlen
            ident = nmatch / alen

            if alen < args.min_aln_len:
                continue
            if qcov < args.min_cov or ident < args.min_id:
                continue

            hits[qname].append((tname, tstart, tend))

            score = nmatch
            prev = best.get(qname)
            if prev is None or score > prev["score"]:
                best[qname] = {
                    "score": score,
                    "cov": qcov,
                    "id": ident,
                    "alen": alen,
                    "tname": tname,
                    "tstart": tstart,
                    "tend": tend,
                }

    qc_path = args.out_prefix + "_qc.tsv"
    pass_path = args.out_prefix + "_pass.txt"

    with open(qc_path, "w") as qc, open(pass_path, "w") as outp:
        qc.write("query\tbest_cov\tbest_id\tbest_alen\thit_count\tbest_tname\tbest_tstart\tbest_tend\n")
        for qname in sorted(set(list(hits.keys()) + list(best.keys()))):
            b = best.get(qname)
            hit_count = non_overlapping_hits(hits.get(qname, []), args.overlap_frac)
            if b is None:
                qc.write(f"{qname}\t0\t0\t0\t{hit_count}\tNA\tNA\tNA\n")
                continue
            qc.write(
                f"{qname}\t{b['cov']:.4f}\t{b['id']:.4f}\t{b['alen']}\t{hit_count}\t"
                f"{b['tname']}\t{b['tstart']}\t{b['tend']}\n"
            )
            if hit_count >= 1 and hit_count <= args.max_hits:
                outp.write(f"{qname}\n")


if __name__ == "__main__":
    sys.exit(main())
