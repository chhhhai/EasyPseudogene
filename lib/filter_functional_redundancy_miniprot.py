#!/usr/bin/env python3
import argparse


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--paf", required=True)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--min-aln-len", type=int, default=30)
    ap.add_argument("--min-cov", type=float, default=0.3)
    ap.add_argument("--min-id", type=float, default=0.3)
    ap.add_argument("--min-mapq", type=int, default=0)
    ap.add_argument("--copy-cov", type=float, default=0.8)
    ap.add_argument("--target-gap", type=int, default=10000)
    args = ap.parse_args()

    # Assemble hits into locus "copies" on target and compute query coverage per copy.
    hits_by_query = {}
    qlen_by_query = {}

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
            mapq = int(parts[11])
            if alen <= 0 or qlen <= 0:
                continue
            if alen < args.min_aln_len:
                continue
            if mapq < args.min_mapq:
                continue

            qcov = (qend - qstart) / qlen
            ident = nmatch / alen
            if qcov < args.min_cov or ident < args.min_id:
                continue

            qlen_by_query[qname] = qlen
            hits_by_query.setdefault(qname, []).append(
                (tname, tstart, tend, qstart, qend, qcov, ident)
            )

    qc_path = args.out_prefix + "_qc.tsv"
    redundant_path = args.out_prefix + "_redundant.txt"
    pass_path = args.out_prefix + "_pass.txt"

    with open(qc_path, "w") as qc, open(redundant_path, "w") as r, open(pass_path, "w") as p:
        qc.write("query\tbest_cov\tbest_id\tfull_copy_count\tredundant\n")
        for qname in sorted(hits_by_query.keys()):
            hits = hits_by_query[qname]
            qlen = qlen_by_query.get(qname, 0)
            if qlen <= 0:
                continue

            best_cov = 0.0
            best_id = 0.0
            full_copy = 0

            # Group by target chromosome/scaffold
            by_tname = {}
            for h in hits:
                by_tname.setdefault(h[0], []).append(h)

            for tname, thits in by_tname.items():
                thits.sort(key=lambda x: (x[1], x[2]))
                # Cluster into loci by target gap
                loci = []
                current = []
                cur_end = None
                for h in thits:
                    ts, te = h[1], h[2]
                    if cur_end is None or ts <= cur_end + args.target_gap:
                        current.append(h)
                        cur_end = max(cur_end or te, te)
                    else:
                        loci.append(current)
                        current = [h]
                        cur_end = te
                if current:
                    loci.append(current)

                for locus_hits in loci:
                    # Merge query intervals
                    intervals = sorted((h[3], h[4]) for h in locus_hits)
                    merged = []
                    for s, e in intervals:
                        if not merged or s > merged[-1][1]:
                            merged.append([s, e])
                        else:
                            merged[-1][1] = max(merged[-1][1], e)
                    cov = sum(e - s for s, e in merged) / qlen

                    if cov > best_cov:
                        best_cov = cov
                    # approximate identity using best hit in locus
                    best_locus_id = max(h[6] for h in locus_hits)
                    if best_locus_id > best_id:
                        best_id = best_locus_id

                    if cov >= args.copy_cov:
                        full_copy += 1

            redundant = 1 if full_copy >= 2 else 0
            qc.write(f"{qname}\t{best_cov:.4f}\t{best_id:.4f}\t{full_copy}\t{redundant}\n")
            if redundant:
                r.write(f"{qname}\n")
            else:
                p.write(f"{qname}\n")


if __name__ == "__main__":
    main()
