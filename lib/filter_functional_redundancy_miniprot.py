#!/usr/bin/env python3
import argparse
import csv


def interval_coverage(intervals, query_length):
    merged = []
    for start, end in sorted(intervals):
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return sum(end - start for start, end in merged) / query_length


def same_locus(copy, candidate):
    if copy["chrom"] != candidate["chrom"]:
        return False
    overlap = max(0, min(copy["end"], candidate["end0"])
                  - max(copy["start"], candidate["start0"]))
    if not overlap:
        return False
    copy_len = max(1, copy["end"] - copy["start"])
    candidate_len = max(1, candidate["end0"] - candidate["start0"])
    return overlap / min(copy_len, candidate_len) >= 0.5


def load_loci(path):
    loci = {}
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            loci[row["locus_id"]] = {"query": row["query_id"],
                                      "chrom": row["chrom"],
                                      "start0": int(row["start0"]),
                                      "end0": int(row["end0"])}
    return loci


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--paf", required=True)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--min-aln-len", type=int, default=30)
    ap.add_argument("--min-cov", type=float, default=0.05)
    ap.add_argument("--min-id", type=float, default=0.3)
    ap.add_argument("--min-mapq", type=int, default=0)
    ap.add_argument("--copy-cov", type=float, default=0.8)
    ap.add_argument("--copy-min-id", type=float, default=0.5)
    ap.add_argument("--copy-policy", choices=("intact", "potential"), default="intact",
                    help="Whether an additional >=copy-cov hit must appear undisrupted")
    ap.add_argument("--target-gap", type=int, default=10000)
    ap.add_argument("--locus-map", help="Locus map; exclude the candidate locus itself from full-copy counts")
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
            tags = {}
            for item in parts[12:]:
                tag = item.split(":", 2)
                if len(tag) == 3 and tag[0] in {"fs", "st"} and tag[1] == "i":
                    tags[tag[0]] = int(tag[2])
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
                (tname, tstart, tend, qstart, qend, qcov, ident,
                 tags.get("fs"), tags.get("st"))
            )

    qc_path = args.out_prefix + "_qc.tsv"
    redundant_path = args.out_prefix + "_redundant.txt"
    pass_path = args.out_prefix + "_pass.txt"
    copies_path = args.out_prefix + "_copies.tsv"

    loci = load_loci(args.locus_map) if args.locus_map else {
        query: None for query in hits_by_query
    }
    copies_by_query = {}
    for qname, hits in hits_by_query.items():
        qlen = qlen_by_query.get(qname, 0)
        if qlen <= 0:
            continue
        by_tname = {}
        for h in hits:
            by_tname.setdefault(h[0], []).append(h)
        copies = []
        for tname, thits in by_tname.items():
            thits.sort(key=lambda x: (x[1], x[2]))
            # Cluster nearby hits into target copies.
            clusters = []
            current = []
            cur_end = None
            for h in thits:
                ts, te = h[1], h[2]
                if cur_end is None or ts <= cur_end + args.target_gap:
                    current.append(h)
                    cur_end = max(cur_end or te, te)
                else:
                    clusters.append(current)
                    current = [h]
                    cur_end = te
            if current:
                clusters.append(current)
            for cluster_hits in clusters:
                aligned = sum(h[4] - h[3] for h in cluster_hits)
                copies.append({"chrom": tname,
                               "start": min(h[1] for h in cluster_hits),
                               "end": max(h[2] for h in cluster_hits),
                               "coverage": interval_coverage(
                                   ((h[3], h[4]) for h in cluster_hits), qlen),
                               "identity": (sum(h[6] * (h[4] - h[3]) for h in cluster_hits)
                                            / max(1, aligned)),
                               "frameshifts": sum(h[7] or 0 for h in cluster_hits),
                               "stops": sum(h[8] or 0 for h in cluster_hits),
                               "tags_complete": all(h[7] is not None and h[8] is not None
                                                    for h in cluster_hits),
                               "undisrupted": all(h[7] == 0 and h[8] == 0
                                                  for h in cluster_hits)})
        copies_by_query[qname] = copies

    with open(qc_path, "w") as qc, open(redundant_path, "w") as r, \
            open(pass_path, "w") as p, open(copies_path, "w", newline="") as detail:
        qc.write("locus_id\tquery\tbest_cov\tbest_id\tfull_copy_count\tother_potential_copy_count\tother_full_copy_count\tredundant\tcopy_policy\n")
        writer = csv.writer(detail, delimiter="\t", lineterminator="\n")
        writer.writerow(("locus_id", "query_id", "chrom", "start0", "end0",
                         "query_coverage", "identity", "frameshifts", "stops",
                         "mutation_tags_complete", "candidate_locus",
                         "potential_full_copy", "undisrupted_full_copy_proxy"))
        for locus_id, candidate in sorted(loci.items()):
            qname = candidate["query"] if candidate else locus_id
            copies = copies_by_query.get(qname, [])
            for copy in copies:
                same = bool(candidate and same_locus(copy, candidate))
                potential = copy["coverage"] >= args.copy_cov and not same
                intact_proxy = (potential and copy["identity"] >= args.copy_min_id
                                and copy["undisrupted"])
                writer.writerow((locus_id, qname, copy["chrom"], copy["start"],
                                 copy["end"], f'{copy["coverage"]:.4f}',
                                 f'{copy["identity"]:.4f}', copy["frameshifts"],
                                 copy["stops"], int(copy["tags_complete"]),
                                 int(same), int(potential), int(intact_proxy)))
            best_cov = max((copy["coverage"] for copy in copies), default=0)
            best_id = max((copy["identity"] for copy in copies), default=0)
            full_copy = sum(copy["coverage"] >= args.copy_cov for copy in copies)
            if candidate:
                other_potential_copy = sum(copy["coverage"] >= args.copy_cov and
                                           not same_locus(copy, candidate) for copy in copies)
                other_full_copy = sum(copy["coverage"] >= args.copy_cov and
                                      copy["identity"] >= args.copy_min_id and
                                      copy["undisrupted"] and
                                      not same_locus(copy, candidate) for copy in copies)
                redundant = bool(other_potential_copy if args.copy_policy == "potential"
                                 else other_full_copy)
            else:
                # Backward-compatible query-level mode for callers without a locus map.
                other_potential_copy = max(0, full_copy - 1)
                other_full_copy = max(0, full_copy - 1)
                redundant = bool(other_full_copy)
            qc.write(f"{locus_id}\t{qname}\t{best_cov:.4f}\t{best_id:.4f}\t{full_copy}\t{other_potential_copy}\t{other_full_copy}\t{int(redundant)}\t{args.copy_policy}\n")
            if redundant:
                r.write(f"{locus_id}\n")
            else:
                p.write(f"{locus_id}\n")


if __name__ == "__main__":
    main()
