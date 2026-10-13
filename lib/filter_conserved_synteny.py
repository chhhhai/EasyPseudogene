#!/usr/bin/env python3
"""Require candidate loci to overlap independently inferred target synteny intervals.

Input TSV columns: query_id, chrom, start0, end0. Coordinates are 0-based,
half-open and may contain multiple intervals per query. The intervals should
come from an independent pairwise genome alignment, not the candidate PAF.
"""
import argparse
import csv
from collections import defaultdict


def overlap_fraction(candidate, expected):
    left = max(candidate[1], expected[1])
    right = min(candidate[2], expected[2])
    overlap = max(0, right - left)
    return (overlap / max(1, candidate[2] - candidate[1]),
            overlap / max(1, expected[2] - expected[1]))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--synteny-tsv", required=True)
    parser.add_argument("--locus-map", required=True)
    parser.add_argument("--in-list", required=True)
    parser.add_argument("--out-pass", required=True)
    parser.add_argument("--out-fail", required=True)
    parser.add_argument("--out-report", required=True)
    parser.add_argument("--min-overlap", type=float, default=0.1)
    args = parser.parse_args()

    expected = defaultdict(list)
    with open(args.synteny_tsv, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            expected[row["query_id"]].append(
                (row["chrom"], int(row["start0"]), int(row["end0"])))
    with open(args.locus_map, newline="") as handle:
        loci = {row["locus_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
    with open(args.in_list) as handle:
        candidates = sorted(set(handle.read().splitlines()))

    with (open(args.out_pass, "w") as passed,
          open(args.out_fail, "w") as failed,
          open(args.out_report, "w", newline="") as report):
        writer = csv.writer(report, delimiter="\t", lineterminator="\n")
        writer.writerow(("locus_id", "query_id", "status", "best_candidate_overlap",
                         "best_expected_overlap"))
        for locus_id in candidates:
            if locus_id not in loci:
                raise ValueError(f"Candidate {locus_id} absent from locus map")
            row = loci[locus_id]
            query = row["query_id"]
            candidate = (row["chrom"], int(row["start0"]), int(row["end0"]))
            overlaps = [overlap_fraction(candidate, interval)
                        for interval in expected[query] if interval[0] == candidate[0]]
            match = [pair for pair in overlaps if min(pair) >= args.min_overlap]
            best = max(overlaps, key=lambda pair: min(pair), default=(0.0, 0.0))
            status = "pass" if match else "no_expected_interval" if not expected[query] else "discordant"
            (passed if status == "pass" else failed).write(locus_id + "\n")
            writer.writerow((locus_id, query, status, f"{best[0]:.4f}", f"{best[1]:.4f}"))


if __name__ == "__main__":
    main()
