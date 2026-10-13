#!/usr/bin/env python3
"""Flag GeneWise disruptions contradicted by a strong intact miniprot locus hit.

The two aligners inspect the same assembly. This is an alignment conflict
filter, not independent validation of a biological genotype.
"""
import argparse
import csv
from collections import defaultdict


def load_mutations(path):
    with open(path, newline="") as handle:
        return {row[0]: (int(row[1]), int(row[2]))
                for row in csv.reader(handle, delimiter="\t") if len(row) >= 3}


def load_loci(path):
    with open(path, newline="") as handle:
        return {row["locus_id"]: row for row in csv.DictReader(handle, delimiter="\t")}


def overlap_fraction(hit, locus):
    if hit["chrom"] != locus["chrom"]:
        return 0.0
    start, end = int(locus["start0"]), int(locus["end0"])
    overlap = max(0, min(hit["end"], end) - max(hit["start"], start))
    return overlap / max(1, min(hit["end"] - hit["start"], end - start))


def load_hits(path):
    result = defaultdict(list)
    with open(path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            try:
                qlen, qstart, qend = map(int, fields[1:4])
                start, end = int(fields[7]), int(fields[8])
                matches, aligned = int(fields[9]), int(fields[10])
            except ValueError:
                continue
            if qlen <= 0 or aligned <= 0:
                continue
            tags = {}
            for item in fields[12:]:
                parts = item.split(":", 2)
                if len(parts) == 3 and parts[0] in {"st", "fs"} and parts[1] == "i":
                    try:
                        tags[parts[0]] = int(parts[2])
                    except ValueError:
                        pass
            result[fields[0]].append({"chrom": fields[5], "start": start,
                                      "end": end, "coverage": (qend - qstart) / qlen,
                                      "identity": matches / aligned,
                                      "aligned": aligned, "stops": tags.get("st"),
                                      "frameshifts": tags.get("fs")})
    return result


def classify(locus, gw, hits, min_cov, min_id):
    nearby = [hit for hit in hits if overlap_fraction(hit, locus) >= 0.5]
    if not nearby:
        return "not_assessed", None
    hit = max(nearby, key=lambda row: (row["aligned"], row["identity"]))
    if hit["stops"] is None or hit["frameshifts"] is None:
        return "missing_miniprot_mutation_tags", hit
    if (gw[0] and hit["stops"]) or (gw[1] and hit["frameshifts"]):
        return "same_type_disruption", hit
    if hit["stops"] or hit["frameshifts"]:
        return "different_disruption_type_review", hit
    if hit["coverage"] >= min_cov and hit["identity"] >= min_id:
        return "intact_alignment_conflict", hit
    return "intact_but_low_alignment_support", hit


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--screening-paf", required=True)
    ap.add_argument("--locus-map", required=True)
    ap.add_argument("--mutation-counts", required=True)
    ap.add_argument("--in-list", required=True)
    ap.add_argument("--out-exclude", required=True)
    ap.add_argument("--out-report", required=True)
    ap.add_argument("--min-cov", type=float, default=0.8)
    ap.add_argument("--min-id", type=float, default=0.5)
    args = ap.parse_args()
    if not 0 <= args.min_cov <= 1 or not 0 <= args.min_id <= 1:
        ap.error("coverage and identity thresholds must be in [0, 1]")
    loci = load_loci(args.locus_map)
    mutations = load_mutations(args.mutation_counts)
    hits = load_hits(args.screening_paf)
    with open(args.in_list) as handle:
        candidates = sorted(set(handle.read().splitlines()))
    fields = ("locus_id", "query_id", "geneWise_stops", "geneWise_frameshifts",
              "miniprot_stops", "miniprot_frameshifts", "miniprot_coverage",
              "miniprot_identity", "status", "excluded")
    with open(args.out_exclude, "w") as excluded, open(args.out_report, "w", newline="") as report:
        writer = csv.DictWriter(report, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for locus_id in candidates:
            locus = loci[locus_id]
            query_id = locus["query_id"]
            gw = mutations.get(locus_id, (0, 0))
            status, hit = classify(locus, gw, hits.get(query_id, []), args.min_cov, args.min_id)
            reject = status == "intact_alignment_conflict"
            if reject:
                excluded.write(locus_id + "\n")
            writer.writerow({"locus_id": locus_id, "query_id": query_id,
                             "geneWise_stops": gw[0], "geneWise_frameshifts": gw[1],
                             "miniprot_stops": "" if hit is None else hit["stops"],
                             "miniprot_frameshifts": "" if hit is None else hit["frameshifts"],
                             "miniprot_coverage": "" if hit is None else f'{hit["coverage"]:.4f}',
                             "miniprot_identity": "" if hit is None else f'{hit["identity"]:.4f}',
                             "status": status, "excluded": int(reject)})


if __name__ == "__main__":
    main()
