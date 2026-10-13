#!/usr/bin/env python3
"""Annotate candidate loci with flanking gene order from the two species.

This is a protein-neighborhood proxy, not independent whole-genome synteny.
Unresolved loci are retained by the main pipeline.
"""
import argparse
import bisect
import csv
import gzip
from collections import defaultdict

from annotation_io import Annotation, stable_id


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def reference_genes(path):
    model = Annotation(path)
    genes = {stable_id(gene): value for gene, value in model.gene_locations.items()}
    by_chrom = defaultdict(list)
    for gene, (chrom, start, end, strand) in genes.items():
        by_chrom[chrom].append((start, end, gene, strand))
    for chrom in by_chrom:
        by_chrom[chrom].sort()
    return genes, by_chrom


def annotation_queries(path):
    with open(path, newline="") as handle:
        return {row["query_id"]: row["gene_id"]
                for row in csv.DictReader(handle, delimiter="\t") if row["gene_id"]}


def best_anchors(path, query_genes):
    anchors = {}
    with open(path, newline="") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if len(row) < 17 or row[0] not in query_genes:
                continue
            try:
                qlen, qstart, qend = int(row[8]), int(row[9]), int(row[10])
                alen, matches = int(row[6]), int(row[11])
                coverage = (qend - qstart) / qlen
                identity = matches / alen
            except (ValueError, ZeroDivisionError):
                continue
            if coverage < 0.7 or identity < 0.5:
                continue
            gene = query_genes[row[0]]
            score = (coverage, identity, alen)
            if gene not in anchors or score > anchors[gene][0]:
                anchors[gene] = (score, row[2], int(row[3]), int(row[4]))
    return anchors


def assess(row, genes, by_chrom, anchors, query_genes, neighbors, reference_window, target_window):
    gene_id = query_genes.get(row["query_id"])
    base = {"locus_id": row["locus_id"], "query_id": row["query_id"],
            "reference_gene_id": gene_id or "", "method": "protein_neighborhood_proxy",
            "status": "unresolved", "reason": "reference_gene_unmapped",
            "left_expected": 0, "right_expected": 0, "wrong_side": 0, "same_chrom_anchors": 0}
    if not gene_id or gene_id not in genes:
        return base
    ref_chrom, ref_start, ref_end, ref_strand = genes[gene_id]
    ordered = by_chrom[ref_chrom]
    positions = [entry[0] for entry in ordered]
    index = bisect.bisect_left(positions, ref_start)
    while index < len(ordered) and ordered[index][2] != gene_id:
        index += 1
    if index == len(ordered):
        return base
    try:
        target_start, target_end = int(row["start0"]), int(row["end0"])
    except ValueError:
        return base
    target_chrom = row["chrom"]
    same_orientation = ref_strand == row["strand"]
    counts = {"left": 0, "right": 0, "wrong": 0, "same_chrom": 0}
    for side, indices in (
        ("left", range(index - 1, max(-1, index - neighbors - 1), -1)),
        ("right", range(index + 1, min(len(ordered), index + neighbors + 1))),
    ):
        for j in indices:
            neighbor_start, neighbor_end, neighbor_gene, _ = ordered[j]
            if side == "left" and ref_start - neighbor_end > reference_window:
                break
            if side == "right" and neighbor_start - ref_end > reference_window:
                break
            anchor = anchors.get(neighbor_gene)
            if not anchor or anchor[1] != target_chrom:
                continue
            counts["same_chrom"] += 1
            _, _, astart, aend = anchor
            if aend < target_start and target_start - aend <= target_window:
                target_side = "left"
            elif astart > target_end and astart - target_end <= target_window:
                target_side = "right"
            else:
                continue
            expected = side if same_orientation else ("right" if side == "left" else "left")
            if target_side == expected:
                counts[side] += 1
            else:
                counts["wrong"] += 1
    base.update(left_expected=counts["left"], right_expected=counts["right"],
                wrong_side=counts["wrong"], same_chrom_anchors=counts["same_chrom"])
    if counts["left"] and counts["right"] and counts["left"] + counts["right"] >= 3:
        base.update(status="supported", reason="flanking_anchors_on_both_sides")
    elif counts["wrong"] >= 3 and counts["wrong"] > counts["left"] + counts["right"]:
        base.update(status="discordant", reason="flanking_anchors_wrong_order")
    else:
        base["reason"] = "insufficient_ordered_flanking_anchors"
    return base


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--reference-annotation", "--reference-gtf", dest="reference_gtf", required=True)
    ap.add_argument("--reference-annotations", required=True)
    ap.add_argument("--best-hits", required=True)
    ap.add_argument("--locus-map", required=True)
    ap.add_argument("--candidate-ids", required=True)
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--neighbors", type=int, default=10)
    ap.add_argument("--reference-window", type=int, default=5_000_000)
    ap.add_argument("--target-window", type=int, default=10_000_000)
    args = ap.parse_args()
    genes, by_chrom = reference_genes(args.reference_gtf)
    query_genes = annotation_queries(args.reference_annotations)
    anchors = best_anchors(args.best_hits, query_genes)
    with open(args.candidate_ids) as handle:
        candidates = set(handle.read().splitlines())
    with open(args.locus_map, newline="") as handle:
        loci = [row for row in csv.DictReader(handle, delimiter="\t")
                if row["locus_id"] in candidates]
    fields = ("locus_id", "query_id", "reference_gene_id", "method", "status", "reason",
              "left_expected", "right_expected", "wrong_side", "same_chrom_anchors")
    with open(args.out_tsv, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in loci:
            writer.writerow(assess(row, genes, by_chrom, anchors, query_genes,
                                   args.neighbors, args.reference_window, args.target_window))
    print(f"Neighborhood evidence: {len(loci)} candidate loci; output {args.out_tsv}")


if __name__ == "__main__":
    main()
