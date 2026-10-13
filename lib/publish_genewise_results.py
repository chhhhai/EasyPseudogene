#!/usr/bin/env python3
"""Publish a readable GeneWise result view without moving legacy work files."""
import argparse
import csv
import hashlib
import json
import os
import shutil
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

from prefilter_proteins import read_fasta
from result_evidence import EVIDENCE_FIELDS, EVIDENCE_SEMANTICS, genewise_evidence
from run_layout import cache_file, work_directory
from publish_result_assets import publish_assets


ROOT = Path(__file__).resolve().parents[1]
FILTERS = (
    "putative_pseudogene_ids.raw",
    "putative_pseudogene_ids",
    "unitary_pseudogene_ids",
    "high_confidence_pseudogene_ids",
    "large_family_ids",
    "predicted_annotation_ids",
    "miniprot_redundancy_redundant",
    "miniprot_intact_conflict",
    "conserved_synteny_failed",
    "false_positive_ids",
    "splicing-stop",
    "ending-stop",
    "low-identity-stop",
    "splicing-shift",
    "ending-shift",
    "low-identity-shift",
    "end_mutation_90pct",
    "multi_disruption_low_quality",
    "high_confidence_pseudogenes_removed",
)
SINGLE_DISRUPTION_FILTERS = (
    "splicing-stop", "ending-stop", "low-identity-stop",
    "splicing-shift", "ending-shift", "low-identity-shift",
    "end_mutation_90pct",
)
LOCUS_FIELDS = (
    "locus_id", "query_id", "chrom", "start0", "end0", "strand",
    "genewise_status", "genewise_stops", "genewise_frameshifts", "stage", "decision_reason",
    "decision_detail", "neighborhood_synteny_status", "neighborhood_synteny_reason",
) + EVIDENCE_FIELDS
QUERY_FIELDS = (
    "query_id", "selected_loci", "disruption_loci", "putative_loci",
    "family_copy_pass_loci", "final_loci", "status",
)


def read_ids(path):
    return set(path.read_text().splitlines()) if path.exists() else set()


def write_tsv(path, fields, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".partial")
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t",
                                lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".partial")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def copy_file(source, destination):
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(destination.name + ".partial")
    shutil.copyfile(source, temporary)
    temporary.replace(destination)


def classify(locus, query, filters, budget_excluded, stops, frameshifts):
    if locus in budget_excluded or query in budget_excluded:
        return "not_run", "compute_budget_excluded", ""
    if locus in filters["high_confidence_pseudogene_ids"]:
        return "final_candidate", "passed_all_enabled_filters", ""
    if locus in filters["unitary_pseudogene_ids"]:
        if locus in filters["miniprot_intact_conflict"]:
            return "filtered", "intact_alignment_conflict", "miniprot_selected_locus"
        if locus in filters["false_positive_ids"]:
            details = [name for name in SINGLE_DISRUPTION_FILTERS
                       if locus in filters[name]]
            return "filtered", "single_disruption_filter", ";".join(details)
        if locus in filters["multi_disruption_low_quality"]:
            return "filtered", "multi_disruption_alignment_quality", ""
        if locus in filters["high_confidence_pseudogenes_removed"]:
            return "filtered", "global_alignment_quality", ""
        return "filtered", "other_final_filter", ""
    if locus in filters["putative_pseudogene_ids"]:
        if locus in filters["large_family_ids"]:
            return "putative", "large_family_filter", ""
        if locus in filters["predicted_annotation_ids"]:
            return "putative", "predicted_reference_annotation_filter", ""
        if locus in filters["miniprot_redundancy_redundant"]:
            return ("putative", "other_full_length_homolog_filter",
                    "copy_is_potential_not_verified_functional")
        if locus in filters["conserved_synteny_failed"]:
            return "putative", "conserved_synteny_filter", ""
        return "putative", "other_family_copy_filter", ""
    if locus in filters["putative_pseudogene_ids.raw"] or stops + frameshifts:
        return "disruption_detected", "protein_prefilter", ""
    return "aligned", "no_genewise_disruption", ""


def publish(run_dir, proteins, species, analysis_version=None, budget_file=None,
            provenance_file=None, comparison_mode="cross"):
    if comparison_mode not in {"cross", "self"}:
        raise ValueError("comparison_mode must be cross or self")
    run_dir = Path(run_dir)
    work = work_directory(run_dir)
    filter_dir = work / "post_genewise_filter"
    map_path = work / "regions" / "locus_map.tsv"
    required = (map_path, filter_dir / "mutation_counts.tsv",
                filter_dir / "high_confidence_pseudogene_ids.txt",
                filter_dir / "effective_config.txt",
                filter_dir / "high_confidence_loci.bed",
                filter_dir / "high_confidence_queries.txt")
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise FileNotFoundError("Cannot publish incomplete GeneWise run: " + ", ".join(missing))
    filters = {name: read_ids(filter_dir / (name + ".txt")) for name in FILTERS}
    budget = Path(budget_file) if budget_file else work / "excluded_due_compute_budget.tsv"
    budget_excluded = set()
    if budget.is_file():
        with budget.open() as handle:
            budget_excluded = {row[0] for row in csv.reader(handle, delimiter="\t")
                               if row and row[0] not in {"query_id", "protein_id", "locus_id"}}
    with (filter_dir / "mutation_counts.tsv").open() as handle:
        mutations = {row[0]: (int(row[1]), int(row[2]))
                     for row in csv.reader(handle, delimiter="\t") if len(row) >= 3}
    with map_path.open() as handle:
        mapped = list(csv.DictReader(handle, delimiter="\t"))
    neighborhood_path = filter_dir / "neighborhood_synteny.tsv"
    neighborhood = {}
    if neighborhood_path.is_file():
        with neighborhood_path.open() as handle:
            neighborhood = {row["locus_id"]: row
                            for row in csv.DictReader(handle, delimiter="\t")}
    locus_ids = {row["locus_id"] for row in mapped}
    if len(locus_ids) != len(mapped):
        raise ValueError("Duplicate locus_id in regions/locus_map.tsv")
    # A compute cap may omit loci from the GeneWise job map. Recover their
    # miniprot coordinates for reporting, while clearly marking them not run.
    mapped_budget_queries = {row["query_id"] for row in mapped}
    best_dir = work / "miniprot"
    best = best_dir / f"{species}.best_hits.uncapped.tsv"
    if not best.is_file():
        best = best_dir / f"{species}.best_hits.tsv"
    if budget_excluded and best.is_file():
        with best.open() as handle:
            for fields in csv.reader(handle, delimiter="\t"):
                if len(fields) < 5 or fields[0] not in budget_excluded:
                    continue
                query, strand, chrom, start0, end0 = fields[:5]
                if query in mapped_budget_queries:
                    continue
                digest = hashlib.sha256("\t".join(fields[:5]).encode()).hexdigest()[:20]
                locus = query if query not in locus_ids else "EPlocus_budget_" + digest
                if locus in locus_ids:
                    continue
                mapped.append({"locus_id": locus, "query_id": query,
                               "strand": strand, "chrom": chrom,
                               "start0": start0, "end0": end0})
                locus_ids.add(locus)
    for stage in FILTERS[:4]:
        unknown = filters[stage] - locus_ids
        if unknown:
            raise ValueError(f"{stage} has {len(unknown)} loci absent from locus_map.tsv")
    with Path(proteins).open() as handle:
        protein_lengths = {header[1:].split()[0]: len(sequence) for header, sequence in read_fasta(handle)}
    with Path(proteins).open() as handle:
        input_queries = [header[1:].split()[0] for header, _ in read_fasta(handle)]
    input_query_ids = set(input_queries)
    if len(input_query_ids) != len(input_queries):
        raise ValueError("Duplicate query ID in input protein FASTA")

    loci = []
    by_query = defaultdict(list)
    for row in mapped:
        locus = row["locus_id"]
        query = row["query_id"]
        if query not in input_query_ids:
            raise ValueError(f"Locus {locus} references absent input query {query}")
        stops, frameshifts = mutations.get(locus, (0, 0))
        stage, reason, detail = classify(locus, query, filters, budget_excluded,
                                         stops, frameshifts)
        output = {**row, "genewise_stops": stops,
                  "genewise_frameshifts": frameshifts,
                  "genewise_status": "not_run" if stage == "not_run" else "completed",
                  "stage": stage,
                  "decision_reason": reason, "decision_detail": detail,
                  "neighborhood_synteny_status": neighborhood.get(locus, {}).get("status", "not_assessed"),
                  "neighborhood_synteny_reason": neighborhood.get(locus, {}).get("reason", "")}
        output.update(genewise_evidence(work / "genewise" / (locus + ".gw.txt"),
                                        protein_lengths.get(query, 0), stops, frameshifts, stage != "not_run"))
        output["evidence_query_id"] = query if stage != "not_run" else ""
        loci.append(output)
        by_query[query].append(output)

    queries = []
    for query in input_queries:
        rows = by_query[query]
        counts = {
            "selected_loci": len(rows),
            "disruption_loci": sum(int(r["genewise_stops"]) + int(r["genewise_frameshifts"]) > 0 for r in rows),
            "putative_loci": sum(r["locus_id"] in filters["putative_pseudogene_ids"] for r in rows),
            "family_copy_pass_loci": sum(r["locus_id"] in filters["unitary_pseudogene_ids"] for r in rows),
            "final_loci": sum(r["stage"] == "final_candidate" for r in rows),
        }
        status = ("final_candidate" if counts["final_loci"] else
                  "filtered" if counts["putative_loci"] else
                  "disruption_detected" if counts["disruption_loci"] else
                  "not_run_compute_budget" if any(r["stage"] == "not_run" for r in rows) else
                  "aligned_no_disruption" if rows else
                  "not_run_compute_budget" if query in budget_excluded else
                  "no_selected_locus")
        queries.append({"query_id": query, **counts, "status": status})

    results = run_dir / "results"
    qc = run_dir / "qc"
    metadata = run_dir / "metadata"
    write_tsv(results / "loci.tsv", LOCUS_FIELDS, loci)
    write_tsv(results / "queries.tsv", QUERY_FIELDS, queries)
    write_tsv(results / "final_loci.tsv", LOCUS_FIELDS,
              (row for row in loci if row["stage"] == "final_candidate"))
    write_tsv(qc / "nonfinal_loci.tsv", LOCUS_FIELDS,
              (row for row in loci if row["stage"] != "final_candidate"))
    stage_counts = [
        ("input", len(input_queries), len(input_queries)),
        ("selected_alignment", len(loci), sum(bool(by_query[q]) for q in input_queries)),
        ("genewise_completed", sum(r["genewise_status"] == "completed" for r in loci),
         sum(any(r["genewise_status"] == "completed" for r in by_query[q]) for q in input_queries)),
        ("genewise_disruption", sum(int(r["genewise_stops"]) + int(r["genewise_frameshifts"]) > 0 for r in loci),
         sum(q["disruption_loci"] > 0 for q in queries)),
        ("putative", len(filters["putative_pseudogene_ids"]), sum(q["putative_loci"] > 0 for q in queries)),
        ("family_copy_pass", len(filters["unitary_pseudogene_ids"]),
         sum(q["family_copy_pass_loci"] > 0 for q in queries)),
        ("final_candidate", len(filters["high_confidence_pseudogene_ids"]),
         sum(q["final_loci"] > 0 for q in queries)),
    ]
    write_tsv(qc / "stage_counts.tsv", ("stage", "loci_or_inputs", "queries"),
              ({"stage": stage, "loci_or_inputs": n, "queries": q}
               for stage, n, q in stage_counts))
    for source_name, target_name in (
        ("high_confidence_loci.bed", "final_loci.bed"),
        ("high_confidence_queries.txt", "final_queries.txt"),
        ("neighborhood_synteny.tsv", "neighborhood_synteny.tsv"),
    ):
        source = filter_dir / source_name
        if source.is_file():
            copy_file(source, results / target_name)
    for source_name in ("miniprot_intact_conflict_report.tsv",
                        "two_frameshift_only_review.txt"):
        source = filter_dir / source_name
        if source.is_file():
            copy_file(source, qc / source_name)
    copy_file(filter_dir / "effective_config.txt", metadata / "effective_config.txt")
    evidence = filter_dir / "evidence_status.tsv"
    if evidence.is_file():
        copy_file(evidence, metadata / "evidence_status.tsv")
    if budget.is_file():
        copy_file(budget, qc / "compute_budget_exclusions.tsv")
    failure = run_dir / "logs" / "genewise_failed.tsv"
    if not failure.is_file():
        failure = work / "logs" / "genewise_failed.tsv"
    if failure.is_file():
        copy_file(failure, qc / "genewise_failed.tsv")
    summary = {
        "input_queries": len(input_queries),
        "selected_loci": len(loci),
        "compute_budget_excluded_queries": len(budget_excluded & input_query_ids),
        "final_loci": len(filters["high_confidence_pseudogene_ids"]),
        "final_queries": sum(q["final_loci"] > 0 for q in queries),
        "genewise_stop_events": sum(int(r["genewise_stops"]) for r in loci),
        "genewise_frameshift_events": sum(int(r["genewise_frameshifts"]) for r in loci),
        "stop_events": sum(int(r["genewise_stops"]) for r in loci),
        "frameshift_events": sum(int(r["genewise_frameshifts"]) for r in loci),
        "workflow": "standard",
        "decision_reasons": dict(sorted(Counter(r["decision_reason"] for r in loci).items())),
        "evidence_levels": dict(sorted(Counter(r["evidence_level"] for r in loci).items())),
    }
    input_manifest_path = metadata / "input_manifest.json"
    if input_manifest_path.is_file():
        quality = json.loads(input_manifest_path.read_text()).get("protein_summary", {})
        summary["input_total_records"] = quality.get("input_records", len(input_queries))
        summary["input_excluded_records"] = quality.get("excluded_proteins", 0)
    write_json(results / "summary.json", summary)
    provenance = Path(provenance_file) if provenance_file else cache_file(run_dir)
    manifest = {
        "layout_version": 3,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "species": species,
        "comparison_mode": comparison_mode,
        "analysis_version": analysis_version or (ROOT / "VERSION").read_text().strip(),
        "publisher_version": (ROOT / "VERSION").read_text().strip(),
        "alignment_engine": "GeneWise",
        "workflow": "standard",
        "evidence_semantics": EVIDENCE_SEMANTICS,
        "coordinate_system": "0-based half-open; same boundaries as miniprot selected hit",
        "unitary_note": "The legacy unitary_pseudogene_ids.txt is a family/copy filter stage, not a verified unitary gene loss",
        "candidate_note": "Final cross-species loci are filtered GeneWise disruption candidates; synteny and relaxed selection require independent evidence",
        "inputs": {"proteins": str(Path(proteins).resolve())},
        "provenance": os.path.relpath(provenance.resolve(), metadata.resolve())
        if provenance.is_file() else None,
        "files": {
            "all_loci": "../results/loci.tsv",
            "all_queries": "../results/queries.tsv",
            "final_locus_table": "../results/final_loci.tsv",
            "final_loci": "../results/final_loci.bed",
            "final_queries": "../results/final_queries.txt",
            "stage_counts": "../qc/stage_counts.tsv",
            "nonfinal_loci": "../qc/nonfinal_loci.tsv",
            "effective_config": "effective_config.txt",
        },
        "counts": summary,
    }
    if budget.is_file():
        manifest["files"]["compute_budget_exclusions"] = "../qc/compute_budget_exclusions.tsv"
    if evidence.is_file():
        manifest["files"]["evidence_status"] = "evidence_status.tsv"
    for name in ("input_manifest.json", "protein_id_mapping.tsv", "reference_annotations.tsv"):
        if (metadata / name).is_file():
            manifest["files"][name.rsplit(".", 1)[0]] = name
    if neighborhood_path.is_file():
        manifest["files"]["neighborhood_synteny"] = "../results/neighborhood_synteny.tsv"
    if (qc / "miniprot_intact_conflict_report.tsv").is_file():
        manifest["files"]["miniprot_intact_conflict_report"] = "../qc/miniprot_intact_conflict_report.tsv"
    if (qc / "two_frameshift_only_review.txt").is_file():
        manifest["files"]["two_frameshift_only_review"] = "../qc/two_frameshift_only_review.txt"
    if failure.is_file():
        manifest["files"]["genewise_failures"] = "../qc/genewise_failed.tsv"
    for name in ("input_qc.json", "protein_exclusions.tsv"):
        if (qc / name).is_file():
            manifest["files"][name.split(".")[0]] = "../qc/" + name
    write_json(metadata / "run_manifest.json", manifest)
    publish_assets(run_dir)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", required=True, type=Path)
    parser.add_argument("--proteins", required=True, type=Path)
    parser.add_argument("--species", required=True)
    parser.add_argument("--comparison-mode", choices=("cross", "self"), default="cross")
    parser.add_argument("--analysis-version", help="Version used for an older run when publishing retrospectively")
    parser.add_argument("--budget-file", type=Path, help="Compute-budget exclusions from an older run")
    parser.add_argument("--provenance-file", type=Path, help="Provenance file from an older run")
    args = parser.parse_args()
    print(json.dumps(publish(args.run_dir, args.proteins, args.species,
                             args.analysis_version, args.budget_file,
                             args.provenance_file, args.comparison_mode), sort_keys=True))


if __name__ == "__main__":
    main()
