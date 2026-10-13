#!/usr/bin/env python3
"""Publish fast catalogs in the same results/qc/metadata layout as GeneWise.

Unvalidated screening hits remain visible in loci.tsv. Only structurally
validated disruptions with medium/high sequence evidence enter final_loci.tsv.
This rule does not reproduce the GeneWise family/copy filtering workflow.
"""
import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

from prefilter_proteins import read_fasta
from publish_genewise_results import LOCUS_FIELDS, QUERY_FIELDS, copy_file, write_json, write_tsv
from result_evidence import EVIDENCE_SEMANTICS, evidence
from publish_result_assets import publish_assets


def load(path):
    with Path(path).open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def number(row, key, kind=float):
    value = row.get(key, "")
    return kind(value) if value not in ("", None) else None


def publish(run_dir, species, mode, engine, proteins, catalog, validated=None, config=None, invocation=None):
    run_dir = Path(run_dir)
    raw = load(catalog)
    source = load(validated) if validated else raw
    raw_ids = [r["candidate_id"] for r in raw]
    ids = [r["candidate_id"] for r in source]
    if len(set(ids)) != len(ids) or len(set(raw_ids)) != len(raw_ids) or set(ids) != set(raw_ids):
        raise ValueError("Validated/raw catalogs must have identical unique candidate IDs")
    with open(proteins) as handle:
        query_ids = [h[1:].split()[0] for h, _ in read_fasta(handle)]
    if len(set(query_ids)) != len(query_ids):
        raise ValueError("Duplicate query IDs in normalized input")
    query_set = set(query_ids)
    loci, by_query = [], defaultdict(list)
    for row in source:
        query = row["query_id"]
        if query not in query_set:
            raise ValueError("Catalog query absent from input: " + query)
        start0, end0 = int(row["start"]) - 1, int(row["end"])
        if start0 < 0 or end0 <= start0:
            raise ValueError("Invalid catalog interval: " + row["candidate_id"])
        assessed = bool(validated and row.get("validator_status") == "validated_miniprot"
                        and row.get("alignment_source") == "miniprot" and row.get("structural_tags_complete") == "1")
        selected_parent = row.get("selected_parent", "")
        if assessed and selected_parent not in query_set:
            raise ValueError("Validated parent absent from protein input: " + selected_parent)
        validated_start = number(row, "validated_genomic_start", int)
        validated_end = number(row, "validated_genomic_end", int)
        interval_matches = bool(validated_start is not None and validated_end is not None and
                                validated_start > 0 and validated_end >= validated_start and
                                validated_start - 1 < end0 and validated_end > start0)
        assessed = assessed and interval_matches
        ev = evidence(number(row, "validated_stop_count", int) if assessed else None,
                      number(row, "validated_frameshift_count", int) if assessed else None,
                      number(row, "validated_query_coverage") if assessed else None,
                      number(row, "validated_identity") if assessed else None,
                      number(row, "validated_aligned_length", int) if assessed else None,
                      "miniprot" if assessed else row.get("alignment_source") or engine,
                      assessed, "PAF_match_fraction_proxy" if assessed else "")
        ev["evidence_query_id"] = selected_parent if assessed else ""
        if assessed:
            start0, end0 = validated_start - 1, validated_end
        final = assessed and ev["evidence_level"] in ("high", "medium") and (
            int(ev["stop_count"] or 0) + int(ev["frameshift_count"] or 0) > 0)
        locus = dict(locus_id=row["candidate_id"], query_id=query, chrom=row["chr"],
                     start0=start0, end0=end0, strand=row.get("validated_strand") if assessed else row["strand"], genewise_status="not_run",
                     genewise_stops="", genewise_frameshifts="",
                     stage="final_candidate" if final else "screening_candidate",
                     decision_reason="validated_disruption_support" if final else ev["evidence_reason"],
                     decision_detail="fast_sequence_evidence_selection; family/copy filters not applied",
                     neighborhood_synteny_status="not_assessed", neighborhood_synteny_reason="", **ev)
        loci.append(locus)
        by_query[query].append(locus)
    final = [r for r in loci if r["stage"] == "final_candidate"]
    queries = []
    for query in query_ids:
        items = by_query[query]
        disruptions = sum(bool(int(r["stop_count"] or 0) + int(r["frameshift_count"] or 0)) for r in items)
        final_count = sum(r["stage"] == "final_candidate" for r in items)
        queries.append(dict(query_id=query, selected_loci=len(items), disruption_loci=disruptions,
                            putative_loci=disruptions, family_copy_pass_loci="", final_loci=final_count,
                            status="final_candidate" if final_count else "screening_candidate" if items else "no_selected_locus"))
    results, qc, metadata = (run_dir / name for name in ("results", "qc", "metadata"))
    write_tsv(results / "loci.tsv", LOCUS_FIELDS, loci)
    write_tsv(results / "queries.tsv", QUERY_FIELDS, queries)
    write_tsv(results / "final_loci.tsv", LOCUS_FIELDS, final)
    write_tsv(qc / "nonfinal_loci.tsv", LOCUS_FIELDS, (r for r in loci if r["stage"] != "final_candidate"))
    results.mkdir(parents=True, exist_ok=True)
    with (results / "final_loci.bed").open("w") as handle:
        csv.writer(handle, delimiter="\t", lineterminator="\n").writerows(
            (r["chrom"], r["start0"], r["end0"], r["locus_id"], 0, r["strand"]) for r in final)
    (results / "final_queries.txt").write_text("".join(q + "\n" for q in sorted({r["query_id"] for r in final})))
    summary = dict(input_queries=len(query_ids), selected_loci=len(loci), final_loci=len(final),
                   final_queries=sum(q["final_loci"] > 0 for q in queries),
                   stop_events=sum(int(r["stop_count"] or 0) for r in loci),
                   frameshift_events=sum(int(r["frameshift_count"] or 0) for r in loci),
                   evidence_levels=dict(Counter(r["evidence_level"] for r in loci)),
                   workflow="fast", structural_validation=bool(validated))
    manifest_path = metadata / "input_manifest.json"
    if manifest_path.is_file():
        input_summary = json.loads(manifest_path.read_text()).get("protein_summary", {})
        summary.update(input_total_records=input_summary.get("input_records", len(query_ids)),
                       input_excluded_records=input_summary.get("excluded_proteins", 0))
    write_json(results / "summary.json", summary)
    stages = [("input", len(query_ids), len(query_ids)), ("selected_alignment", len(loci), len(by_query)),
              ("structural_validation", sum(r["evidence_level"] != "unassessed" for r in loci),
               sum(any(r["evidence_level"] != "unassessed" for r in by_query[q]) for q in query_ids)),
              ("disruption_detected", sum(bool(int(r["stop_count"] or 0) + int(r["frameshift_count"] or 0)) for r in loci),
               sum(q["disruption_loci"] > 0 for q in queries)), ("final_candidate", len(final), summary["final_queries"])]
    write_tsv(qc / "stage_counts.tsv", ("stage", "loci_or_inputs", "queries"),
              (dict(stage=s, loci_or_inputs=n, queries=q) for s, n, q in stages))
    write_tsv(metadata / "evidence_status.tsv", ("evidence", "status", "source"), [
        dict(evidence="sequence_disruption", status="assessed" if validated else "not_assessed", source=str(validated or catalog)),
        *[dict(evidence=name, status="not_assessed", source="") for name in
          ("family_copy_filter", "conserved_synteny", "relaxed_selection")]])
    if config:
        copy_file(Path(config), metadata / "config_overrides.conf")
    write_json(metadata / "run_manifest.json", dict(layout_version=3, species=species,
               comparison_mode=mode, workflow="fast", alignment_engine=engine,
               analysis_version=(Path(__file__).resolve().parents[1] / "VERSION").read_text().strip(),
               coordinate_system="0-based half-open; validated alignment intervals when assessed, otherwise original catalog intervals",
               evidence_semantics=EVIDENCE_SEMANTICS,
               final_selection_rule="validated miniprot disruption with medium/high sequence evidence; not GeneWise family/copy filtering",
               candidate_note="Computational candidates requiring biological review",
               files=dict(all_loci="../results/loci.tsv", all_queries="../results/queries.tsv",
                          final_locus_table="../results/final_loci.tsv", final_loci="../results/final_loci.bed",
                          final_queries="../results/final_queries.txt", stage_counts="../qc/stage_counts.tsv",
                          input_manifest="input_manifest.json", protein_mapping="protein_id_mapping.tsv",
                          effective_config="effective_config.txt", evidence_status="evidence_status.tsv",
                          input_qc="../qc/input_qc.json", protein_exclusions="../qc/protein_exclusions.tsv",
                          raw_catalog=str(Path(catalog).resolve()), validated_catalog=str(Path(validated).resolve()) if validated else None),
               expanded_configuration=dict(catalog=str(catalog), validated=str(validated) if validated else None,
                                           engine=engine, config_file=str(config) if config else None), counts=summary))
    manifest = json.loads((metadata / "run_manifest.json").read_text())
    manifest["arguments"] = invocation or []
    write_json(metadata / "run_manifest.json", manifest)
    publish_assets(run_dir)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("run-dir", "species", "engine", "proteins", "catalog"):
        parser.add_argument("--" + name, required=True)
    parser.add_argument("--mode", choices=("cross", "self"), required=True)
    parser.add_argument("--validated")
    parser.add_argument("--config")
    parser.add_argument("--invocation", default="[]", help="Original CLI arguments as JSON")
    args = parser.parse_args()
    print(json.dumps(publish(args.run_dir, args.species, args.mode, args.engine, args.proteins,
                             args.catalog, args.validated, args.config, json.loads(args.invocation)), sort_keys=True))


if __name__ == "__main__":
    main()
