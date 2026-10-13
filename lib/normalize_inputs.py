"""Stage generic annotation/FASTA inputs while retaining original query IDs."""
import argparse
import csv
import hashlib
import json
import re
import shutil
from collections import Counter
from pathlib import Path

from annotation_io import Annotation, header_field, header_id, is_gzip, open_text, read_id_map, stable_id
from run_layout import cache_file


def sha(path):
    value = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(chunk)
    return value.hexdigest()


def install(temporary, destination):
    """Keep existing indexes valid when staging byte-identical inputs again."""
    if destination.exists() and sha(temporary) == sha(destination):
        temporary.unlink()
        return
    temporary.replace(destination)
    for suffix in (".fai", ".gzi"):
        Path(str(destination) + suffix).unlink(missing_ok=True)


def plain_genome(source, destination):
    if not is_gzip(source):
        if Path(source).resolve() == destination.resolve():
            return destination.resolve()
        temporary = destination.with_suffix(destination.suffix + ".partial")
        with open(source, "rb") as handle:
            has_bom = handle.read(3) == b"\xef\xbb\xbf"
        if has_bom:
            with open_text(source) as handle, open(temporary, "w") as output:
                shutil.copyfileobj(handle, output)
        else:
            shutil.copyfile(source, temporary)
        install(temporary, destination)
        return destination.resolve()
    temporary = destination.with_suffix(destination.suffix + ".partial")
    with open_text(source) as handle, open(temporary, "w") as output:
        shutil.copyfileobj(handle, output)
    install(temporary, destination)
    return destination.resolve()


def genome_lengths(path):
    lengths = {}
    identifier = None
    with open_text(path) as handle:
        for number, line in enumerate(handle, 1):
            if line.startswith(">"):
                fields = line[1:].split()
                if not fields or fields[0] in lengths:
                    raise ValueError(f"{path}:{number}: empty or duplicate genome sequence ID")
                identifier = fields[0]
                lengths[identifier] = 0
            elif line.strip():
                if identifier is None:
                    raise ValueError(f"{path}:{number}: sequence before a FASTA header")
                sequence = "".join(line.split())
                if not re.fullmatch(r"[ACGTRYSWKMBDHVNacgtryswkmbdhvn]+", sequence):
                    raise ValueError(f"{path}:{number}: invalid DNA symbols in genome FASTA")
                lengths[identifier] += len(sequence)
    if not lengths or any(length == 0 for length in lengths.values()):
        raise ValueError(f"Empty genome FASTA record in {path}")
    return lengths


def validate_annotation(model, lengths):
    for feature in model.features:
        if feature.seqid not in lengths:
            raise ValueError(f"Annotation contig {feature.seqid!r} is absent from its genome; use files from the same assembly (names are not guessed)")
        if feature.end > lengths[feature.seqid]:
            raise ValueError(f"Annotation interval {feature.seqid}:{feature.start}-{feature.end} exceeds its genome sequence")


def fasta_records(path):
    header, chunks = None, []
    with open_text(path) as handle:
        for number, line in enumerate(handle, 1):
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header, chunks = line[1:].strip(), []
                if not header:
                    raise ValueError(f"{path}:{number}: empty protein FASTA header")
            elif line.strip():
                if header is None:
                    raise ValueError(f"{path}:{number}: sequence before a FASTA header")
                chunks.append("".join(line.split()))
    if header is not None:
        yield header, "".join(chunks)


def normalize_proteins(source, destination, model, id_map, report, invalid_policy="error", exclusions=None):
    seen, statuses = set(), Counter()
    excluded = []
    retained = 0
    temporary = destination.with_suffix(destination.suffix + ".partial")
    fields = ("query_id", "original_header", "annotation_gene_id", "annotation_transcript_id",
              "normalized_gene_id", "normalized_transcript_id", "mapping_status", "terminal_stop_removed")
    try:
        with open(temporary, "w") as output, open(report, "w", newline="") as table:
            writer = csv.DictWriter(table, fieldnames=fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for header, raw_sequence in fasta_records(source):
                query = header.split()[0]
                if query in seen:
                    raise ValueError(f"Duplicate protein FASTA ID: {query}; retain unique transcript/protein IDs")
                seen.add(query)
                sequence = raw_sequence.upper().rstrip("*")
                invalid = set(sequence) - set("ACDEFGHIKLMNPQRSTVWYBXZJUO")
                if not sequence or invalid:
                    reason = "empty_sequence" if not sequence else "internal_stop" if "*" in invalid else "invalid_amino_acid_symbols"
                    excluded.append(dict(query_id=query, reason=reason, invalid_symbols=",".join(sorted(invalid)),
                                         internal_stops=sequence.count("*"), original_header=header))
                    if invalid_policy == "error":
                        raise ValueError(f"Protein {query}: {reason}; supply unaligned protein FASTA or use --invalid-proteins skip (exclusions are reported)")
                    continue
                retained += 1
                if model:
                    gene, tx, status = model.resolve(header, id_map.get(query))
                elif query in id_map:
                    gene, tx, status = id_map[query].get("gene_id", ""), id_map[query].get("transcript_id", ""), "explicit_map"
                else:
                    gene, tx, status = header_field(header, "gene") or header_field(header, "gene_id"), header_field(header, "transcript") or header_field(header, "transcript_id"), "header_only"
                statuses[status] += 1
                normalized = header
                if gene:
                    # Sidecar retains the original header; first-token query ID never changes.
                    parts = header.split(maxsplit=1)
                    description = parts[1] if len(parts) == 2 else ""
                    description = re.sub(r'(?:^|\s)(?:gene_id|gene|transcript_id|transcript)[=:](?:"[^"]*"|[^\s]+)', "", description).strip()
                    normalized = query + (" " + description if description else "")
                    normalized += " gene:" + header_id(gene)
                    if tx:
                        normalized += " transcript:" + header_id(tx)
                output.write(">" + normalized + "\n")
                for start in range(0, len(sequence), 70):
                    output.write(sequence[start:start + 70] + "\n")
                writer.writerow(dict(query_id=query, original_header=header, annotation_gene_id=gene,
                                     annotation_transcript_id=tx, normalized_gene_id=stable_id(gene),
                                     normalized_transcript_id=stable_id(tx), mapping_status=status,
                                     terminal_stop_removed=int(len(sequence) != len(raw_sequence))))
        if not retained:
            raise ValueError(f"No valid protein FASTA records in {source}; see the protein exclusion report")
        install(temporary, destination)
    except Exception:
        temporary.unlink(missing_ok=True)
        raise
    finally:
        if exclusions:
            temporary_report = exclusions.with_suffix(".tsv.partial")
            with open(temporary_report, "w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=("query_id", "reason", "invalid_symbols", "internal_stops", "original_header"),
                                        delimiter="\t", lineterminator="\n")
                writer.writeheader()
                writer.writerows(excluded)
            install(temporary_report, exclusions)
    return dict(input_records=len(seen), proteins=retained, excluded_proteins=len(excluded),
                invalid_protein_policy=invalid_policy, mapping_statuses=dict(statuses),
                unused_id_map_entries=len(set(id_map) - seen))


def prepare(args):
    directory = args.directory.resolve()
    work = Path(getattr(args, 'work_dir', None) or directory).resolve()
    inputs = work / "inputs"
    if work != directory and (directory / 'inputs').exists():
        raise ValueError('Legacy output layout exists; use a new --output directory')
    sources = {}
    for key in ("proteins", "genome", "reference_genome", "reference_annotation", "target_annotation", "protein_id_map"):
        path = getattr(args, key)
        if path:
            sources[key] = dict(path=str(Path(path).resolve()), sha256=sha(path))
    library = Path(__file__).parent
    signature = dict(sources=sources, mode=args.mode, invalid_protein_policy=args.invalid_proteins,
                     analysis_version=(library.parent / "VERSION").read_text().strip(),
                     normalizer_source={name: sha(library / name) for name in (
                         "normalize_inputs.py", "annotation_io.py", "build_reference_proteins.py", "run_layout.py")})
    manifest = inputs / "input_manifest.json"
    if cache_file(directory).exists() or (directory / "reference_cache.json").exists() or (work / "fast_catalog").exists() or (directory / "metadata/run_manifest.json").exists():
        old = json.loads(manifest.read_text()) if manifest.exists() else {}
        if old.get("signature") != signature:
            raise ValueError("Input normalization differs from this existing run; use a new --output directory")
    inputs.mkdir(parents=True, exist_ok=True)
    genome = plain_genome(args.genome, inputs / "target.genome.fa")
    target_lengths = genome_lengths(genome)
    annotation_path = args.reference_annotation if args.mode == "cross" else args.target_annotation
    model = Annotation(annotation_path) if annotation_path else None
    normalized_annotation = ""
    if model:
        normalized_annotation = inputs / ("reference.annotation.gff3" if args.mode == "cross" else "target.annotation.gff3")
        if args.mode == "self":
            validate_annotation(model, target_lengths)
        temporary = normalized_annotation.with_suffix(".gff3.partial")
        model.write_gff3(temporary)
        install(temporary, normalized_annotation)
    reference_genome = ""
    protein_source = args.proteins
    if args.reference_genome:
        from build_reference_proteins import build
        reference_genome = plain_genome(args.reference_genome, inputs / "reference.genome.fa")
        validate_annotation(model, genome_lengths(reference_genome))
        protein_source = inputs / "translated.reference.proteins.fa"
        build(reference_genome, normalized_annotation, protein_source, 30)
    proteins = inputs / "reference.proteins.fa"
    report = inputs / "protein_id_mapping.tsv"
    protein_summary = normalize_proteins(protein_source, proteins, model,
                                         read_id_map(args.protein_id_map), report, args.invalid_proteins,
                                         inputs / "protein_exclusions.tsv")
    paths = dict(proteins=str(proteins), genome=str(genome),
                 reference_genome=str(reference_genome),
                 reference_annotation=str(normalized_annotation) if args.mode == "cross" else "",
                 target_annotation=str(normalized_annotation) if args.mode == "self" else "")
    payload = dict(signature=signature, normalized_paths=paths, protein_summary=protein_summary,
                   target_contigs=len(target_lengths), annotation_warnings=model.warnings if model else [],
                   original_query_ids_preserved=True,
                   mapping_note="Unmatched/ambiguous annotation joins remain in the protein input; use --protein-id-map for explicit associations")
    temporary = manifest.with_suffix(".json.partial")
    temporary.write_text(json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n")
    install(temporary, manifest)
    quality = dict(status="passed_with_warnings" if protein_summary["excluded_proteins"] or model and model.warnings or
                   any(protein_summary["mapping_statuses"].get(s) for s in ("unmatched", "ambiguous", "ambiguous_transcript")) else "passed",
                   protein_summary=protein_summary, target_contigs=len(target_lengths),
                   target_bases=sum(target_lengths.values()), annotation_warnings=model.warnings if model else [],
                   protein_exclusions="protein_exclusions.tsv", protein_mapping="protein_id_mapping.tsv")
    temporary = inputs / "input_qc.json.partial"
    temporary.write_text(json.dumps(quality, indent=2, sort_keys=True) + "\n")
    install(temporary, inputs / "input_qc.json")
    print("Input preparation: " + json.dumps(protein_summary, ensure_ascii=False))
    if protein_summary["excluded_proteins"]:
        print(f"Excluded {protein_summary['excluded_proteins']} invalid reference proteins; see inputs/protein_exclusions.tsv")
    if model and any(protein_summary["mapping_statuses"].get(status) for status in ("ambiguous", "ambiguous_transcript", "unmatched")):
        print("Some proteins lack a unique annotation join; see inputs/protein_id_mapping.tsv")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=("cross", "self"), required=True)
    parser.add_argument("--directory", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, help="Directory for staged inputs; defaults to --directory")
    parser.add_argument("--genome", type=Path, required=True)
    parser.add_argument("--invalid-proteins", choices=("error", "skip"), default="error")
    for name in ("proteins", "reference-genome", "reference-annotation", "target-annotation", "protein-id-map"):
        parser.add_argument("--" + name, type=Path)
    args = parser.parse_args()
    if bool(args.proteins) == bool(args.reference_genome):
        parser.error("provide proteins or an annotated reference genome")
    if args.reference_genome and not args.reference_annotation:
        parser.error("reference genome translation requires a reference annotation")
    if args.mode == "self" and not args.target_annotation:
        parser.error("self mode requires a target annotation")
    if args.mode == "self" and (args.reference_genome or args.reference_annotation):
        parser.error("self mode uses proteins, genome and target annotation")
    try:
        prepare(args)
    except (OSError, ValueError) as exc:
        try:
            args.directory.mkdir(parents=True, exist_ok=True)
            failure = (args.work_dir or args.directory) / "input_error.json"
            failure.parent.mkdir(parents=True, exist_ok=True)
            temporary = failure.with_suffix(".json.partial")
            temporary.write_text(json.dumps(dict(status="failed", error=str(exc)), indent=2) + "\n")
            temporary.replace(failure)
        except OSError:
            pass
        raise SystemExit("ERROR: " + str(exc)) from exc


if __name__ == "__main__":
    main()
