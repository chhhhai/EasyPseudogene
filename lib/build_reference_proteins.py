#!/usr/bin/env python3
"""Translate an annotated reference genome and retain one protein per gene.

gffread handles CDS splicing, strand and phase. A matching GTF/GFF3 is required;
two unannotated genome FASTAs alone cannot define reference coding genes.
"""
import argparse
import gzip
import shutil
import subprocess
from collections import defaultdict
from pathlib import Path

from annotation_io import Annotation, header_id, open_text


def fasta_records(path):
    name, chunks = None, []
    with open_text(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def build(genome, gtf, output, min_length):
    model = Annotation(gtf)
    genes, transcripts = model.genes, model.transcripts
    tx_protein = {}
    tx_protein_candidates = defaultdict(set)
    for feature in model.features:
        protein_id = feature.attrs.get("protein_id")
        if protein_id and len(model.proteins[protein_id]) == 1:
            tx_protein_candidates[next(iter(model.proteins[protein_id]))].add(protein_id)
    for tx_id, candidates in tx_protein_candidates.items():
        if len(candidates) > 1:
            raise ValueError(f"Transcript {tx_id} has multiple protein IDs; supply its explicit protein FASTA rather than choosing an ORF implicitly")
        tx_protein[tx_id] = next(iter(candidates))
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    raw = output.with_suffix(output.suffix + ".gffread.partial")
    gtf_source = Path(gtf)
    unpacked = None
    if str(gtf_source).endswith(".gz"):
        unpacked = output.with_suffix(output.suffix + ".input.gtf.partial")
        with gzip.open(gtf_source, "rb") as source, open(unpacked, "wb") as target:
            shutil.copyfileobj(source, target)
        gtf_source = unpacked
    command = ["gffread", "-g", str(genome), "-S", "-y", str(raw), str(gtf_source)]
    try:
        subprocess.run(command, stdout=subprocess.DEVNULL, check=True)
    finally:
        if unpacked is not None:
            unpacked.unlink(missing_ok=True)
    selected = {}
    translated = 0
    for tx_id, sequence in fasta_records(raw):
        translated += 1
        gene_id = transcripts.get(tx_id, ("", ""))[0]
        if not gene_id:
            continue
        gene_biotype = genes.get(gene_id, ("", ""))[1]
        if gene_biotype and gene_biotype != "protein_coding":
            continue
        sequence = sequence.rstrip("*.").upper()
        if len(sequence) < min_length or any(base not in "ACDEFGHIKLMNPQRSTVWY" for base in sequence):
            continue
        feature = model.transcript_features.get(tx_id)
        original_tx = feature.attrs.get("transcript_id", tx_id) if feature else tx_id
        protein_id = tx_protein.get(tx_id, original_tx)
        old = selected.get(gene_id)
        if old is None or (len(sequence), protein_id) > (len(old[2]), old[0]):
            selected[gene_id] = (protein_id, tx_id, sequence)
    if not selected:
        raise ValueError("No valid translated proteins; check reference genome/GTF assembly compatibility")
    ids = [protein_id for protein_id, _, _ in selected.values()]
    if len(ids) != len(set(ids)):
        raise ValueError("GTF maps the same protein identifier to multiple selected genes")
    temporary = output.with_suffix(output.suffix + ".partial")
    with open(temporary, "w") as handle:
        for gene_id, (protein_id, tx_id, sequence) in sorted(selected.items()):
            handle.write(f">{protein_id} gene:{header_id(gene_id)} transcript:{header_id(tx_id)}\n")
            for i in range(0, len(sequence), 70):
                handle.write(sequence[i:i + 70] + "\n")
    temporary.replace(output)
    raw.unlink()
    print(f"Reference translation: {translated} GTF transcripts, {len(selected)} representative proteins")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--annotation", "--gtf", dest="gtf", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--min-length", type=int, default=30)
    args = ap.parse_args()
    build(args.genome, args.gtf, args.output, args.min_length)


if __name__ == "__main__":
    main()
