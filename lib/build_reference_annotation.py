#!/usr/bin/env python3
"""Join reference protein IDs to GTF/GFF3 gene/transcript relationships."""
import argparse
import csv

from annotation_io import Annotation, header_field, open_text, parse_attributes, read_id_map, stable_id


FIELDS = ("query_id", "gene_id", "gene_symbol", "gene_biotype",
          "transcript_id", "transcript_biotype", "exon_count",
          "mapping_status", "annotation_gene_id", "annotation_transcript_id", "exon_count_source")


def read_gtf(path):
    """Compatibility interface; accepts both GTF and GFF3."""
    model = Annotation(path)
    genes = {stable_id(gene): value for gene, value in model.genes.items()}
    transcripts = {stable_id(tx): (stable_id(gene), biotype)
                   for tx, (gene, biotype) in model.transcripts.items()}
    exons = {stable_id(tx): len(values) for tx, values in model.exons.items()}
    proteins = {}
    for protein, txs in model.proteins.items():
        if len(txs) == 1:
            tx = next(iter(txs))
            proteins[stable_id(protein)] = (stable_id(model.transcripts[tx][0]), stable_id(tx))
    return genes, transcripts, exons, proteins


def build(fasta, gtf, output, id_map=None):
    model = Annotation(gtf)
    explicit = read_id_map(id_map)
    total = matched = 0
    seen = set()
    with open_text(fasta) as source, open(output, "w", newline="") as target:
        writer = csv.DictWriter(target, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for line in source:
            if not line.startswith(">"):
                continue
            header = line[1:].strip()
            if not header:
                raise ValueError("Empty protein FASTA header")
            query_id = header.split()[0]
            if query_id in seen:
                raise ValueError("Duplicate protein FASTA ID: " + query_id)
            seen.add(query_id)
            gene, tx, status = model.resolve(header, explicit.get(query_id))
            symbol, gene_biotype = model.genes.get(gene, ("", ""))
            _, tx_biotype = model.transcripts.get(tx, ("", ""))
            exon_count = len(model.exons.get(tx, ())) if tx else ""
            exon_source = "annotated_exons" if exon_count else ""
            if tx and not exon_count and model.cds[tx]:
                exon_count = len({(f.seqid, f.start, f.end, f.strand) for f in model.cds[tx]})
                exon_source = "CDS_blocks"
            writer.writerow(dict(query_id=query_id, gene_id=stable_id(gene),
                                 gene_symbol=symbol or header_field(header, "gene_symbol"),
                                 gene_biotype=gene_biotype, transcript_id=stable_id(tx),
                                 transcript_biotype=tx_biotype, exon_count=exon_count,
                                 mapping_status=status, annotation_gene_id=gene,
                                 annotation_transcript_id=tx, exon_count_source=exon_source))
            total += 1
            matched += bool(gene)
    return total, matched


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--proteins", required=True)
    parser.add_argument("--annotation", "--gtf", dest="annotation", required=True)
    parser.add_argument("--protein-id-map")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    total, matched = build(args.proteins, args.annotation, args.output, args.protein_id_map)
    print(f"Reference annotation: {matched}/{total} proteins linked to genes")


if __name__ == "__main__":
    main()
