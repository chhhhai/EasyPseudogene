#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Ensembl human: pick the longest transcript per gene (by CDS length), then output the corresponding protein.

You currently have:
  - Homo_sapiens.GRCh38.cds.all.fa   (FASTA records keyed by ENST...)
  - Homo_sapiens.GRCh38.pep.all.fa   (FASTA records keyed by ENSP..., header contains transcript:ENST...)

This script:
  1) Reads pep FASTA once to find, for each ENST, the "best" protein record (max AA length; tie-break by ENSP).
  2) Reads CDS FASTA to compute CDS length per ENST and pick, for each ENSG gene, the ENST with longest CDS
     that also has a protein in step (1). (tie-break: longer protein; then ENST)
  3) Reads pep FASTA again to extract only the chosen proteins and write a reduced FASTA + TSV report.

Important:
  - With CDS FASTA we select "longest CDS", not full transcript with UTR.
    If you want "longest transcript by cDNA length", download Ensembl cdna.all.fa and we can switch to that.
"""

from __future__ import annotations

import argparse
import re
from dataclasses import dataclass
from typing import Dict, Iterable, Iterator, List, Optional, Tuple


RE_ENSG = re.compile(r"\bgene:(ENSG[0-9]+(?:\.[0-9]+)?)\b")
RE_ENST = re.compile(r"\btranscript:(ENST[0-9]+(?:\.[0-9]+)?)\b")
RE_GENE_SYMBOL = re.compile(r"\bgene_symbol:([^\s]+)\b")


def iter_fasta(path: str) -> Iterator[Tuple[str, str, str]]:
    """Yield (first_token_id, full_header, seq). header excludes leading '>'."""
    seq_id: Optional[str] = None
    header: Optional[str] = None
    chunks: List[str] = []

    def flush():
        nonlocal seq_id, header, chunks
        if seq_id is not None and header is not None:
            yield (seq_id, header, "".join(chunks))
        seq_id, header, chunks = None, None, []

    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                yield from flush()
                header = line[1:]
                seq_id = header.split()[0]
            else:
                if seq_id is None:
                    raise ValueError(f"FASTA parse error in {path}: sequence line before header")
                chunks.append(line.strip())
        yield from flush()


def wrap_seq(seq: str, width: int = 60) -> str:
    return "\n".join(seq[i : i + width] for i in range(0, len(seq), width))


@dataclass(frozen=True)
class PepBest:
    ensp: str
    pep_len: int


@dataclass(frozen=True)
class GenePick:
    ensg: str
    gene_symbol: str
    enst: str
    cds_len: int
    ensp: str
    pep_len: int


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cds", required=True, help="Ensembl CDS FASTA (ENST...)")
    ap.add_argument("--pep", required=True, help="Ensembl peptide FASTA (ENSP..., header has transcript:ENST...)")
    ap.add_argument("--out-prefix", default="ensembl_longest", help="Output prefix")
    args = ap.parse_args()

    # 1) pep: ENST -> best protein (length, ENSP)
    pep_best_by_enst: Dict[str, PepBest] = {}
    pep_total = 0
    pep_with_enst = 0
    for ensp, header, seq in iter_fasta(args.pep):
        pep_total += 1
        m_tx = RE_ENST.search(header)
        if not m_tx:
            continue
        pep_with_enst += 1
        enst = m_tx.group(1)
        pep_len = len(seq)
        prev = pep_best_by_enst.get(enst)
        cand = PepBest(ensp=ensp, pep_len=pep_len)
        if prev is None or (cand.pep_len > prev.pep_len) or (cand.pep_len == prev.pep_len and cand.ensp < prev.ensp):
            pep_best_by_enst[enst] = cand

    # 2) cds: pick per gene (only if protein exists for that ENST)
    best_by_gene: Dict[str, GenePick] = {}
    cds_total = 0
    cds_with_gene = 0
    cds_with_protein = 0

    for enst, header, seq in iter_fasta(args.cds):
        cds_total += 1
        m_g = RE_ENSG.search(header)
        if not m_g:
            continue
        cds_with_gene += 1
        ensg = m_g.group(1)
        m_sym = RE_GENE_SYMBOL.search(header)
        gene_symbol = m_sym.group(1) if m_sym else ""
        cds_len = len(seq)

        best_pep = pep_best_by_enst.get(enst)
        if best_pep is None:
            continue
        cds_with_protein += 1

        cand = GenePick(
            ensg=ensg,
            gene_symbol=gene_symbol,
            enst=enst,
            cds_len=cds_len,
            ensp=best_pep.ensp,
            pep_len=best_pep.pep_len,
        )

        prev = best_by_gene.get(ensg)
        if prev is None:
            best_by_gene[ensg] = cand
        else:
            # choose longest CDS; tie-break by protein length; then ENST
            if (
                (cand.cds_len > prev.cds_len)
                or (cand.cds_len == prev.cds_len and cand.pep_len > prev.pep_len)
                or (cand.cds_len == prev.cds_len and cand.pep_len == prev.pep_len and cand.enst < prev.enst)
            ):
                best_by_gene[ensg] = cand

    # build quick lookup: (ENST -> chosen ENSP) for output extraction
    chosen_enst_to_ensp = {p.enst: p.ensp for p in best_by_gene.values()}

    # 3) pep second pass: extract only chosen proteins
    out_fa = f"{args.out_prefix}.longest_cds.proteins.fa"
    out_tsv = f"{args.out_prefix}.longest_cds.tsv"

    # write TSV
    with open(out_tsv, "w", encoding="utf-8") as w:
        w.write("gene\tenst\tensp\tgene_symbol\tcds_len\tpep_len\n")
        for ensg in sorted(best_by_gene.keys()):
            p = best_by_gene[ensg]
            w.write(f"{p.ensg}\t{p.enst}\t{p.ensp}\t{p.gene_symbol}\t{p.cds_len}\t{p.pep_len}\n")

    # write FASTA (ID uses ENST, since your downstream is "longest transcript -> protein")
    enst_to_pick: Dict[str, GenePick] = {p.enst: p for p in best_by_gene.values()}
    with open(out_fa, "w", encoding="utf-8") as w:
        wrote = 0
        for ensp, header, seq in iter_fasta(args.pep):
            m_tx = RE_ENST.search(header)
            if not m_tx:
                continue
            enst = m_tx.group(1)
            pick = enst_to_pick.get(enst)
            if pick is None:
                continue
            if pick.ensp != ensp:
                continue
            desc = f"gene={pick.ensg} ensp={pick.ensp} gene_symbol={pick.gene_symbol} cds_len={pick.cds_len} pep_len={pick.pep_len}"
            w.write(f">{pick.enst} {desc}\n")
            w.write(wrap_seq(seq) + "\n")
            wrote += 1

    print(f"pep_total_records: {pep_total}")
    print(f"pep_with_transcript_tag: {pep_with_enst}")
    print(f"cds_total_records: {cds_total}")
    print(f"cds_with_gene_tag: {cds_with_gene}")
    print(f"cds_with_protein: {cds_with_protein}")
    print(f"genes_picked: {len(best_by_gene)}")
    print(f"proteins_written: {wrote}")
    print(f"wrote: {out_tsv}")
    print(f"wrote: {out_fa}")


if __name__ == "__main__":
    main()

