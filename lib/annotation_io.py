"""Read GTF/GFF3 relationships without relying on a database ID namespace.

Exact IDs are retained. Only unambiguous aliases may join annotations to
protein headers; a dot in a BRAKER/MAKER identifier is not a version suffix.
"""
import csv
import gzip
import hashlib
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from urllib.parse import quote, unquote


def is_gzip(path):
    with open(path, "rb") as handle:
        return handle.read(2) == b"\x1f\x8b"


def open_text(path):
    return gzip.open(path, "rt", encoding="utf-8-sig") if is_gzip(path) else open(
        path, encoding="utf-8-sig")


def stable_id(identifier):
    # Keep g1.t1, maker IDs and arbitrary dotted local IDs intact.
    return re.sub(r"^(ENS[A-Z]*[GPT]\d+)\.\d+$", r"\1", identifier)


def aliases(identifier):
    values = {identifier, stable_id(identifier)}
    for prefix in ("gene:", "transcript:", "mRNA:", "protein:", "rna-", "cds-"):
        if identifier.startswith(prefix):
            value = identifier[len(prefix):]
            values.update((value, stable_id(value)))
    if identifier.startswith("gnl|"):
        parts = identifier.split("|")
        if len(parts) >= 3:
            values.add(parts[2])  # The middle token is a database name, not an ID.
    elif identifier.startswith(("ref|", "gb|", "emb|", "sp|", "tr|")):
        values.update(value for value in identifier.split("|")[1:] if value)
    return {value for value in values if value}


def parse_attributes(text):
    result = {}
    # GTF attributes are quoted; GFF3 values use URL escaping.
    for part in re.split(r';(?=(?:[^"]*"[^"]*")*[^"]*$)', text.strip()):
        part = part.strip()
        if not part or part == ".":
            continue
        match = re.match(r'^([^\s=]+)\s+"(.*)"$', part)
        if match:
            result[match[1]] = match[2]
        elif "=" in part:
            key, value = part.split("=", 1)
            result[key.strip()] = unquote(value.strip())
        elif " " in part:
            key, value = part.split(None, 1)
            result[key] = value.strip().strip('"')
    return result


def header_field(header, key):
    pattern = r'(?:^|[\s\[])' + re.escape(key) + r'[=:](?:"([^"]*)"|([^\s\]]+))'
    match = re.search(pattern, header)
    return unquote(match[1] if match[1] is not None else match[2]) if match else ""


def header_id(identifier):
    return quote(stable_id(identifier), safe=':@|+-_.')


@dataclass
class Feature:
    seqid: str
    source: str
    kind: str
    start: int
    end: int
    score: str
    strand: str
    phase: str
    attrs: dict
    parents: tuple


def features(path):
    with open_text(path) as handle:
        for number, line in enumerate(handle, 1):
            if line.startswith("##FASTA") or line.startswith(">"):
                break  # MAKER GFF3 may embed sequence after the annotation.
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"{path}:{number}: expected nine tab-separated annotation columns")
            try:
                start, end = int(fields[3]), int(fields[4])
            except ValueError as exc:
                raise ValueError(f"{path}:{number}: invalid annotation coordinates") from exc
            if start < 1 or end < start or fields[6] not in ("+", "-", ".", "?"):
                raise ValueError(f"{path}:{number}: invalid interval or strand")
            if fields[2] == "CDS" and fields[7] not in ("0", "1", "2", "."):
                raise ValueError(f"{path}:{number}: invalid CDS phase")
            attrs = parse_attributes(fields[8])
            # Split before URL decoding so an escaped comma remains part of an ID.
            parent = re.search(r'(?:^|;)\s*Parent=([^;]+)', fields[8])
            parents = tuple(unquote(x) for x in parent[1].split(",")) if parent else ()
            yield Feature(fields[0], fields[1], fields[2], start, end, fields[5],
                          fields[6], fields[7], attrs, parents)


class Annotation:
    def __init__(self, path):
        self.path = str(path)
        self.features = list(features(path))
        self.genes = {}
        self.transcripts = {}
        self.exons = defaultdict(set)
        self.cds = defaultdict(list)
        self.children = defaultdict(list)
        self.transcript_features = {}
        self.gene_features = {}
        self.proteins = defaultdict(set)
        self.tx_genes = defaultdict(set)
        self.gene_locations = {}
        self.gene_aliases = defaultdict(set)
        self.tx_aliases = defaultdict(set)
        self.protein_aliases = defaultdict(set)
        self.feature_transcripts = {}
        self.warnings = []
        self._inferred_gene_transcripts = set()
        if not self.features:
            raise ValueError(f"No annotation features in {path}")
        self._build()

    def _location(self, gene, feature):
        old = self.gene_locations.get(gene)
        value = (feature.seqid, feature.start - 1, feature.end, feature.strand)
        if old:
            if old[0] != value[0] or (old[3] in ("+", "-") and value[3] in ("+", "-") and old[3] != value[3]):
                raise ValueError(f"Gene {gene} occurs on incompatible contigs/strands in {self.path}")
            value = (old[0], min(old[1], value[1]), max(old[2], value[2]),
                     old[3] if old[3] in ("+", "-") else value[3])
        self.gene_locations[gene] = value

    def _gene(self, gene, feature):
        if not gene:
            return
        attrs = feature.attrs
        name = attrs.get("gene_name") or attrs.get("gene_symbol") or (
            attrs.get("Name", "") if feature.kind.lower() in ("gene", "pseudogene") else "")
        biotype = attrs.get("gene_biotype") or attrs.get("gene_type") or (
            attrs.get("biotype", "") if feature.kind.lower() in ("gene", "pseudogene") else "")
        if feature.kind.lower() == "pseudogene" and not biotype:
            biotype = "pseudogene"
        old = self.genes.get(gene, ("", ""))
        self.genes[gene] = (name or old[0], biotype or old[1])
        self._location(gene, feature)

    def _build(self):
        tx_features = []
        for feature in self.features:
            kind, attrs = feature.kind.lower(), feature.attrs
            if kind in ("gene", "pseudogene", "protein_coding_gene"):
                gene = attrs.get("ID") or attrs.get("gene_id")
                if gene:
                    self._gene(gene, feature)
                    self.gene_features[gene] = feature
                    for value in aliases(gene):
                        self.gene_aliases[value].add(gene)
                    if attrs.get("gene_id"):
                        for value in aliases(attrs["gene_id"]):
                            self.gene_aliases[value].add(gene)
            elif kind in ("mrna", "transcript", "pseudogenic_transcript") or kind.endswith("rna"):
                tx = attrs.get("ID") or attrs.get("transcript_id")
                if tx:
                    tx_features.append((tx, feature))
                    self.transcript_features[tx] = feature
                    if feature.parents:
                        self.tx_genes[tx].update(feature.parents)
                    elif attrs.get("gene_id"):
                        self.tx_genes[tx].add(attrs["gene_id"])

        # GTF often lacks gene/transcript rows; CDS/exon attributes supply them.
        for index, feature in enumerate(self.features):
            if feature.kind.lower() not in ("cds", "exon", "start_codon", "stop_codon", "utr", "five_prime_utr", "three_prime_utr"):
                continue
            attrs = feature.attrs
            txs = feature.parents or ((attrs["transcript_id"],) if attrs.get("transcript_id") else ())
            if feature.parents:
                inferred = []
                for parent in txs:
                    if parent in self.genes and parent not in self.transcript_features:
                        tx = parent + ".EPG_transcript"
                        if tx in self.transcript_features:
                            raise ValueError(f"Inferred transcript {tx} conflicts with an existing ID; provide transcript parents")
                        self.tx_genes[tx].add(parent)
                        inferred.append(tx)
                        message = f"Gene {parent}: child features directly parented to gene; inferred transcript {tx}"
                        if tx not in self._inferred_gene_transcripts:
                            self.warnings.append(message)
                            self._inferred_gene_transcripts.add(tx)
                    else:
                        inferred.append(parent)
                txs = tuple(inferred)
            if not txs:
                self.warnings.append(f"{feature.seqid}:{feature.start}-{feature.end}: {feature.kind} has no transcript parent")
            self.feature_transcripts[index] = txs
            for tx in txs:
                if attrs.get("gene_id") and (not feature.parents or not self.tx_genes[tx]):
                    self.tx_genes[tx].add(attrs["gene_id"])
                self.transcripts.setdefault(tx, ("", ""))
                self.children[tx].append(feature)
                if feature.kind.lower() == "exon":
                    self.exons[tx].add((feature.seqid, feature.start, feature.end, feature.strand))
                elif feature.kind.lower() == "cds":
                    self.cds[tx].append(feature)
                    protein = attrs.get("protein_id")
                    if protein:
                        self.proteins[protein].add(tx)
                    if attrs.get("ID"):
                        self.proteins[attrs["ID"]].add(tx)
                    if attrs.get("Name"):
                        self.proteins[attrs["Name"]].add(tx)

        for tx, feature in tx_features:
            biotype = feature.attrs.get("transcript_biotype") or feature.attrs.get("transcript_type") or feature.attrs.get("biotype", "")
            self.transcripts[tx] = ("", biotype)
            if feature.attrs.get("protein_id"):
                self.proteins[feature.attrs["protein_id"]].add(tx)

        for feature in self.features:
            if feature.kind.lower() not in ("polypeptide", "protein"):
                continue
            txs = set()
            parents = feature.parents or tuple(feature.attrs.get("Derives_from", "").split(","))
            for parent in parents:
                if parent in self.transcripts:
                    txs.add(parent)
                txs.update(self.proteins.get(parent, ()))
            for key in ("ID", "Name", "protein_id"):
                if feature.attrs.get(key):
                    self.proteins[feature.attrs[key]].update(txs)

        for tx in self.transcripts:
            parents = set()
            for parent in self.tx_genes[tx]:
                known = self._lookup(parent, self.gene_aliases)
                parents.update(known or {parent})
            if len(parents) > 1:
                raise ValueError(f"Transcript {tx} has multiple gene parents in {self.path}")
            gene = next(iter(parents), "")
            if not gene:
                # Keep an orphan transcript usable, but do not invent an isoform family.
                gene = "EPG_orphan_gene_" + hashlib.sha256(tx.encode()).hexdigest()[:20]
                if gene in self.genes:
                    raise ValueError(f"Inferred gene ID {gene} conflicts with an annotation gene")
                self.warnings.append(f"Transcript {tx}: no gene parent; inferred separate container {gene}")
            self.tx_genes[tx] = {gene}
            biotype = self.transcripts[tx][1]
            self.transcripts[tx] = (gene, biotype or ("protein_coding" if self.cds[tx] else ""))
            for feature in self.children[tx]:
                self._gene(gene, feature)
            for value in aliases(tx):
                self.tx_aliases[value].add(tx)

        for tx, feature in tx_features:
            self._gene(self.transcripts[tx][0], feature)
            for value in (feature.attrs.get("Name", ""), feature.attrs.get("transcript_id", "")):
                if value:
                    for alias in aliases(value):
                        self.tx_aliases[alias].add(tx)
        for gene in self.genes:
            for alias in aliases(gene):
                self.gene_aliases[alias].add(gene)
        for protein, txs in self.proteins.items():
            for alias in aliases(protein):
                self.protein_aliases[alias].update(txs)
        for tx, (gene, _) in self.transcripts.items():
            name, biotype = self.genes.get(gene, ("", ""))
            self.genes[gene] = (name, biotype or ("protein_coding" if self.cds[tx] else ""))
        if not self.transcripts:
            raise ValueError(f"No transcript/CDS relationships in {self.path}; provide GTF gene_id/transcript_id or GFF3 ID/Parent")

    def _lookup(self, identifier, table):
        found = set()
        if identifier in table:
            found.update(table[identifier])
        else:
            for alias in aliases(identifier):
                found.update(table.get(alias, ()))
        return found

    def resolve(self, header, explicit=None):
        query = header.split()[0]
        explicit = explicit or {}
        gene_hint = explicit.get("gene_id") or header_field(header, "gene_id") or header_field(header, "gene")
        tx_hint = explicit.get("transcript_id") or header_field(header, "transcript_id") or header_field(header, "transcript")
        query_txs = self._lookup(query, self.protein_aliases) | self._lookup(query, self.tx_aliases)
        txs = set()
        if tx_hint:
            txs = self._lookup(tx_hint, self.tx_aliases)
            if not txs and explicit:
                raise ValueError(f"ID map transcript {tx_hint} for {query} is absent from annotation")
            if txs and query_txs:
                overlap = txs & query_txs
                if not overlap:
                    raise ValueError(f"Protein/annotation transcript conflict for {query}: {tx_hint}")
                txs = overlap
        if not txs:
            txs = query_txs
        genes = self._lookup(gene_hint, self.gene_aliases) if gene_hint else set()
        if gene_hint and not genes and explicit:
            raise ValueError(f"ID map gene {gene_hint} for {query} is absent from annotation")
        if genes and txs:
            matching = {tx for tx in txs if self.transcripts[tx][0] in genes}
            if not matching:
                raise ValueError(f"Protein/annotation gene conflict for {query}: {gene_hint}")
            txs = matching
        if len(txs) == 1:
            tx = next(iter(txs))
            return self.transcripts[tx][0], tx, "explicit_map" if explicit else "annotation_relationship"
        if len(txs) > 1:
            possible_genes = {self.transcripts[tx][0] for tx in txs}
            if len(possible_genes) == 1:
                return next(iter(possible_genes)), "", "ambiguous_transcript"
            return "", "", "ambiguous"
        if len(genes) > 1:
            return "", "", "ambiguous"
        if len(genes) == 1:
            return next(iter(genes)), "", "gene_only"
        # A protein may intentionally use a gene ID rather than a transcript ID.
        genes = self._lookup(query, self.gene_aliases)
        if len(genes) == 1:
            return next(iter(genes)), "", "gene_only"
        return "", "", "unmatched"

    def write_gff3(self, path):
        def attributes(values):
            return ";".join(f"{key}={quote(str(value), safe=':@|+-_.')}" for key, value in values.items() if value)
        tx_ids = {}
        reserved = set(self.genes) | set(self.transcripts)
        for tx in self.transcripts:
            normalized = tx
            if tx in self.genes:
                normalized = "EPG_transcript_" + hashlib.sha256(tx.encode()).hexdigest()[:20]
                if normalized in reserved:
                    raise ValueError(f"Normalized transcript ID {normalized} conflicts with annotation IDs")
                reserved.add(normalized)
                self.warnings.append(f"Gene and transcript share ID {tx}; transcript container normalized to {normalized}")
            tx_ids[tx] = normalized
        with open(path, "w") as handle:
            handle.write("##gff-version 3\n")
            for gene, (chrom, start, end, strand) in sorted(self.gene_locations.items()):
                name, biotype = self.genes[gene]
                feature = self.gene_features.get(gene)
                values = dict(feature.attrs) if feature else {}
                values.pop("gene_id", None)
                values.update(ID=gene, Name=name, gene_biotype=biotype)
                handle.write("\t".join((chrom, "EasyPseudogene", "gene", str(start + 1), str(end), ".", strand, ".",
                                        attributes(values))) + "\n")
            for tx, (gene, biotype) in sorted(self.transcripts.items()):
                children = self.children[tx]
                explicit = self.transcript_features.get(tx)
                if explicit:
                    loc = (explicit.seqid, explicit.start, explicit.end, explicit.strand)
                    if children:
                        loc = (loc[0], min(loc[1], min(f.start for f in children)),
                               max(loc[2], max(f.end for f in children)), loc[3])
                elif children:
                    loc = (children[0].seqid, min(f.start for f in children), max(f.end for f in children), children[0].strand)
                else:
                    continue
                values = dict(explicit.attrs) if explicit else {}
                values.pop("gene_id", None)
                values.pop("transcript_id", None)
                values.update(ID=tx_ids[tx], Parent=gene, transcript_biotype=biotype)
                if tx_ids[tx] != tx:
                    values["transcript_id"] = tx
                kind = explicit.kind if explicit else "mRNA" if self.cds[tx] else "transcript"
                handle.write("\t".join((loc[0], "EasyPseudogene", kind, str(loc[1]), str(loc[2]), ".", loc[3], ".",
                                        attributes(values))) + "\n")
            for index, feature in enumerate(self.features):
                txs = self.feature_transcripts.get(index, ())
                if not txs:
                    continue
                values = {key: value for key, value in feature.attrs.items() if key not in ("Parent", "gene_id", "transcript_id")}
                encoded = attributes(values)
                parents = ",".join(quote(tx_ids[tx], safe=':@|+-_.') for tx in txs)
                encoded = (encoded + ";" if encoded else "") + "Parent=" + parents
                handle.write("\t".join((feature.seqid, feature.source, feature.kind, str(feature.start), str(feature.end),
                                        feature.score, feature.strand, feature.phase, encoded)) + "\n")


def read_id_map(path):
    if not path:
        return {}
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = set(reader.fieldnames or ())
        if "query_id" not in fields or not fields.intersection(("gene_id", "transcript_id")):
            raise ValueError("Protein ID map needs query_id and transcript_id and/or gene_id TSV columns")
        result = {}
        for row in reader:
            query = row["query_id"]
            if not query or query in result:
                raise ValueError(f"Empty or duplicate protein ID map query: {query}")
            result[query] = {key: (row.get(key) or "").strip() for key in ("gene_id", "transcript_id")}
            if not any(result[query].values()):
                raise ValueError(f"Protein ID map has no gene/transcript association for {query}")
        return result
