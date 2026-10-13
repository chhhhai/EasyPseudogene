#!/usr/bin/env python3
"""Search target-genome copies for GeneWise-disrupted reference queries.

Discovery still uses the original best hit. This PAF is evidence for the
existing functional-copy filter and does not introduce new candidate loci.
"""
import argparse
import csv
import re
import shlex
import subprocess
from pathlib import Path

from genewise_alignment import mutation_counts


def disrupted_loci(genewise):
    result = set()
    block = []
    with open(genewise) as handle:
        for line in handle:
            if line.strip() == "//":
                content = "".join(block)
                match = re.search(r"Query protein:\s+(\S+)", content)
                if match and sum(mutation_counts(content)) > 0:
                    result.add(match.group(1))
                block.clear()
            else:
                block.append(line)
    if block:
        content = "".join(block)
        match = re.search(r"Query protein:\s+(\S+)", content)
        if match and sum(mutation_counts(content)) > 0:
            result.add(match.group(1))
    return result


def locus_queries(locus_map, loci):
    result = set()
    found = set()
    with open(locus_map, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["locus_id"] in loci:
                found.add(row["locus_id"])
                result.add(row["query_id"])
    missing = loci - found
    if missing:
        raise ValueError(f"{len(missing)} disrupted loci are absent from locus_map.tsv")
    return result


def write_subset(source, output, queries):
    found = set()
    keep = False
    with open(source) as src, open(output, "w") as dst:
        for line in src:
            if line.startswith(">"):
                query = line[1:].split()[0]
                keep = query in queries
                if keep:
                    found.add(query)
            if keep:
                dst.write(line)
    missing = queries - found
    if missing:
        raise ValueError(f"{len(missing)} disrupted queries are absent from protein FASTA")
    return len(found)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--proteins", required=True)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--genewise", required=True)
    ap.add_argument("--locus-map", required=True)
    ap.add_argument("--output-paf", required=True)
    ap.add_argument("--threads", type=int, required=True)
    ap.add_argument("--outn", type=int, default=500)
    ap.add_argument("--outs", type=float, default=0.1)
    args = ap.parse_args()
    if args.threads < 1 or args.outn < 1 or not 0 < args.outs <= 1:
        ap.error("threads and outn must be positive; outs must be in (0, 1]")
    output = Path(args.output_paf)
    output.parent.mkdir(parents=True, exist_ok=True)
    queries = locus_queries(args.locus_map, disrupted_loci(args.genewise))
    subset = output.with_suffix(".queries.fa")
    count = write_subset(args.proteins, subset, queries)
    command = ["miniprot", "-t", str(args.threads), f"--outn={args.outn}",
               f"--outs={args.outs}", args.genome, str(subset)]
    output.with_suffix(".command.txt").write_text(shlex.join(command) + "\n")
    temporary = output.with_suffix(output.suffix + ".partial")
    if count:
        with open(temporary, "w") as stdout, open(output.with_suffix(".stderr.log"), "w") as stderr:
            subprocess.run(command, stdout=stdout, stderr=stderr, check=True)
    else:
        temporary.touch()
    temporary.replace(output)
    print(f"Candidate copy scan: {count} reference queries, {output.stat().st_size} PAF bytes")


if __name__ == "__main__":
    main()
