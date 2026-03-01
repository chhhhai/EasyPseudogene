#!/usr/bin/env python3
import argparse
import urllib.parse
import urllib.request
import sys


DEFAULT_URL = "https://www.ensembl.org/biomart/martservice"


def build_query(dataset):
    return f"""<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query virtualSchemaName="default" formatter="TSV" header="0" uniqueRows="1" count="" datasetConfigVersion="0.6">
  <Dataset name="{dataset}" interface="default">
    <Attribute name="ensembl_gene_id"/>
    <Attribute name="ensembl_transcript_id"/>
    <Attribute name="gene_biotype"/>
    <Attribute name="transcript_biotype"/>
  </Dataset>
</Query>
"""


def main():
    ap = argparse.ArgumentParser(description="Download Ensembl gene/transcript biotypes via BioMart.")
    ap.add_argument("--out", required=True, help="Output TSV path")
    ap.add_argument("--mart-url", default=DEFAULT_URL, help="BioMart martservice URL")
    ap.add_argument("--dataset", default="hsapiens_gene_ensembl", help="BioMart dataset name")
    args = ap.parse_args()

    query = build_query(args.dataset)
    data = urllib.parse.urlencode({"query": query}).encode("utf-8")

    req = urllib.request.Request(args.mart_url, data=data, method="POST")
    try:
        with urllib.request.urlopen(req) as resp, open(args.out, "wb") as out:
            out.write(resp.read())
    except Exception as e:
        sys.stderr.write(f"ERROR: failed to download BioMart data: {e}\n")
        sys.exit(1)


if __name__ == "__main__":
    main()
