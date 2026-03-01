#!/usr/bin/env python3
import argparse


def parse_fasta_lengths(path):
    lengths = {}
    with open(path, "r") as fh:
        tid = None
        seq = []
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if tid is not None:
                    lengths[tid] = sum(len(s) for s in seq)
                tid = line[1:].split()[0]
                seq = []
            else:
                seq.append(line)
        if tid is not None:
            lengths[tid] = sum(len(s) for s in seq)
    return lengths


def parse_alignment(path):
    records = {}
    with open(path, "r") as fh:
        cur_id = None
        seqs = []
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if cur_id and len(seqs) >= 2:
                    records[cur_id] = (seqs[0], seqs[1])
                cur_id = line[1:].split()[0]
                seqs = []
            else:
                seqs.append(line)
        if cur_id and len(seqs) >= 2:
            records[cur_id] = (seqs[0], seqs[1])
    return records


def identity(a, b):
    if not a or not b:
        return 0.0
    L = min(len(a), len(b))
    if L == 0:
        return 0.0
    return sum(1 for i in range(L) if a[i] == b[i]) / L


def main():
    ap = argparse.ArgumentParser(description="Filter candidates by overall alignment quality.")
    ap.add_argument("--protein-alignment", required=True)
    ap.add_argument("--proteins", required=True)
    ap.add_argument("--in-list", required=True)
    ap.add_argument("--min-cov", type=float, default=0.7)
    ap.add_argument("--min-id", type=float, default=0.3)
    ap.add_argument("--min-aln-len", type=int, default=200)
    ap.add_argument("--out-pass", required=True)
    ap.add_argument("--out-fail", required=True)
    ap.add_argument("--out-report", required=True)
    args = ap.parse_args()

    ids = []
    with open(args.in_list, "r") as fh:
        for line in fh:
            tid = line.strip()
            if tid:
                ids.append(tid)

    lengths = parse_fasta_lengths(args.proteins)
    aln = parse_alignment(args.protein_alignment)

    with open(args.out_pass, "w") as outp, open(args.out_fail, "w") as outf, open(args.out_report, "w") as rep:
        rep.write("transcript_id\talign_len\tprotein_len\tcoverage\tidentity\tstatus\n")
        for tid in ids:
            q, t = aln.get(tid, ("", ""))
            aln_len = min(len(q), len(t))
            prot_len = lengths.get(tid, 0)
            cov = aln_len / max(1, prot_len) if prot_len else 0.0
            ident = identity(q, t) if aln_len else 0.0
            if aln_len >= args.min_aln_len and cov >= args.min_cov and ident >= args.min_id:
                outp.write(f"{tid}\n")
                rep.write(f"{tid}\t{aln_len}\t{prot_len}\t{cov:.4f}\t{ident:.4f}\tpass\n")
            else:
                outf.write(f"{tid}\n")
                rep.write(f"{tid}\t{aln_len}\t{prot_len}\t{cov:.4f}\t{ident:.4f}\tfail\n")


if __name__ == "__main__":
    main()
