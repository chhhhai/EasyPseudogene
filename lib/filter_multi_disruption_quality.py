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


def parse_mutation_counts(path):
    counts = {}
    with open(path, "r") as fh:
        for line in fh:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            tid = parts[0]
            stop = int(parts[1])
            shift = int(parts[2])
            counts[tid] = (stop, shift, stop + shift)
    return counts


def calc_identity(qseq, tseq):
    if not qseq or not tseq:
        return 0.0
    matches = 0
    length = min(len(qseq), len(tseq))
    for i in range(length):
        if qseq[i] == tseq[i]:
            matches += 1
    return matches / max(1, length)


def main():
    ap = argparse.ArgumentParser(description="Filter multi-disruption candidates by alignment quality.")
    ap.add_argument("--protein-alignment", required=True)
    ap.add_argument("--proteins", required=True)
    ap.add_argument("--mutation-counts", required=True)
    ap.add_argument("--in-list", required=True)
    ap.add_argument("--min-cov", type=float, default=0.5)
    ap.add_argument("--min-id", type=float, default=0.4)
    ap.add_argument("--min-aln-len", type=int, default=100)
    ap.add_argument("--out-list", required=True)
    ap.add_argument("--out-report", required=True)
    args = ap.parse_args()

    lengths = parse_fasta_lengths(args.proteins)
    alignments = parse_alignment(args.protein_alignment)
    counts = parse_mutation_counts(args.mutation_counts)
    with open(args.in_list, "r") as fh:
        multi_ids = [line.strip() for line in fh if line.strip()]

    with open(args.out_list, "w") as outp, open(args.out_report, "w") as rep:
        rep.write("transcript_id\tstop\tshift\ttotal\talign_len\tprotein_len\tcoverage\tidentity\treason\n")
        for tid in multi_ids:
            stop, shift, total = counts.get(tid, (0, 0, 0))
            if tid not in alignments:
                rep.write(f"{tid}\t{stop}\t{shift}\t{total}\t0\t0\t0\t0\tmissing_alignment\n")
                outp.write(f"{tid}\n")
                continue
            if tid not in lengths:
                rep.write(f"{tid}\t{stop}\t{shift}\t{total}\t0\t0\t0\t0\tmissing_length\n")
                outp.write(f"{tid}\n")
                continue

            qseq, tseq = alignments[tid]
            aln_len = min(len(qseq), len(tseq))
            prot_len = lengths[tid]
            cov = aln_len / max(1, prot_len)
            ident = calc_identity(qseq, tseq)

            reason = []
            if aln_len < args.min_aln_len:
                reason.append("short_alignment")
            if cov < args.min_cov:
                reason.append("low_coverage")
            if ident < args.min_id:
                reason.append("low_identity")

            if reason:
                outp.write(f"{tid}\n")
                rep.write(f"{tid}\t{stop}\t{shift}\t{total}\t{aln_len}\t{prot_len}\t{cov:.4f}\t{ident:.4f}\t{';'.join(reason)}\n")
            else:
                rep.write(f"{tid}\t{stop}\t{shift}\t{total}\t{aln_len}\t{prot_len}\t{cov:.4f}\t{ident:.4f}\tpass\n")


if __name__ == "__main__":
    main()
