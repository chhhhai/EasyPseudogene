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


def parse_list(path):
    if not path:
        return set()
    try:
        return set([line.strip() for line in open(path) if line.strip()])
    except FileNotFoundError:
        return set()


def grid_from_arg(arg, default):
    if not arg:
        return default
    vals = []
    for p in arg.split(","):
        p = p.strip()
        if not p:
            continue
        vals.append(float(p))
    return vals if vals else default


def main():
    ap = argparse.ArgumentParser(description="Compare overlap with reference list and suggest thresholds.")
    ap.add_argument("--ref-list", required=True)
    ap.add_argument("--filter-dir", required=True)
    ap.add_argument("--proteins", required=True)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--cov-grid", default="")
    ap.add_argument("--id-grid", default="")
    ap.add_argument("--min-aln-len-grid", default="50,100,150,200")
    args = ap.parse_args()

    ref = parse_list(args.ref_list)
    lengths = parse_fasta_lengths(args.proteins)
    alignments = parse_alignment(f"{args.filter_dir}/protein-alignment.txt")
    counts = parse_mutation_counts(f"{args.filter_dir}/mutation_counts.tsv")

    multi_list = parse_list(f"{args.filter_dir}/multi_mutation_for_manual_review.txt")
    unitary = parse_list(f"{args.filter_dir}/unitary_candidates.txt")

    cov_grid = grid_from_arg(args.cov_grid, [0.3, 0.4, 0.5, 0.6, 0.7, 0.8])
    id_grid = grid_from_arg(args.id_grid, [0.3, 0.4, 0.5, 0.6, 0.7])
    aln_grid = grid_from_arg(args.min_aln_len_grid, [50, 100, 150, 200])

    feature_path = f"{args.out_prefix}.features.tsv"
    with open(feature_path, "w") as out:
        out.write("transcript_id\tin_ref\tin_unitary\tstop\tshift\ttotal\talign_len\tprotein_len\tcoverage\tidentity\tmulti_disruption\n")
        for tid in sorted(unitary):
            stop, shift, total = counts.get(tid, (0, 0, 0))
            qseq, tseq = alignments.get(tid, ("", ""))
            aln_len = min(len(qseq), len(tseq))
            prot_len = lengths.get(tid, 0)
            cov = aln_len / max(1, prot_len) if prot_len else 0.0
            ident = calc_identity(qseq, tseq) if aln_len else 0.0
            multi = 1 if tid in multi_list else 0
            out.write(f"{tid}\t{int(tid in ref)}\t1\t{stop}\t{shift}\t{total}\t{aln_len}\t{prot_len}\t{cov:.4f}\t{ident:.4f}\t{multi}\n")

    # Suggest thresholds for multi-disruption candidates
    multi_ref = [tid for tid in multi_list if tid in ref]
    multi_nonref = [tid for tid in multi_list if tid not in ref]

    def get_metrics(tid):
        qseq, tseq = alignments.get(tid, ("", ""))
        aln_len = min(len(qseq), len(tseq))
        prot_len = lengths.get(tid, 0)
        cov = aln_len / max(1, prot_len) if prot_len else 0.0
        ident = calc_identity(qseq, tseq) if aln_len else 0.0
        return aln_len, cov, ident

    sugg_path = f"{args.out_prefix}.threshold_suggestions.tsv"
    with open(sugg_path, "w") as out:
        out.write("min_aln_len\tmin_cov\tmin_id\tkeep_ref\tdrop_nonref\tref_kept\tref_total\tnonref_dropped\tnonref_total\tscore\n")
        for min_aln in aln_grid:
            for min_cov in cov_grid:
                for min_id in id_grid:
                    ref_total = len(multi_ref)
                    nonref_total = len(multi_nonref)
                    ref_kept = 0
                    nonref_dropped = 0
                    for tid in multi_ref:
                        aln_len, cov, ident = get_metrics(tid)
                        if aln_len >= min_aln and cov >= min_cov and ident >= min_id:
                            ref_kept += 1
                    for tid in multi_nonref:
                        aln_len, cov, ident = get_metrics(tid)
                        if not (aln_len >= min_aln and cov >= min_cov and ident >= min_id):
                            nonref_dropped += 1
                    keep_ref = ref_kept / ref_total if ref_total else 0.0
                    drop_nonref = nonref_dropped / nonref_total if nonref_total else 0.0
                    score = keep_ref + drop_nonref
                    out.write(f"{min_aln}\t{min_cov:.2f}\t{min_id:.2f}\t{keep_ref:.3f}\t{drop_nonref:.3f}\t{ref_kept}\t{ref_total}\t{nonref_dropped}\t{nonref_total}\t{score:.3f}\n")


if __name__ == "__main__":
    main()
