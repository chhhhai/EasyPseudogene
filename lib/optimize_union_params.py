#!/usr/bin/env python3
import argparse
from pathlib import Path


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


def identity(a, b):
    if not a or not b:
        return 0.0
    L = min(len(a), len(b))
    if L == 0:
        return 0.0
    return sum(1 for i in range(L) if a[i] == b[i]) / L


def flank_identity(q, t, pos, flank):
    q_b = q[max(0, pos - flank) : pos]
    t_b = t[max(0, pos - flank) : pos]
    q_a = q[pos + 1 : pos + 1 + flank]
    t_a = t[pos + 1 : pos + 1 + flank]
    return identity(q_b, t_b), identity(q_a, t_a)


def load_set(path):
    if not path.exists():
        return set()
    return set([l.strip() for l in path.read_text().splitlines() if l.strip()])


def compute_features(unitary, one_mut, multi_mut, alignments, lengths, flank_size):
    feats = {}
    for tid in unitary:
        q, t = alignments.get(tid, ("", ""))
        aln_len = min(len(q), len(t))
        prot_len = lengths.get(tid, 0)
        cov = aln_len / max(1, prot_len) if prot_len else 0.0
        ident = identity(q, t) if aln_len else 0.0

        pos_x = t.find("X")
        pos_b = t.find("!")
        positions = [(pos_x, "stop"), (pos_b, "shift")]
        positions = [p for p in positions if p[0] >= 0]
        if positions:
            pos, mtype = min(positions, key=lambda x: x[0])
            bef, aft = flank_identity(q, t, pos, flank_size)
            frac = (pos + 1) / max(1, prot_len) if prot_len else 0.0
        else:
            pos = -1
            mtype = ""
            bef = 1.0
            aft = 1.0
            frac = 0.0

        feats[tid] = {
            "aln_len": aln_len,
            "cov": cov,
            "ident": ident,
            "mut_pos": pos,
            "mut_type": mtype,
            "flank_bef": bef,
            "flank_aft": aft,
            "mut_frac": frac,
            "is_one": tid in one_mut,
            "is_multi": tid in multi_mut,
        }
    return feats


def apply_params(unitary, feats, splicing_sets, params):
    low_mode = params["low_mode"]
    low_flank = params["low_flank"]
    low_overall = params["low_overall"]
    end_frac = params["end_frac"]
    multi_aln = params["multi_aln"]
    multi_cov = params["multi_cov"]
    multi_id = params["multi_id"]
    glob = params.get("global")

    false_pos = set()
    false_pos |= splicing_sets

    for tid in unitary:
        f = feats.get(tid)
        if f is None:
            continue
        if f["is_one"]:
            if f["mut_pos"] >= 0:
                if low_mode == "or":
                    low_flag = (f["flank_bef"] < low_flank or f["flank_aft"] < low_flank)
                else:
                    low_flag = (f["flank_bef"] < low_flank and f["flank_aft"] < low_flank)
                if low_flag and f["ident"] < low_overall:
                    false_pos.add(tid)
            if f["mut_frac"] >= end_frac:
                false_pos.add(tid)

    high_conf = set([t for t in unitary if t not in false_pos])

    # Multi-disruption quality filter
    low_quality_multi = set()
    for tid in list(high_conf):
        f = feats.get(tid)
        if not f or not f["is_multi"]:
            continue
        if not (f["aln_len"] >= multi_aln and f["cov"] >= multi_cov and f["ident"] >= multi_id):
            low_quality_multi.add(tid)
    high_conf = high_conf - low_quality_multi

    # Global alignment filter
    if glob is not None:
        g_aln, g_cov, g_id = glob
        high_conf = set(
            t for t in high_conf
            if feats.get(t) and feats[t]["aln_len"] >= g_aln
            and feats[t]["cov"] >= g_cov and feats[t]["ident"] >= g_id
        )

    return high_conf


def parse_list_arg(arg, cast=float):
    return [cast(x) for x in arg.split(",") if x.strip()]


def parse_global_arg(arg):
    options = []
    if not arg:
        return [None]
    for token in arg.split(","):
        token = token.strip()
        if not token:
            continue
        if token.lower() == "off":
            options.append(None)
            continue
        parts = token.split(":")
        if len(parts) != 3:
            continue
        options.append((int(parts[0]), float(parts[1]), float(parts[2])))
    # de-duplicate
    seen = set()
    uniq = []
    for opt in options:
        key = "off" if opt is None else opt
        if key in seen:
            continue
        seen.add(key)
        uniq.append(opt)
    return uniq


def main():
    ap = argparse.ArgumentParser(description="Optimize union vs ref overlap by parameter grid search.")
    ap.add_argument("--base-dir", required=True)
    ap.add_argument("--proteins", required=True)
    ap.add_argument("--ref-list", required=True)
    ap.add_argument("--species", default="dolphin,sperm_whale")
    ap.add_argument("--flank-size", type=int, default=30)
    ap.add_argument("--low-mode", default="or,and")
    ap.add_argument("--low-flank", default="0.3,0.4")
    ap.add_argument("--low-overall", default="0.2,1.0")
    ap.add_argument("--end-frac", default="0.9,0.95")
    ap.add_argument("--multi-aln", default="200,300")
    ap.add_argument("--multi-cov", default="0.6,0.7")
    ap.add_argument("--multi-id", default="0.05,0.1")
    ap.add_argument("--global-filter", default="off,150:0.8:0.3,200:0.8:0.3")
    ap.add_argument("--top", type=int, default=20)
    ap.add_argument("--max-union", type=int, default=0, help="Only report combinations with union <= this value")
    ap.add_argument("--target-union", type=int, default=0, help="Sort by closeness to this union size")
    ap.add_argument("--min-recall", type=float, default=0.0, help="Only report combinations with recall >= this value")
    args = ap.parse_args()

    base = Path(args.base_dir)
    ref = set([l.strip() for l in Path(args.ref_list).read_text().splitlines() if l.strip()])
    lengths = parse_fasta_lengths(args.proteins)

    species_list = [s.strip() for s in args.species.split(",") if s.strip()]

    data = {}
    for sp in species_list:
        sp_dir = base / sp / "post_genewise_filter"
        unitary = load_set(sp_dir / "unitary_candidates.txt")
        one_mut = load_set(sp_dir / "one_mutation_list.txt")
        multi_mut = load_set(sp_dir / "multi_mutation_for_manual_review.txt")
        splicing_sets = set()
        for name in ["splicing-stop.txt","ending-stop.txt","splicing-shift.txt","ending-shift.txt"]:
            splicing_sets |= load_set(sp_dir / name)
        alignments = parse_alignment(sp_dir / "protein-alignment.txt")
        feats = compute_features(unitary, one_mut, multi_mut, alignments, lengths, args.flank_size)
        data[sp] = (unitary, feats, splicing_sets)

    low_mode_list = [m.strip() for m in args.low_mode.split(",") if m.strip()]
    low_flank_list = parse_list_arg(args.low_flank, float)
    low_overall_list = parse_list_arg(args.low_overall, float)
    end_frac_list = parse_list_arg(args.end_frac, float)
    multi_aln_list = parse_list_arg(args.multi_aln, int)
    multi_cov_list = parse_list_arg(args.multi_cov, float)
    multi_id_list = parse_list_arg(args.multi_id, float)
    global_list = parse_global_arg(args.global_filter)

    rows = []
    for low_mode in low_mode_list:
        for low_flank in low_flank_list:
            for low_overall in low_overall_list:
                for end_frac in end_frac_list:
                    for multi_aln in multi_aln_list:
                        for multi_cov in multi_cov_list:
                            for multi_id in multi_id_list:
                                for glob in global_list:
                                    params = {
                                        "low_mode": low_mode,
                                        "low_flank": low_flank,
                                        "low_overall": low_overall,
                                        "end_frac": end_frac,
                                        "multi_aln": multi_aln,
                                        "multi_cov": multi_cov,
                                        "multi_id": multi_id,
                                        "global": glob,
                                    }
                                    final_sets = {}
                                    for sp in species_list:
                                        unitary, feats, splicing_sets = data[sp]
                                        final_sets[sp] = apply_params(unitary, feats, splicing_sets, params)
                                    union = set().union(*final_sets.values())
                                    overlap = union & ref
                                    recall = len(overlap) / len(ref) if ref else 0.0
                                    precision = len(overlap) / len(union) if union else 0.0
                                    rows.append({
                                        "low_mode": low_mode,
                                        "low_flank": low_flank,
                                        "low_overall": low_overall,
                                        "end_frac": end_frac,
                                        "multi_aln": multi_aln,
                                        "multi_cov": multi_cov,
                                        "multi_id": multi_id,
                                        "global": glob,
                                        "union": len(union),
                                        "overlap": len(overlap),
                                        "recall": recall,
                                        "precision": precision,
                                        "dolphin": len(final_sets.get("dolphin", set())),
                                        "sperm_whale": len(final_sets.get("sperm_whale", set())),
                                    })

    if args.max_union > 0:
        rows = [r for r in rows if r["union"] <= args.max_union]
    if args.min_recall > 0:
        rows = [r for r in rows if r["recall"] >= args.min_recall]

    if args.target_union > 0:
        rows.sort(key=lambda r: (abs(r["union"] - args.target_union), -r["recall"], -r["precision"]))
    elif args.min_recall > 0:
        rows.sort(key=lambda r: (r["union"], -r["precision"], -r["recall"]))
    else:
        # Sort: maximize recall, then minimize union
        rows.sort(key=lambda r: (-r["recall"], r["union"], -r["precision"]))

    print("low_mode\tlow_flank\tlow_overall\tend_frac\tmulti_aln\tmulti_cov\tmulti_id\tglobal\tunion\toverlap\trecall\tprecision\tdolphin\tsperm_whale")
    for r in rows[: args.top]:
        glob = "off" if r["global"] is None else f"{r['global'][0]}:{r['global'][1]}:{r['global'][2]}"
        print(
            f"{r['low_mode']}\t{r['low_flank']}\t{r['low_overall']}\t{r['end_frac']}\t"
            f"{r['multi_aln']}\t{r['multi_cov']}\t{r['multi_id']}\t{glob}\t"
            f"{r['union']}\t{r['overlap']}\t{r['recall']:.3f}\t{r['precision']:.3f}\t"
            f"{r['dolphin']}\t{r['sperm_whale']}"
        )


if __name__ == "__main__":
    main()
