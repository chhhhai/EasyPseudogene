#!/usr/bin/env python3
import argparse


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


def flank_identities(q, t, pos, flank):
    q_b = q[max(0, pos - flank) : pos]
    t_b = t[max(0, pos - flank) : pos]
    q_a = q[pos + 1 : pos + 1 + flank]
    t_a = t[pos + 1 : pos + 1 + flank]
    return identity(q_b, t_b), identity(q_a, t_a)


def main():
    ap = argparse.ArgumentParser(description="Filter single-mutation candidates with low-identity flanks.")
    ap.add_argument("--protein-alignment", required=True)
    ap.add_argument("--one-mutation-list", required=True)
    ap.add_argument("--out-stop", required=True)
    ap.add_argument("--out-shift", required=True)
    ap.add_argument("--flank-size", type=int, default=30)
    ap.add_argument("--flank-th", type=float, default=0.2)
    ap.add_argument("--overall-th", type=float, default=0.1)
    ap.add_argument("--mode", choices=["and", "or"], default="and")
    args = ap.parse_args()

    one_mut = set()
    with open(args.one_mutation_list, "r") as fh:
        for line in fh:
            line = line.strip()
            if line:
                one_mut.add(line.split("\t")[0])

    aln = parse_alignment(args.protein_alignment)

    stop_ids = set()
    shift_ids = set()

    for tid in one_mut:
        if tid not in aln:
            continue
        q, t = aln[tid]
        pos_x = t.find("X")
        pos_b = t.find("!")
        positions = [(pos_x, "stop"), (pos_b, "shift")]
        positions = [p for p in positions if p[0] >= 0]
        if not positions:
            continue
        pos, mtype = min(positions, key=lambda x: x[0])

        bef, aft = flank_identities(q, t, pos, args.flank_size)
        overall = identity(q, t)

        if args.mode == "and":
            flag = (bef < args.flank_th and aft < args.flank_th and overall < args.overall_th)
        else:
            flag = ((bef < args.flank_th or aft < args.flank_th) and overall < args.overall_th)

        if flag:
            if mtype == "stop":
                stop_ids.add(tid)
            else:
                shift_ids.add(tid)

    with open(args.out_stop, "w") as out:
        for tid in sorted(stop_ids):
            out.write(f"{tid}\n")
    with open(args.out_shift, "w") as out:
        for tid in sorted(shift_ids):
            out.write(f"{tid}\n")


if __name__ == "__main__":
    main()
