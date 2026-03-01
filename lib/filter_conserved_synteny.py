#!/usr/bin/env python3
import argparse


def load_synteny(path):
    synteny = {}
    with open(path, "r") as fh:
        for line in fh:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            gid = parts[0]
            chr_ = parts[1]
            start = int(parts[2])
            end = int(parts[3])
            if end < start:
                start, end = end, start
            synteny[gid] = (chr_, start, end)
    return synteny


def overlap_ratio(a_start, a_end, b_start, b_end):
    overlap = max(0, min(a_end, b_end) - max(a_start, b_start))
    if overlap <= 0:
        return 0.0, 0.0
    a_len = max(1, a_end - a_start)
    b_len = max(1, b_end - b_start)
    return overlap / a_len, overlap / b_len


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--synteny-tsv", required=True)
    ap.add_argument("--best-hits", required=True)
    ap.add_argument("--out-pass", required=True)
    ap.add_argument("--out-fail", required=True)
    ap.add_argument("--min-overlap", type=float, default=0.1)
    args = ap.parse_args()

    syn = load_synteny(args.synteny_tsv)
    pass_ids = set()
    fail_ids = set()

    with open(args.best_hits, "r") as fh:
        for line in fh:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            gid = parts[0]
            tchr = parts[2]
            tstart = int(parts[3])
            tend = int(parts[4])
            if tend < tstart:
                tstart, tend = tend, tstart

            if gid not in syn:
                fail_ids.add(gid)
                continue
            schr, sstart, send = syn[gid]
            if schr != tchr:
                fail_ids.add(gid)
                continue
            r1, r2 = overlap_ratio(tstart, tend, sstart, send)
            if r1 >= args.min_overlap and r2 >= args.min_overlap:
                pass_ids.add(gid)
            else:
                fail_ids.add(gid)

    with open(args.out_pass, "w") as outp:
        for gid in sorted(pass_ids):
            outp.write(f"{gid}\n")
    with open(args.out_fail, "w") as outp:
        for gid in sorted(fail_ids):
            if gid not in pass_ids:
                outp.write(f"{gid}\n")


if __name__ == "__main__":
    main()
