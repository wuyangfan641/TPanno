#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Convert a CD-HIT .clstr file into a representative<TAB>member mapping."""
import argparse
import re
import sys


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--clstr", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    rows = []
    cluster = []
    rep = None

    def flush():
        # The representative is marked with '*' and may appear on any line of
        # the cluster, not necessarily the first one.
        if cluster:
            representative = rep or cluster[0]
            for member in cluster:
                rows.append((representative, member))

    with open(args.clstr, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">Cluster"):
                flush()
                cluster = []
                rep = None
                continue
            m = re.search(r">(\S+?)\.\.\.", line)
            if not m:
                continue
            sid = m.group(1)
            cluster.append(sid)
            if line.rstrip().endswith("*"):
                rep = sid
        flush()

    with open(args.out, "w", encoding="utf-8", newline="\n") as out:
        out.write("representative\tmember\n")
        for r, s in rows:
            out.write(f"{r}\t{s}\n")
    print(f"[cdhit] {len(rows)} sequences -> {len({r for r, _ in rows})} representatives",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
