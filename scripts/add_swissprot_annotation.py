#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
add_swissprot_annotation.py
===========================
Attach the best Swiss-Prot BLASTP hit to each protein FASTA header.

The BLAST table must use this column order:

    qseqid sseqid pident length evalue bitscore stitle

For example:

    blastp -query protein.faa -db uniprot_sprot \
        -evalue 1e-5 -max_target_seqs 5 -max_hsps 1 -num_threads 50 \
        -outfmt "6 qseqid sseqid pident length evalue bitscore stitle" \
        -out protein.blastp.tsv

For each query, the hit with the lowest E-value is retained; ties are broken
by higher bit score and then by higher percent identity. If several hits share
the same E-value but differ in score, the highest-scoring hit is kept.

Usage:
    python3 scripts/add_swissprot_annotation.py \
        --query protein.faa \
        --blast protein.blastp.tsv \
        --out protein.swissprot.faa \
        --summary protein.swissprot.tsv
"""
from __future__ import annotations

import argparse
import gzip
import sys
from collections import Counter
from pathlib import Path


def open_text(path: Path):
    # utf-8-sig removes a UTF-8 BOM if present (common in Windows-edited files)
    # and behaves like utf-8 otherwise.
    if path.suffix.lower() == ".gz":
        return gzip.open(path, "rt", encoding="utf-8-sig", errors="replace")
    return open(path, "r", encoding="utf-8-sig", errors="replace")


def parse_blast(path: Path, max_evalue: float, min_pident: float) -> dict[str, tuple]:
    """Return query -> (sseqid, pident, length, evalue, bitscore, stitle)."""
    best: dict[str, tuple] = {}
    with open_text(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t", 6)
            if len(fields) < 7:
                continue
            qseqid, sseqid, pident_s, length_s, evalue_s, bitscore_s, stitle = fields
            stitle = stitle.replace("\t", " ").strip()
            try:
                pident = float(pident_s)
                length = int(float(length_s))
                evalue = float(evalue_s)
                bitscore = float(bitscore_s)
            except ValueError:
                continue
            if evalue > max_evalue or pident < min_pident:
                continue
            candidate = (evalue, -bitscore, -pident)
            previous = best.get(qseqid)
            if previous is None or candidate < previous[0]:
                best[qseqid] = (candidate, sseqid, pident, length, evalue, bitscore, stitle)
    return {q: value[1:] for q, value in best.items()}


def read_fasta(path: Path):
    records = []
    header = None
    seq = []
    with open_text(path) as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                if header is not None and header.strip():
                    records.append((header, "".join(seq)))
                header = line[1:].strip()
                seq = []
            elif header is not None:
                seq.append(line.strip())
    if header is not None and header.strip():
        records.append((header, "".join(seq)))
    return records


def wrap(seq: str, width: int = 60) -> str:
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--query", type=Path, required=True, help="input protein FASTA")
    ap.add_argument("--blast", type=Path, required=True, help="BLASTP tabular output")
    ap.add_argument("--out", type=Path, required=True, help="annotated output FASTA")
    ap.add_argument("--summary", type=Path, default=None, help="optional TSV summary")
    ap.add_argument("--max-evalue", type=float, default=1e-5,
                    help="maximum E-value retained (default: 1e-5)")
    ap.add_argument("--min-pident", type=float, default=0.0,
                    help="minimum percent identity retained (default: 0)")
    ap.add_argument("--mark-missing", action="store_true",
                    help="add swissprot=no_hit to headers without a retained hit")
    args = ap.parse_args()

    if args.max_evalue < 0:
        raise SystemExit("ERROR: --max-evalue must be non-negative")
    if args.min_pident < 0 or args.min_pident > 100:
        raise SystemExit("ERROR: --min-pident must be between 0 and 100")

    for path in (args.query, args.blast):
        if not path.exists():
            raise SystemExit(f"ERROR: input file does not exist: {path}")

    hits = parse_blast(args.blast, args.max_evalue, args.min_pident)
    records = read_fasta(args.query)
    if not records:
        raise SystemExit(f"ERROR: no FASTA records found in {args.query}")

    query_ids = [header.split()[0] for header, _seq in records]
    counts = Counter(query_ids)
    duplicated = sorted(qid for qid, n in counts.items() if n > 1)
    if duplicated:
        preview = ", ".join(duplicated[:10])
        raise SystemExit(f"ERROR: duplicate FASTA sequence IDs: {preview}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    n_annotated = 0
    with open(args.out, "w", encoding="utf-8", newline="\n") as out:
        for header, seq in records:
            qid = header.split()[0]
            hit = hits.get(qid)
            if hit is None:
                new_header = f"{header} swissprot=no_hit" if args.mark_missing else header
            else:
                sseqid, pident, length, evalue, bitscore, stitle = hit
                new_header = (f"{header} swissprot={sseqid} pident={pident:g} "
                              f"length={length} evalue={evalue:g} bitscore={bitscore:g} {stitle}")
                n_annotated += 1
            out.write(f">{new_header}\n{wrap(seq)}\n")

    if args.summary is not None:
        args.summary.parent.mkdir(parents=True, exist_ok=True)
        with open(args.summary, "w", encoding="utf-8", newline="\n") as out:
            out.write("query\tseq_id\tpident\tlength\tevalue\tbitscore\tdescription\n")
            for header, _seq in records:
                qid = header.split()[0]
                hit = hits.get(qid)
                if hit is None:
                    out.write(f"{qid}\t\t\t\t\t\t\n")
                else:
                    sseqid, pident, length, evalue, bitscore, stitle = hit
                    out.write(f"{qid}\t{sseqid}\t{pident:g}\t{length}\t{evalue:g}\t"
                              f"{bitscore:g}\t{stitle}\n")

    print(f"[swissprot] {n_annotated}/{len(records)} proteins annotated -> {args.out}",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
