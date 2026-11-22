#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
00_prepare_inputs.py
====================
Scan the protein FASTA files in 00_raw/, perform input quality control (QC)
and normalization, and write:

  results/00_prep/clean_fasta/<sample>.faa   cleaned protein sequences (used downstream)
  results/00_prep/samples.clean.tsv          sample -> cleaned FASTA mapping
  results/00_prep/input_qc.tsv               per-sample QC metrics
  results/00_prep/input_qc.md                human-readable QC report

Processing rules (all changes are recorded explicitly, never silently):
  1. The sequence ID is the first whitespace-delimited field of the header.
  2. Sequences are uppercased; whitespace, '-' and '.' are removed.
  3. Terminal '*' characters are removed; internal '*' are converted to X by
     default (use --drop-internal-stop to discard such sequences instead).
  4. U/O are converted to X by default (HMMER's 20-letter models do not model
     these rare amino acids).
  5. Other non-standard characters are converted to X and counted.
  6. Sequences shorter than --min-len or with an X fraction above --max-x are
     discarded and counted.
  7. Duplicate IDs abort the run by default (to avoid ambiguity in the
     downstream sample|protein mapping).

Usage:
  python bin/00_prepare_inputs.py --raw-dir 00_raw --outdir results/00_prep
"""
from __future__ import annotations

import argparse
import gzip
import re
import sys
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

try:
    from Bio import SeqIO
except ImportError:  # pragma: no cover
    sys.exit("ERROR: Biopython is required. Run `conda env create -f env/environment.yml` first.")


FASTA_SUFFIXES = (
    ".pep.fasta", ".pep.fa", ".pep.faa",
    ".protein.fasta", ".protein.fa", ".protein.faa",
    ".fasta", ".faa", ".fa", ".fas",
)
VALID_AA = set("ACDEFGHIKLMNPQRSTVWYBZXJ")


@dataclass
class SampleQC:
    sample: str
    raw_fasta: str
    n_raw: int = 0
    n_clean: int = 0
    n_short: int = 0
    n_high_x: int = 0
    n_internal_stop: int = 0
    n_stop_dropped: int = 0
    n_uo: int = 0
    n_other: int = 0
    n_dup: int = 0
    total_aa: int = 0
    min_len: int = 0
    max_len: int = 0
    clean_fasta: str = ""
    warnings: list[str] = field(default_factory=list)


def open_maybe_gzip(path: Path):
    # utf-8-sig strips a UTF-8 BOM if one was added by a Windows editor.
    if path.suffix.lower() == ".gz":
        return gzip.open(path, "rt", encoding="utf-8-sig", errors="replace")
    return open(path, "r", encoding="utf-8-sig", errors="replace")


def infer_sample_name(path: Path) -> str:
    name = path.name
    if name.lower().endswith(".gz"):
        name = name[:-3]
    low = name.lower()
    for suf in FASTA_SUFFIXES:  # longest suffix first
        if low.endswith(suf):
            name = name[: -len(suf)]
            break
    else:
        name = path.stem
    name = re.sub(r"[^A-Za-z0-9_.-]+", "_", name).strip("_")
    if not name:
        raise ValueError(f"cannot infer a sample name from file name: {path}")
    return name


def discover_samples(raw_dir: Path) -> list[tuple[str, Path]]:
    files = []
    for p in sorted(raw_dir.iterdir()):
        if not p.is_file():
            continue
        low = p.name.lower()
        if low.endswith(".gz"):
            low = low[:-3]
        if low.endswith(FASTA_SUFFIXES):
            files.append(p)
    if not files:
        raise SystemExit(f"ERROR: no FASTA files found in {raw_dir}")
    seen: dict[str, Path] = {}
    out = []
    for p in files:
        s = infer_sample_name(p)
        if s in seen:
            raise SystemExit(
                f"ERROR: sample name collision: {seen[s].name} and {p.name} both map to '{s}'. "
                "Specify the mapping manually in config/samples.tsv."
            )
        seen[s] = p
        out.append((s, p))
    return out


def read_samples_tsv(path: Path, raw_dir: Path) -> list[tuple[str, Path]]:
    rows = []
    with open(path, "r", encoding="utf-8-sig") as fh:
        header = None
        for line in fh:
            line = line.rstrip("\r\n")
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.split("\t")
            if header is None:
                header = [x.strip().lower() for x in parts]
                for req in ("sample", "fasta"):
                    if req not in header:
                        raise SystemExit(f"ERROR: {path} is missing the column '{req}' (header: {header})")
                continue
            rec = dict(zip(header, [x.strip() for x in parts]))
            if "sample" not in rec or "fasta" not in rec:
                raise SystemExit(f"ERROR: malformed row in {path}: {line.rstrip()}")
            p = Path(rec["fasta"])
            if not p.is_absolute():
                p = raw_dir / p
            rows.append((rec["sample"], p))
    if not rows:
        raise SystemExit(f"ERROR: no sample records found in {path}")
    return rows


def check_unique_ids(ids: list[str], sample: str, path: Path) -> None:
    dup = [k for k, v in Counter(ids).items() if v > 1]
    if dup:
        preview = ", ".join(dup[:10])
        raise SystemExit(
            f"ERROR: sample {sample} ({path.name}) contains {len(dup)} duplicate protein IDs, "
            f"for example: {preview}\n"
            "       Add a unique suffix to the IDs, or split the file via config/samples.tsv.\n"
            "       Downstream code uses sample|protein as the unique key, and duplicate IDs "
            "would cause domain misassignment."
        )


def clean_sequence(raw: str, args) -> tuple[str, dict]:
    seq = re.sub(r"\s+", "", raw).upper()
    n_dash = seq.count("-") + seq.count(".")
    seq = seq.replace("-", "").replace(".", "")

    n_trailing_stop = len(seq) - len(seq.rstrip("*"))
    n_internal_stop = seq.count("*") - n_trailing_stop
    seq = seq.rstrip("*")
    stop_dropped = False
    if n_internal_stop:
        if args.drop_internal_stop:
            stop_dropped = True
        else:
            seq = seq.replace("*", "X")

    n_uo = seq.count("U") + seq.count("O")
    if args.uo_to_x:
        seq = seq.replace("U", "X").replace("O", "X")

    n_other = sum(1 for c in seq if c not in VALID_AA and c not in ("U", "O"))
    if n_other:
        seq = "".join(c if c in VALID_AA else "X" for c in seq)

    return seq, {
        "n_dash": n_dash,
        "n_internal_stop": n_internal_stop,
        "stop_dropped": stop_dropped,
        "n_uo": n_uo,
        "n_other": n_other,
    }


def wrap(seq: str, width: int = 60) -> str:
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def process_sample(sample: str, path: Path, out_faa: Path, args, write_output: bool = True) -> SampleQC:
    qc = SampleQC(sample=sample, raw_fasta=str(path))
    ids: list[str] = []
    records: list[tuple[str, str]] = []

    with open_maybe_gzip(path) as fh:
        for rec in SeqIO.parse(fh, "fasta"):
            qc.n_raw += 1
            rid = rec.id
            if not rid:
                raise SystemExit(f"ERROR: {path} contains a record with an empty ID")
            if re.search(r"\s", rid):
                rid = rid.split()[0]
            ids.append(rid)

            seq, info = clean_sequence(str(rec.seq), args)
            qc.n_internal_stop += 1 if info["n_internal_stop"] else 0
            qc.n_uo += info["n_uo"]
            qc.n_other += info["n_other"]

            if info["stop_dropped"]:
                qc.n_stop_dropped += 1
                continue
            if not seq or len(seq) < args.min_len:
                qc.n_short += 1
                continue
            x_frac = seq.count("X") / len(seq)
            if x_frac > args.max_x:
                qc.n_high_x += 1
                continue

            records.append((rid, seq))
            qc.total_aa += len(seq)

    qc.n_dup = sum(1 for _, n in Counter(ids).items() if n > 1)
    check_unique_ids(ids, sample, path)

    qc.n_clean = len(records)
    if qc.n_clean == 0:
        raise SystemExit(f"ERROR: sample {sample} has no sequences left after cleaning; check {path}")
    lengths = [len(s) for _, s in records]
    qc.min_len, qc.max_len = min(lengths), max(lengths)

    if qc.n_short:
        qc.warnings.append(f"{qc.n_short} sequences shorter than {args.min_len} aa were discarded")
    if qc.n_high_x:
        qc.warnings.append(f"{qc.n_high_x} sequences with X fraction > {args.max_x} were discarded")
    if qc.n_internal_stop:
        if args.drop_internal_stop:
            qc.warnings.append(f"{qc.n_internal_stop} sequences contained internal stop codons and were discarded")
        else:
            qc.warnings.append(f"{qc.n_internal_stop} sequences contained internal stop codons (converted to X)")
    if qc.n_uo:
        if args.uo_to_x:
            qc.warnings.append(f"{qc.n_uo} U/O residues were converted to X")
        else:
            qc.warnings.append(f"{qc.n_uo} U/O residues were retained")

    qc.clean_fasta = str(out_faa)
    if write_output:
        out_faa.parent.mkdir(parents=True, exist_ok=True)
        with open(out_faa, "w", encoding="utf-8", newline="\n") as out:
            for rid, seq in records:
                out.write(f">{rid}\n{wrap(seq)}\n")
    return qc


def write_markdown(qcs: list[SampleQC], out_md: Path, args) -> None:
    total_raw = sum(q.n_raw for q in qcs)
    total_clean = sum(q.n_clean for q in qcs)
    lines = [
        "# Input data QC report",
        "",
        f"- Samples: **{len(qcs)}**",
        f"- Raw sequences: **{total_raw:,}**",
        f"- Sequences after cleaning: **{total_clean:,}**",
        f"- Cleaning rules: minimum length {args.min_len} aa; X fraction limit {args.max_x:.0%}; "
        f"U/O to X = {'yes' if args.uo_to_x else 'no'}; internal stop to X = "
        f"{'no (discard)' if args.drop_internal_stop else 'yes'}",
        "",
        "## Per-sample statistics",
        "",
        "| Sample | Raw | Clean | Discarded (short) | Discarded (high X) | Internal stops | U/O | Other non-standard | Mean length | Range |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---|",
    ]
    for q in qcs:
        mean_len = q.total_aa / q.n_clean if q.n_clean else 0
        lines.append(
            f"| {q.sample} | {q.n_raw:,} | {q.n_clean:,} | {q.n_short:,} | {q.n_high_x:,} | "
            f"{q.n_internal_stop:,} | {q.n_uo:,} | {q.n_other:,} | {mean_len:.1f} | "
            f"{q.min_len}-{q.max_len} |"
        )
    warns = [(q.sample, w) for q in qcs for w in q.warnings]
    if warns:
        lines += ["", "## Warnings", ""]
        lines += [f"- **{s}**: {w}" for s, w in warns]
    out_md.parent.mkdir(parents=True, exist_ok=True)
    out_md.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--raw-dir", type=Path, default=Path("00_raw"))
    ap.add_argument("--samples", type=Path, default=None,
                    help="optional samples.tsv; if omitted, the raw directory is scanned")
    ap.add_argument("--outdir", type=Path, default=Path("results/00_prep"))
    ap.add_argument("--write-samples", type=Path, default=None,
                    help="write the discovered sample -> FASTA mapping as a samples.tsv template")
    ap.add_argument("--min-len", type=int,
                    default=int(__import__("os").environ.get("MIN_PROTEIN_LEN", 30)))
    ap.add_argument("--max-x", type=float,
                    default=float(__import__("os").environ.get("MAX_X_FRAC", 0.5)))
    ap.add_argument("--uo-to-x", action=argparse.BooleanOptionalAction,
                    default=__import__("os").environ.get("CLEAN_UO_TO_X", "1") == "1")
    ap.add_argument("--drop-internal-stop", action="store_true",
                    default=__import__("os").environ.get("DROP_INTERNAL_STOP", "0") == "1")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    if not args.raw_dir.is_dir():
        raise SystemExit(f"ERROR: raw directory does not exist: {args.raw_dir}")

    if args.samples is not None:
        pairs = read_samples_tsv(args.samples, args.raw_dir)
    else:
        pairs = discover_samples(args.raw_dir)

    names = [s for s, _ in pairs]
    if len(set(names)) != len(names):
        raise SystemExit("ERROR: duplicate sample names in samples.tsv")
    for s in names:
        if not re.fullmatch(r"[A-Za-z0-9_.-]+", s):
            raise SystemExit(
                f"ERROR: sample name '{s}' contains unsafe characters; only A-Z a-z 0-9 _ . - "
                "are allowed (the name is used in file names)"
            )

    if args.write_samples is not None:
        args.write_samples.parent.mkdir(parents=True, exist_ok=True)
        with open(args.write_samples, "w", encoding="utf-8", newline="\n") as out:
            out.write("# sample\tfasta (relative to 00_raw/ or an absolute path)\n")
            out.write("# Sample names can be edited manually; rerun this script with --force after editing.\n")
            out.write("sample\tfasta\n")
            for s, p in pairs:
                try:
                    rel = p.relative_to(args.raw_dir)
                except ValueError:
                    rel = p
                out.write(f"{s}\t{rel}\n")
        print(f"[write] sample table template -> {args.write_samples}", file=sys.stderr)

    clean_dir = args.outdir / "clean_fasta"
    clean_dir.mkdir(parents=True, exist_ok=True)

    qcs = []
    for sample, path in pairs:
        if not path.exists():
            raise SystemExit(f"ERROR: FASTA for sample {sample} does not exist: {path}")
        out_faa = clean_dir / f"{sample}.faa"
        write_output = args.force or not out_faa.exists()
        if not write_output:
            print(f"[skip] {out_faa} already exists; QC is recomputed and the file is kept "
                  "(use --force to rewrite)", file=sys.stderr)
        qc = process_sample(sample, path, out_faa, args, write_output=write_output)
        qcs.append(qc)
        print(f"[ok] {sample}: {qc.n_raw} -> {qc.n_clean} sequences", file=sys.stderr)

    with open(args.outdir / "input_qc.tsv", "w", encoding="utf-8", newline="\n") as out:
        cols = ["sample", "raw_fasta", "clean_fasta", "n_raw", "n_clean", "n_short",
                "n_high_x", "n_internal_stop", "n_stop_dropped", "n_uo", "n_other",
                "n_dup", "total_aa", "min_len", "max_len", "warnings"]
        out.write("\t".join(cols) + "\n")
        for q in qcs:
            out.write("\t".join(str(x) for x in [
                q.sample, q.raw_fasta, q.clean_fasta, q.n_raw, q.n_clean, q.n_short,
                q.n_high_x, q.n_internal_stop, q.n_stop_dropped, q.n_uo, q.n_other,
                q.n_dup, q.total_aa, q.min_len, q.max_len, "; ".join(q.warnings),
            ]) + "\n")

    with open(args.outdir / "samples.clean.tsv", "w", encoding="utf-8", newline="\n") as out:
        out.write("sample\tfasta\tn_raw\tn_clean\n")
        for q in qcs:
            out.write(f"{q.sample}\t{q.clean_fasta}\t{q.n_raw}\t{q.n_clean}\n")

    write_markdown(qcs, args.outdir / "input_qc.md", args)
    print(f"[done] QC report: {args.outdir / 'input_qc.md'}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
