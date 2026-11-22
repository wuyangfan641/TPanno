#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
03_extract_domains.py
=====================
Parse HMMER domtblout files, group hits into protein families by HMM model,
filter domain hits, and extract domain subsequences for MAFFT/trimAl/IQ-TREE.

Key scientific conventions
--------------------------
* Domain boundaries: by default, use the HMMER ali coordinates (the region that
  was actually aligned), not the whole protein. Use --coord env to switch to
  envelope coordinates. The phylogenetic signal then comes from the homologous
  region modeled by the HMM, avoiding the influence of non-homologous N- and
  C-terminal segments, alternative splicing, and low-complexity regions.
* Domain completeness: require (hmm_to - hmm_from + 1) / qlen >= --hmm-cov to
  prevent fragmented hits covering only a small part of the model from
  contaminating the alignment.
* By default, keep only the best domain copy per protein per family
  (--all-domains retains tandem duplicates), so that multiple copies from the
  same protein are not counted as independent taxonomic units.
* All filtering thresholds, coordinates, and E-values are written to
  metadata.tsv, so every result is traceable.

Usage:
  python bin/03_extract_domains.py \
      --search-dir results/02_search \
      --samples results/00_prep/samples.clean.tsv \
      --outdir results/03_domains
"""
from __future__ import annotations

import argparse
import itertools
import os
import re
import sqlite3
import sys
from pathlib import Path

try:
    from Bio import SeqIO
except ImportError:  # pragma: no cover
    sys.exit("ERROR: Biopython is required.")

# domtblout columns (0-based); the HMMER3 format has exactly 23 columns
C_TARGET, C_TACC, C_TLEN = 0, 1, 2
C_QUERY, C_QACC, C_QLEN = 3, 4, 5
C_FULL_E, C_FULL_SCORE = 6, 7
C_DOM_N, C_DOM_OF = 9, 10
C_CE, C_IE, C_DOM_SCORE = 11, 12, 13
C_HMM_FROM, C_HMM_TO = 15, 16
C_ALI_FROM, C_ALI_TO = 17, 18
C_ENV_FROM, C_ENV_TO = 19, 20
C_ACC, C_DESC = 21, 22

COLUMNS = [
    "family", "sample", "protein", "tlen", "qlen", "i_evalue", "score",
    "hmm_from", "hmm_to", "ali_from", "ali_to", "env_from", "env_to",
    "hmm_cov", "description", "search_file",
]


def eprint(*a):
    print(*a, file=sys.stderr)


def safe_token(s: str, maxlen: int = 120) -> str:
    """Convert an arbitrary string into a Newick- and filename-safe token."""
    s = re.sub(r"[^A-Za-z0-9_.-]+", "_", s).strip("_")
    if not s:
        s = "x"
    if len(s) > maxlen:
        import hashlib
        s = s[:maxlen - 9] + "_" + hashlib.sha1(s.encode()).hexdigest()[:8]
    return s


def parse_domtbl(path: Path):
    with open(path, "r", encoding="utf-8-sig", errors="replace") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split(None, C_DESC)
            if len(f) <= C_DESC:
                continue
            yield f


def as_float(x: str, default: float = float("inf")) -> float:
    try:
        return float(x)
    except (TypeError, ValueError):
        return default


def as_int(x: str, default: int = 0) -> int:
    try:
        return int(x)
    except (TypeError, ValueError):
        return default


def load_samples(path: Path) -> dict[str, Path]:
    out = {}
    with open(path, "r", encoding="utf-8-sig") as fh:
        header = None
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\r\n").split("\t")
            if header is None:
                header = [x.strip().lower() for x in parts]
                for req in ("sample", "fasta"):
                    if req not in header:
                        raise SystemExit(f"ERROR: {path} is missing the column '{req}' (header: {header})")
                continue
            rec = dict(zip(header, [x.strip() for x in parts]))
            if "sample" not in rec or "fasta" not in rec:
                raise SystemExit(f"ERROR: malformed row in {path}: {line.rstrip()}")
            if rec["sample"] in out:
                raise SystemExit(f"ERROR: duplicate sample '{rec['sample']}' in {path}")
            out[rec["sample"]] = Path(rec["fasta"])
    if not out:
        raise SystemExit(f"ERROR: no sample records found in {path}")
    return out


def build_db(search_dir: Path, db_path: Path, args) -> sqlite3.Connection:
    if db_path.exists():
        db_path.unlink()
    conn = sqlite3.connect(db_path)
    conn.executescript("""
        PRAGMA journal_mode=OFF;
        PRAGMA synchronous=OFF;
        PRAGMA temp_store=MEMORY;
        CREATE TABLE hits (
            family TEXT, sample TEXT, protein TEXT,
            tlen INTEGER, qlen INTEGER,
            i_evalue REAL, score REAL,
            hmm_from INTEGER, hmm_to INTEGER,
            ali_from INTEGER, ali_to INTEGER,
            env_from INTEGER, env_to INTEGER,
            hmm_cov REAL, description TEXT, search_file TEXT
        );
    """)
    insert_sql = f"INSERT INTO hits VALUES ({','.join('?' * len(COLUMNS))})"

    allowed = None
    if args.family_list:
        allowed = set()
        with open(args.family_list, "r", encoding="utf-8-sig") as fh:
            for line in fh:
                line = line.strip()
                if line and not line.startswith("#"):
                    allowed.add(line)
                    allowed.add(line.split(".")[0])

    n_in = n_keep = 0
    rows: list[tuple] = []
    files = sorted(search_dir.glob("*.domtblout"))
    if not files:
        raise SystemExit(f"ERROR: no *.domtblout files found in {search_dir}")

    for path in files:
        for f in parse_domtbl(path):
            n_in += 1
            family = f[C_QUERY].strip()
            if allowed is not None and family not in allowed and family.split(".")[0] not in allowed:
                continue
            protein = f[C_TARGET].strip()
            sample = args.sample_of.get(path.name)
            if sample is None:
                # hmmsearch: <hmmbase>__<sample>.domtblout; hmmscan: <sample>.domtblout
                stem = path.name[: -len(".domtblout")]
                sample = None
                for s in args.sample_names:
                    if stem.endswith("__" + s):
                        sample = s
                        break
                if sample is None and "__" in stem:
                    sample = stem.rsplit("__", 1)[1]
                if sample is None:
                    sample = stem
            tlen = as_int(f[C_TLEN])
            qlen = as_int(f[C_QLEN])
            i_e = as_float(f[C_IE])
            hf, ht = as_int(f[C_HMM_FROM]), as_int(f[C_HMM_TO])
            af, at = as_int(f[C_ALI_FROM]), as_int(f[C_ALI_TO])
            ef, et = as_int(f[C_ENV_FROM]), as_int(f[C_ENV_TO])

            if tlen <= 0 or qlen <= 0 or not protein:
                continue
            if i_e > args.max_evalue:
                continue
            hmm_cov = (ht - hf + 1) / qlen if ht >= hf else 0.0
            if hmm_cov < args.hmm_cov:
                continue
            if not (1 <= af <= at <= tlen):
                continue
            if (at - af + 1) < args.min_domain_len:
                continue

            desc = f[C_DESC].strip() if len(f) > C_DESC else ""
            rows.append((family, sample, protein, tlen, qlen, i_e, as_float(f[C_DOM_SCORE], 0.0),
                         hf, ht, af, at, ef, et, hmm_cov, desc, path.name))
            n_keep += 1
        if len(rows) >= 5000:
            conn.executemany(insert_sql, rows)
            rows.clear()

    if rows:
        conn.executemany(insert_sql, rows)
    conn.commit()
    conn.execute("CREATE INDEX idx_family ON hits(family)")
    conn.commit()
    eprint(f"[db] read {n_in:,} domain hits; {n_keep:,} passed domain-level filters")
    return conn


def extract_family(conn, family: str, rows, indexes, args, outdir: Path):
    """rows are already sorted by sample, protein, and i_evalue."""
    records = []
    n_domains_total = 0
    n_missing = 0

    for (sample, protein), grp in itertools.groupby(rows, key=lambda r: (r["sample"], r["protein"])):
        grp = list(grp)
        n_all = len(grp)
        n_domains_total += n_all
        if not args.all_domains:
            # Prioritize domain completeness (HMM coverage), then break ties by
            # i-Evalue and score. Ranking by i-Evalue alone favours "local but
            # highly significant" fragments and biases the phylogenetic alignment.
            grp = [max(grp, key=lambda r: (round(r["hmm_cov"], 3), -r["i_evalue"], r["score"]))]

        idx = indexes.get(sample)
        if idx is None:
            eprint(f"WARNING: no FASTA index for sample {sample}; skipping {protein}")
            n_missing += 1
            continue
        try:
            full = str(idx[protein].seq)
        except KeyError:
            n_missing += 1
            continue

        # Number tandem duplicates by coordinate for reproducible output
        grp = sorted(grp, key=lambda r: (r["ali_from"], r["ali_to"]))
        for k, r in enumerate(grp, 1):
            if args.coord == "env":
                start, end = r["env_from"], r["env_to"]
            else:
                start, end = r["ali_from"], r["ali_to"]
            start = max(1, start - args.flank)
            end = min(len(full), end + args.flank)
            if end < start:
                continue
            sub = full[start - 1:end].upper()
            if not sub or len(sub) < args.min_domain_len:
                continue
            xfrac = sub.count("X") / len(sub)
            if xfrac > args.max_x:
                continue

            base = f"{sample}__{protein}"
            uid = safe_token(base) if not args.all_domains else safe_token(f"{base}__d{k}")
            records.append({
                "seq_id": uid,
                "family": family,
                "sample": sample,
                "protein_id": protein,
                "domain_no": k if args.all_domains else 1,
                "n_domains": n_all,
                "i_evalue": r["i_evalue"],
                "score": r["score"],
                "hmm_cov": r["hmm_cov"],
                "hmm_from": r["hmm_from"], "hmm_to": r["hmm_to"],
                "ali_from": r["ali_from"], "ali_to": r["ali_to"],
                "env_from": r["env_from"], "env_to": r["env_to"],
                "protein_len": len(full),
                "domain_len": len(sub),
                "x_frac": xfrac,
                "description": re.sub(r"\s+", " ", r["description"] or ""),
                "search_file": r["search_file"],
                "seq": sub,
            })

    # Uniqueness fallback (safe_token can theoretically map different IDs to the same token)
    seen = {}
    for rec in records:
        sid = rec["seq_id"]
        if sid in seen:
            seen[sid] += 1
            rec["seq_id"] = f"{sid}_{seen[sid]}"
        else:
            seen[sid] = 1

    return records, n_domains_total, n_missing


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--search-dir", type=Path, required=True)
    ap.add_argument("--samples", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--hmm-cov", type=float, default=float(os.environ.get("HMM_COV_MIN", 0.70)))
    ap.add_argument("--max-evalue", type=float, default=float(os.environ.get("HMM_EVALUE_MAX", 1e-3)))
    ap.add_argument("--min-domain-len", type=int, default=int(os.environ.get("DOMAIN_MIN_LEN", 30)))
    ap.add_argument("--max-x", type=float, default=float(os.environ.get("DOMAIN_MAX_X_FRAC", 0.30)))
    ap.add_argument("--coord", choices=["ali", "env"], default=os.environ.get("DOMAIN_COORD", "ali"))
    ap.add_argument("--flank", type=int, default=int(os.environ.get("DOMAIN_FLANK", 0)))
    ap.add_argument("--all-domains", action="store_true",
                    default=os.environ.get("KEEP_ALL_DOMAINS", "0") == "1")
    ap.add_argument("--min-seqs", type=int, default=int(os.environ.get("MIN_SEQS", 4)))
    ap.add_argument("--min-samples", type=int, default=int(os.environ.get("MIN_SAMPLES", 2)))
    ap.add_argument("--family-list", type=Path, default=None)
    ap.add_argument("--keep-db", action="store_true")
    args = ap.parse_args()

    if not args.samples.exists():
        raise SystemExit(f"ERROR: sample table does not exist: {args.samples}")
    samples = load_samples(args.samples)
    if args.min_samples > len(samples):
        eprint(f"WARNING: MIN_SAMPLES={args.min_samples} exceeds the number of samples "
               f"({len(samples)}); all families will be skipped.")
    args.sample_of = {}
    args.sample_names = list(samples)
    # Infer the sample from the file name (needed only in hmmscan mode); also
    # build an index for each cleaned FASTA file
    indexes = {}
    for s, fa in samples.items():
        if not fa.exists():
            raise SystemExit(f"ERROR: cleaned FASTA file does not exist: {fa}")
        indexes[s] = SeqIO.index(str(fa), "fasta")
        args.sample_of[f"{s}.domtblout"] = s

    args.outdir.mkdir(parents=True, exist_ok=True)
    db_path = args.outdir / "hits.sqlite"
    conn = build_db(args.search_dir, db_path, args)
    conn.row_factory = sqlite3.Row

    families = [r[0] for r in conn.execute(
        "SELECT family FROM hits GROUP BY family ORDER BY family")]

    summary = []
    skipped = []
    total_written = 0

    for family in families:
        rows = conn.execute(
            "SELECT * FROM hits WHERE family=? ORDER BY sample, protein, i_evalue, score DESC",
            (family,)).fetchall()
        proteins = {(r["sample"], r["protein"]) for r in rows}
        sample_set = {r["sample"] for r in rows}
        if len(proteins) < args.min_seqs or len(sample_set) < args.min_samples:
            skipped.append((family, len(rows), len(proteins), len(sample_set),
                            f"proteins<{args.min_seqs} or samples<{args.min_samples}"))
            continue

        records, n_dom, n_missing = extract_family(conn, family, rows, indexes, args, args.outdir)
        if len(records) < args.min_seqs:
            skipped.append((family, len(rows), len(records), len(sample_set),
                            "too few sequences after extraction"))
            continue

        fdir_name = safe_token(family, 80)
        fdir = args.outdir / fdir_name
        fdir.mkdir(parents=True, exist_ok=True)
        faa = fdir / f"{fdir_name}.domains.faa"
        meta = fdir / "metadata.tsv"

        with open(faa, "w", encoding="utf-8", newline="\n") as out:
            for r in records:
                out.write(f">{r['seq_id']}\n")
                s = r["seq"]
                for i in range(0, len(s), 60):
                    out.write(s[i:i + 60] + "\n")

        with open(meta, "w", encoding="utf-8", newline="\n") as out:
            cols = ["seq_id", "family", "sample", "protein_id", "domain_no", "n_domains",
                    "i_evalue", "score", "hmm_cov", "hmm_from", "hmm_to", "ali_from",
                    "ali_to", "env_from", "env_to", "protein_len", "domain_len", "x_frac",
                    "description", "search_file"]
            out.write("\t".join(cols) + "\n")
            for r in records:
                out.write("\t".join(str(r[c]) for c in cols) + "\n")

        summary.append((family, fdir_name, len(rows), n_dom, len(records),
                        len({r["sample"] for r in records}), n_missing, "ok"))
        total_written += len(records)

    conn.close()
    if not args.keep_db:
        db_path.unlink(missing_ok=True)

    with open(args.outdir / "families.tsv", "w", encoding="utf-8", newline="\n") as out:
        out.write("family\tdir\tn_hits\tn_domains\tn_seqs\tn_samples\tn_missing\tstatus\n")
        for row in summary:
            out.write("\t".join(map(str, row)) + "\n")
        for fam, nh, np_, ns, why in skipped:
            out.write(f"{fam}\tNA\t{nh}\tNA\t{np_}\t{ns}\t0\tskipped: {why}\n")

    with open(args.outdir / "families_skipped.tsv", "w", encoding="utf-8", newline="\n") as out:
        out.write("family\tn_hits\tn_proteins\tn_samples\treason\n")
        for row in skipped:
            out.write("\t".join(map(str, row)) + "\n")

    eprint(f"[done] wrote {len(summary)} families and {total_written:,} domain sequences -> {args.outdir}")
    if skipped:
        eprint(f"[skip] {len(skipped)} families did not meet the thresholds; see families_skipped.tsv")
    return 0


if __name__ == "__main__":
    sys.exit(main())
