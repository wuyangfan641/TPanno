#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
03b_screen_families.py
======================
Family screening and prioritization after a full Pfam-A `hmmscan` run.

Why this step is needed
-----------------------
A full Pfam scan yields thousands of domain families. Sending every family to
MAFFT + IQ-TREE is neither feasible nor useful: many families are paralogous,
fragmented, low-coverage, or redundant models from the same clan. This script
runs after domain extraction and before tree inference and screens families
with interpretable metrics:

  * sample coverage: n_samples >= MIN_SAMPLES
  * sequence count: n_seqs >= MIN_SEQS
  * domain quality: median HMM coverage, median domain length
  * copy number: median copies and single-copy sample fraction (paralogy diagnostics)
  * clan redundancy: keep only the highest-scoring model per Pfam clan (optional)
  * computational cap: MAX_FAMILIES (0 = no limit)

Outputs
-------
  <outdir>/family_screen.tsv     metrics for every family and selection status
  <outdir>/family_list.txt       selected family names (one per line)
  <outdir>/family_dirs.txt       selected family directory names (used by 04/05)
  <outdir>/selected_families.tsv selected families only
  <outdir>/screen_report.md      human-readable report

Usage
-----
  python3 bin/03b_screen_families.py \
      --domain-dir results/03_domains \
      --outdir results/03_domains/screen \
      --min-samples 22 --min-seqs 22 --max-families 300 \
      --clan-map 01_hmms/Pfam-A.clans.tsv --collapse-clans
"""
from __future__ import annotations

import argparse
import csv
import os
import statistics
import sys
from collections import Counter, defaultdict
from pathlib import Path


def eprint(*a):
    print(*a, file=sys.stderr)


def read_tsv(path: Path) -> list[dict]:
    with open(path, "r", encoding="utf-8-sig") as fh:
        return list(csv.DictReader(fh, delimiter="\t", quoting=csv.QUOTE_NONE))


def as_float(x, default=0.0):
    try:
        return float(x)
    except (TypeError, ValueError):
        return default


def as_int(x, default=0):
    try:
        return int(float(x))
    except (TypeError, ValueError):
        return default


def load_clan_map(path: Path) -> dict[str, str]:
    """Pfam-A.clans.tsv: family_accession<TAB>clan_id<TAB>clan_name..."""
    out: dict[str, str] = {}
    with open(path, "r", encoding="utf-8-sig", errors="replace") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\r\n").split("\t")
            if len(f) < 2:
                continue
            acc, clan = f[0].strip(), f[1].strip()
            if acc.lower().startswith("family"):
                continue
            if acc and clan:
                out[acc.split(".")[0]] = clan
    return out


def metrics_for_family(domain_dir: Path, row: dict) -> dict:
    family = row["family"]
    dirname = row.get("dir", "NA")
    rec = {
        "family": family,
        "dir": dirname,
        "clan": "",
        "n_hits": as_int(row.get("n_hits")),
        "n_seqs": 0,
        "n_samples": 0,
        "single_copy_samples": 0,
        "multi_copy_samples": 0,
        "single_copy_fraction": 0.0,
        "median_copies": 0.0,
        "max_copies": 0,
        "median_hmm_cov": 0.0,
        "median_domain_len": 0.0,
        "score": 0.0,
        "pass": False,
        "selected": "no",
        "reasons": "",
    }
    if not str(row.get("status", "")).startswith("ok") or dirname in ("", "NA"):
        rec["reasons"] = "03_skipped"
        return rec

    meta_path = domain_dir / dirname / "metadata.tsv"
    if not meta_path.exists():
        rec["reasons"] = "metadata_missing"
        return rec

    rows = read_tsv(meta_path)
    if not rows:
        rec["reasons"] = "metadata_empty"
        return rec

    counts = Counter(r["sample"] for r in rows)
    covs = [as_float(r.get("hmm_cov")) for r in rows]
    lens = [as_int(r.get("domain_len")) for r in rows]
    n_samples = len(counts)
    n_single = sum(1 for v in counts.values() if v == 1)
    n_seqs = len(rows)

    rec.update({
        "n_seqs": n_seqs,
        "n_samples": n_samples,
        "single_copy_samples": n_single,
        "multi_copy_samples": n_samples - n_single,
        "single_copy_fraction": round(n_single / n_samples, 4) if n_samples else 0.0,
        "median_copies": round(float(statistics.median(counts.values())), 3) if counts else 0.0,
        "max_copies": max(counts.values()) if counts else 0,
        "median_hmm_cov": round(float(statistics.median(covs)), 4) if covs else 0.0,
        "median_domain_len": round(float(statistics.median(lens)), 1) if lens else 0.0,
    })
    # Ranking score: cross-sample conservation first, then single-copy fraction,
    # with a small weight on sequence count.
    rec["score"] = round(n_samples + 10.0 * rec["single_copy_fraction"]
                         + min(n_seqs, 200) / 200.0, 4)
    return rec


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--domain-dir", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--min-samples", type=int, default=int(os.environ.get("MIN_SAMPLES", 22)))
    ap.add_argument("--min-seqs", type=int, default=int(os.environ.get("MIN_SEQS", 22)))
    ap.add_argument("--min-median-cov", type=float,
                    default=float(os.environ.get("MIN_MEDIAN_COV", 0.75)))
    ap.add_argument("--min-median-domain-len", type=float,
                    default=float(os.environ.get("MIN_MEDIAN_DOMAIN_LEN", 30)))
    ap.add_argument("--max-median-copies", type=float,
                    default=float(os.environ.get("MAX_MEDIAN_COPIES", 0)),
                    help="0 = no limit; 1 = only families with median copies <= 1")
    ap.add_argument("--min-single-copy-frac", type=float,
                    default=float(os.environ.get("MIN_SINGLE_COPY_FRAC", 0.0)))
    ap.add_argument("--max-families", type=int, default=int(os.environ.get("MAX_FAMILIES", 0)),
                    help="0 = no limit; 200-500 is recommended for the first round")
    ap.add_argument("--clan-map", type=Path, default=None)
    ap.add_argument("--collapse-clans", action=argparse.BooleanOptionalAction,
                    default=os.environ.get("COLLAPSE_CLANS", "0") == "1")
    args = ap.parse_args()

    fam_tsv = args.domain_dir / "families.tsv"
    if not fam_tsv.exists():
        raise SystemExit(f"ERROR: {fam_tsv} is missing; run 03_extract_domains.py first")
    args.outdir.mkdir(parents=True, exist_ok=True)

    clan_map: dict[str, str] = {}
    if args.clan_map:
        if not args.clan_map.exists():
            raise SystemExit(f"ERROR: clan map does not exist: {args.clan_map}")
        clan_map = load_clan_map(args.clan_map)
        eprint(f"[screen] read {len(clan_map):,} Pfam clan assignments")
    if args.collapse_clans and not clan_map:
        eprint("[screen] WARNING: --collapse-clans was requested without --clan-map; "
               "clan dereplication is skipped")

    rows = []
    for raw in read_tsv(fam_tsv):
        rec = metrics_for_family(args.domain_dir, raw)
        rec["clan"] = clan_map.get(rec["family"].split(".")[0], "")
        rows.append(rec)

    # ---------- Metric filters -------------------------------------------------
    for r in rows:
        if r["reasons"]:
            continue
        reasons = []
        if r["n_samples"] < args.min_samples:
            reasons.append(f"samples<{args.min_samples}")
        if r["n_seqs"] < args.min_seqs:
            reasons.append(f"seqs<{args.min_seqs}")
        if r["median_hmm_cov"] < args.min_median_cov:
            reasons.append(f"cov<{args.min_median_cov}")
        if r["median_domain_len"] < args.min_median_domain_len:
            reasons.append(f"len<{args.min_median_domain_len}")
        if args.max_median_copies > 0 and r["median_copies"] > args.max_median_copies:
            reasons.append(f"median_copies>{args.max_median_copies}")
        if r["single_copy_fraction"] < args.min_single_copy_frac:
            reasons.append(f"single_copy_frac<{args.min_single_copy_frac}")
        r["reasons"] = ";".join(reasons)
        r["pass"] = not reasons

    # ---------- Clan dereplication --------------------------------------------
    n_metric_pass = sum(1 for r in rows if r["pass"])
    if args.collapse_clans and clan_map:
        by_clan: dict[str, list[dict]] = defaultdict(list)
        for r in rows:
            if r["pass"] and r["clan"]:
                by_clan[r["clan"]].append(r)
        for clan, members in by_clan.items():
            members.sort(key=lambda x: (-x["score"], -x["n_samples"], -x["n_seqs"], x["family"]))
            for r in members[1:]:
                r["pass"] = False
                r["reasons"] = (r["reasons"] + ";" if r["reasons"] else "") + f"clan_redundant:{clan}"

    # ---------- Ranking and truncation ----------------------------------------
    candidates = [r for r in rows if r["pass"]]
    candidates.sort(key=lambda x: (-x["score"], -x["n_samples"], -x["n_seqs"], x["family"]))
    if args.max_families > 0:
        selected = candidates[:args.max_families]
        for r in candidates[args.max_families:]:
            r["pass"] = False
            r["reasons"] = (r["reasons"] + ";" if r["reasons"] else "") + f"rank>{args.max_families}"
    else:
        selected = candidates
    selected_set = {r["family"] for r in selected}
    for r in rows:
        r["selected"] = "yes" if r["family"] in selected_set else "no"

    # ---------- Outputs --------------------------------------------------------
    cols = ["family", "dir", "clan", "n_hits", "n_seqs", "n_samples",
            "single_copy_samples", "multi_copy_samples", "single_copy_fraction",
            "median_copies", "max_copies", "median_hmm_cov", "median_domain_len",
            "score", "selected", "reasons"]
    with open(args.outdir / "family_screen.tsv", "w", encoding="utf-8", newline="\n") as out:
        out.write("\t".join(cols) + "\n")
        for r in sorted(rows, key=lambda x: (-x["score"], x["family"])):
            out.write("\t".join(str(r[c]) for c in cols) + "\n")

    with open(args.outdir / "selected_families.tsv", "w", encoding="utf-8", newline="\n") as out:
        out.write("\t".join(cols) + "\n")
        for r in selected:
            out.write("\t".join(str(r[c]) for c in cols) + "\n")

    with open(args.outdir / "family_list.txt", "w", encoding="utf-8", newline="\n") as out:
        for r in selected:
            out.write(r["family"] + "\n")
    with open(args.outdir / "family_dirs.txt", "w", encoding="utf-8", newline="\n") as out:
        for r in selected:
            out.write(r["dir"] + "\n")

    # ---------- Report ---------------------------------------------------------
    n_total = len(rows)
    lines = [
        "# Full Pfam scan: family screening report",
        "",
        f"- Families in the 03 output: **{n_total}**",
        f"- Passing all metrics: **{n_metric_pass}**",
        f"- Final selection (after clan dereplication and MAX_FAMILIES): **{len(selected)}**",
        "",
        "## Screening thresholds",
        "",
        f"- Samples >= {args.min_samples}",
        f"- Sequences >= {args.min_seqs}",
        f"- Median HMM coverage >= {args.min_median_cov}",
        f"- Median domain length >= {args.min_median_domain_len} aa",
        f"- Maximum median copies: {args.max_median_copies if args.max_median_copies > 0 else 'no limit'}",
        f"- Minimum single-copy sample fraction: {args.min_single_copy_frac}",
        f"- Clan dereplication: {'yes' if args.collapse_clans and clan_map else 'no'}",
        f"- MAX_FAMILIES: {args.max_families if args.max_families > 0 else 'no limit'}",
        "",
        "## Selected families (ranked by score; at most 50 shown)",
        "",
        "| family | clan | samples | sequences | single-copy fraction | median coverage | median length | score |",
        "|---|---|---:|---:|---:|---:|---:|---:|",
    ]
    for r in selected[:50]:
        lines.append(f"| {r['family']} | {r['clan'] or '-'} | {r['n_samples']} | {r['n_seqs']} | "
                     f"{r['single_copy_fraction']} | {r['median_hmm_cov']} | "
                     f"{r['median_domain_len']} | {r['score']} |")
    lines += [
        "",
        "> Note: `n_samples` is the number of samples in which the family occurs, and "
        "`single_copy_fraction` is the fraction of samples with exactly one copy. A low "
        "value indicates strong gene-duplication or paralogy effects; such trees cannot be "
        "used directly as orthogroups. `clan_redundant` means that a higher-scoring model "
        "from the same Pfam clan was already selected. If `KEEP_ALL_DOMAINS=0` was used in "
        "step 03, tandem duplicates within a protein were collapsed, so the copy numbers "
        "here are protein-level rather than domain-level.",
    ]
    (args.outdir / "screen_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    eprint(f"[screen] {n_total} families -> pass metrics {n_metric_pass} -> selected {len(selected)}")
    eprint(f"[screen] selection list: {args.outdir / 'family_list.txt'}")
    if not selected:
        eprint("ERROR: no family passed the screening. Lower MIN_SAMPLES/MIN_SEQS or check the hmmscan output.")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
