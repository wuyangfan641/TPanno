#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
06_make_report.py
=================
Summarize the results of steps 03/04/05 and write:

  results/06_report/family_summary.tsv   one row per family (tree + cluster metrics)
  results/06_report/report.md            human-readable report with QC alerts
  results/06_report/parameters.tsv       key parameters of this run (for the Methods section)

Usage:
  python bin/06_make_report.py --results results
"""
from __future__ import annotations

import argparse
import csv
import os
import re
import sys
from collections import defaultdict
from pathlib import Path


def load_config_defaults() -> None:
    """Fill parameters.tsv when this script is run outside run_pipeline.sh.

    run_pipeline.sh sources config/config.sh and exports every parameter, so the
    environment is complete in the normal path. When 06 is invoked directly, the
    environment is empty; parse the simple KEY=value assignments in config.sh so
    that the parameter snapshot is still meaningful. Existing environment
    variables always win.
    """
    cfg = Path(__file__).resolve().parents[1] / "config" / "config.sh"
    if not cfg.exists():
        return
    assign = re.compile(r"^(?:export\s+)?([A-Z_][A-Z0-9_]*)=(.*)$")
    ref = re.compile(r"^\$\{([A-Z_][A-Z0-9_]*):-([^}]*)\}$")
    text = cfg.read_text(encoding="utf-8-sig", errors="replace")
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        m = assign.match(line)
        if not m:
            continue
        key, raw = m.group(1), m.group(2).strip()
        # Drop an unquoted trailing comment such as: KEY="value"  # explanation
        raw = re.split(r"\s+#", raw, maxsplit=1)[0].strip()
        if len(raw) >= 2 and raw[0] == raw[-1] and raw[0] in ("'", '"'):
            raw = raw[1:-1]
        r = ref.match(raw)
        if r:
            raw = os.environ.get(r.group(1), r.group(2))
        elif raw.startswith("${") or "$(" in raw:
            # Paths and command substitutions are not simple parameter values.
            continue
        if key not in os.environ:
            os.environ[key] = raw


def read_tsv(path: Path) -> list[dict]:
    if not path.exists():
        return []
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


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", type=Path, default=Path("results"))
    args = ap.parse_args()
    load_config_defaults()
    R = args.results
    outdir = R / "06_report"
    outdir.mkdir(parents=True, exist_ok=True)

    families_path = R / "03_domains" / "families.tsv"
    if not families_path.exists():
        raise SystemExit(f"ERROR: {families_path} is missing; run 03_extract_domains.py first")
    families = read_tsv(families_path)
    # 04 records the sanitized directory name as the family; map it back to the
    # original family name so that the join with 03/05 works for any HMM name.
    dir_to_family = {r["dir"]: r["family"] for r in families
                     if r.get("dir") and r["dir"] != "NA"}
    status = {}
    for r in read_tsv(R / "04_phylo" / "all_status.tsv"):
        r = dict(r)
        r["family"] = dir_to_family.get(r["family"], r["family"])
        status[r["family"]] = r
    clusters = read_tsv(R / "05_clusters" / "all_cluster_summary.tsv")

    by_family: dict[str, list[dict]] = defaultdict(list)
    for c in clusters:
        by_family[c["family"]].append(c)

    rows = []
    for f in families:
        fam = f["family"]
        st = status.get(fam)
        base_status = f.get("status", "NA")
        if st is None:
            # Passed 03 but has no 04 status: usually screened out by 03b
            # (not selected), so it must not count as a successful tree.
            run_status = "not_run" if base_status == "ok" else base_status
        else:
            run_status = st.get("status", "NA")
        cl = by_family.get(fam, [])
        n_clusters = len(cl)
        sizes = [as_int(c["n_seqs"]) for c in cl]
        ge_min = sum(1 for c in cl if c.get("passes_min_size") == "yes")
        total = sum(sizes)
        mean_sil = (sum(as_float(c["mean_silhouette"]) * as_int(c["n_seqs"]) for c in cl) / total
                    if total else 0.0)
        rows.append({
            "family": fam,
            "status": run_status,
            "n_hits": f.get("n_hits", "NA"),
            "n_domains": f.get("n_domains", "NA"),
            "n_seqs": f.get("n_seqs", "NA"),
            "n_samples": f.get("n_samples", "NA"),
            "n_used_tree": st.get("n_used", "NA") if st else "NA",
            "aln_cols": st.get("n_aln_cols", "NA") if st else "NA",
            "trim_cols": st.get("n_trim_cols", "NA") if st else "NA",
            "trim_fraction": st.get("trim_fraction", "NA") if st else "NA",
            "best_model": st.get("best_model", "NA") if st else "NA",
            "n_clusters": n_clusters,
            "n_clusters_ge_min": ge_min,
            "largest_cluster": max(sizes) if sizes else 0,
            "singleton_fraction": (sum(1 for s in sizes if s == 1) / n_clusters) if n_clusters else 0.0,
            "mean_silhouette": round(mean_sil, 4),
            "notes": (st.get("notes", "") if st else
                      ("not_run_or_screened_out" if run_status == "not_run" else base_status)),
        })

    cols = ["family", "status", "n_hits", "n_domains", "n_seqs", "n_samples",
            "n_used_tree", "aln_cols", "trim_cols", "trim_fraction", "best_model",
            "n_clusters", "n_clusters_ge_min", "largest_cluster", "singleton_fraction",
            "mean_silhouette", "notes"]
    with open(outdir / "family_summary.tsv", "w", encoding="utf-8", newline="\n") as out:
        out.write("\t".join(cols) + "\n")
        for r in rows:
            out.write("\t".join(str(r[c]) for c in cols) + "\n")

    # ---------- Parameter snapshot --------------------------------------------
    param_keys = [
        "SEARCH_MODE", "HMM_CUT_MODE", "HMM_EVALUE_MAX", "HMM_COV_MIN", "DOMAIN_COORD",
        "DOMAIN_FLANK", "DOMAIN_MIN_LEN", "DOMAIN_MAX_X_FRAC", "KEEP_ALL_DOMAINS",
        "MIN_PROTEIN_LEN", "MAX_X_FRAC", "CLEAN_UO_TO_X", "DROP_INTERNAL_STOP",
        "MIN_SEQS", "MIN_SAMPLES", "MAX_FAMILIES", "MIN_MEDIAN_COV",
        "MIN_MEDIAN_DOMAIN_LEN", "MAX_MEDIAN_COPIES", "MIN_SINGLE_COPY_FRAC",
        "CLAN_MAP", "COLLAPSE_CLANS", "DEREPLICATE", "DEREP_ID", "DEREP_MAX_SEQS",
        "ALIGNER", "MAFFT_MODE", "MAFFT_MAXITERATE", "MAFFT_EXTRA", "HMM_SOURCE",
        "TRIMAL_MODE", "TRIMAL_MIN_COL", "MIN_COL_FRACTION",
        "IQTREE_MODEL", "IQTREE_MSET", "UFBOOT", "ALRT", "IQTREE_SEED",
        "OUTGROUP", "IQTREE_EXTRA", "CLUSTER_CUTOFF", "CLUSTER_LINKAGE",
        "MIN_CLUSTER_SIZE", "CLUSTER_SCAN",
    ]
    with open(outdir / "parameters.tsv", "w", encoding="utf-8", newline="\n") as out:
        out.write("parameter\tvalue\n")
        for k in param_keys:
            out.write(f"{k}\t{os.environ.get(k, 'NA')}\n")

    # ---------- Markdown report ------------------------------------------------
    n_total = len(rows)
    n_ok = sum(1 for r in rows if r["status"] == "ok")
    n_failed = sum(1 for r in rows if r["status"] == "failed")
    n_skipped = sum(1 for r in rows if r["status"].startswith("skipped"))
    n_not_run = sum(1 for r in rows if r["status"] == "not_run")
    n_cl = sum(r["n_clusters"] for r in rows)
    n_seq_dom = sum(as_int(r["n_seqs"]) for r in rows
                    if not str(r["status"]).startswith("skipped")
                    and r["n_seqs"] not in ("NA", ""))
    n_seq_tree = sum(as_int(r["n_used_tree"]) for r in rows
                     if r["status"] == "ok" and r["n_used_tree"] not in ("NA", ""))
    max_rows = int(os.environ.get("MAX_REPORT_ROWS", "200"))

    # QC alerts (a full Pfam scan can produce thousands; list at most 100)
    alerts = []
    for r in rows:
        if r["status"] == "ok":
            if as_float(r["trim_fraction"], 1.0) < 0.5:
                alerts.append(f"{r['family']}: only {r['trim_fraction']} of the alignment columns "
                              "remain after trimAl; tree uncertainty may be high, inspect the alignment.")
            if as_int(r["n_used_tree"]) < 10:
                alerts.append(f"{r['family']}: only {r['n_used_tree']} sequences entered tree "
                              "inference; topology and support may be unstable.")
            if r["best_model"] in ("NA", ""):
                alerts.append(f"{r['family']}: could not parse the best-fit substitution model "
                              "(a fixed model may have been used).")
        elif r["status"] == "failed":
            alerts.append(f"{r['family']}: pipeline failure -- {r['notes']}")

    lines = [
        "# Domain-based phylogeny and domain-family clustering -- results report",
        "",
        "## 1. Overview",
        "",
        f"- HMM/domain families analyzed: **{n_total}**",
        f"- Successful trees: **{n_ok}**; skipped: **{n_skipped}**; failed: **{n_failed}**; "
        f"not run (screened out or not yet inferred): **{n_not_run}**",
        f"- Domain sequences in the 03 output (including unselected families): **{n_seq_dom:,}**",
        f"- Sequences actually used in domain trees: **{n_seq_tree:,}**",
        f"- Total clusters: **{n_cl}**",
        "",
        "## 2. Family-level results",
        "",
        f"> The table below shows at most {max_rows} rows (sorted by domain sequence count); "
        "the complete results are in `06_report/family_summary.tsv`.",
        "",
        "| Family | Status | Domain sequences | Tree sequences | Alignment columns | After trimming | Best model | Clusters | Largest cluster | Mean silhouette |",
        "|---|---|---:|---:|---:|---:|---|---:|---:|---:|",
    ]
    rows_sorted = sorted(rows, key=lambda x: (-as_int(x["n_seqs"]), x["family"]))
    for r in rows_sorted[:max_rows]:
        lines.append(
            f"| {r['family']} | {r['status']} | {r['n_seqs']} | {r['n_used_tree']} | "
            f"{r['aln_cols']} | {r['trim_cols']} | {r['best_model']} | {r['n_clusters']} | "
            f"{r['largest_cluster']} | {r['mean_silhouette']} |"
        )
    if n_total > max_rows:
        lines.append(f"\n> The remaining {n_total - max_rows} families are not expanded here.")

    lines += ["", "## 3. The 15 largest clusters", "",
              "| Family | Cluster | Members | Samples | Mean within | Min between | Silhouette | Monophyletic | UFBoot |",
              "|---|---|---:|---:|---:|---:|---:|---|---:|"]
    top = sorted(clusters, key=lambda c: -as_int(c.get("n_members", c["n_seqs"])))[:15]
    for c in top:
        lines.append(
            f"| {c['family']} | {c['cluster']} | {c.get('n_members', c['n_seqs'])} | {c['n_samples']} | "
            f"{c['mean_within']} | {c['min_between']} | {c['mean_silhouette']} | "
            f"{c['monophyletic']} | {c['UFBoot']} |"
        )

    lines += ["", "## 4. QC alerts", ""]
    if alerts:
        lines += [f"- {a}" for a in alerts[:100]]
        if len(alerts) > 100:
            lines.append(f"- The remaining {len(alerts) - 100} alerts are not expanded; see each "
                         "family's `status.tsv`.")
    else:
        lines.append("- No preset QC alerts were triggered.")
    lines += [
        "",
        "## 5. Notes for interpretation",
        "",
        "- **Monophyly and support:** `UFBoot` is ultrafast bootstrap support and `SH-aLRT` is the "
        "Shimodaira-Hasegawa approximate likelihood ratio test. Only branches with high values of "
        "both (for example, UFBoot >= 95 and SH-aLRT >= 80) should be treated as strong topological "
        "evidence; low-support branches should not be forced into subfamilies.",
        "- **Thresholds are not objective boundaries:** the patristic distance threshold t is a "
        "descriptive choice. Interpret cluster biology together with the silhouette plateau across "
        "thresholds in `*.cutoff_scan.tsv` and with independent evidence such as domain architecture, "
        "gene synteny, expression, or functional experiments.",
        "- **Domain tree is not a gene tree or a species tree:** when only one tandem copy per protein "
        "is retained, the tree reflects the evolution of domain copies. If all copies are retained, "
        "note that multiple copies from the same protein are not independent. Cross-species comparisons "
        "must also distinguish orthologs from paralogs.",
        "- **A full Pfam scan is not unsupervised clustering:** it relies on existing Pfam models, so "
        "families absent from Pfam do not appear. Several HMMs from the same clan may describe the same "
        "domain; use `03b_screen_families.py --collapse-clans` to remove redundancy, or merge by clan "
        "during downstream interpretation.",
        "- **Long branches and saturation:** very long branches at deep nodes may reflect saturation or "
        "alignment error. Inspect `iqtree.treefile` and the alignment, and if necessary run sensitivity "
        "analyses with long-branch sequences removed or with different trimming parameters.",
        "",
        "## 6. File index",
        "",
        "| Path | Content |",
        "|---|---|",
        "| `00_prep/input_qc.md` | Input sequence QC |",
        "| `03_domains/families.tsv` | Initial statistics for all domain families |",
        "| `03_domains/screen/family_screen.tsv` | Full-Pfam screening and copy-number metrics |",
        "| `03_domains/screen/family_list.txt` | Selected families |",
        "| `03_domains/<family>/metadata.tsv` | Coordinates and E-values of each domain sequence |",
        "| `04_phylo/<family>/iqtree.treefile` | ML tree (with UFBoot/SH-aLRT) |",
        "| `04_phylo/<family>/iqtree.iqtree` | ModelFinder and tree statistics |",
        "| `05_clusters/<family>.clusters.tsv` | Sequence-level cluster assignments |",
        "| `05_clusters/<family>.cluster_summary.tsv` | Cluster-level statistics |",
        "| `05_clusters/<family>.cutoff_scan.tsv` | Threshold sensitivity scan |",
        "| `06_report/family_summary.tsv` | Family-level summary (this table) |",
    ]
    (outdir / "report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[done] report -> {outdir / 'report.md'}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
