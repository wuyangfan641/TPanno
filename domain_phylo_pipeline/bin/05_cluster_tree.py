#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
05_cluster_tree.py
==================
Tree-based clustering on family phylogenies:

  1. Compute pairwise patristic distances between all tips of a rooted or
     unrooted tree.
  2. Perform hierarchical clustering on the distance matrix (default complete
     linkage: the maximum pairwise patristic distance within a cluster is <= t,
     close to TreeCluster's "cluster diameter" definition; average linkage
     tends to chain).
  3. Cut the dendrogram at the user-supplied threshold t to obtain clusters.
  4. Report each cluster's sample composition, within- and between-cluster
     distances, silhouette, monophyly, and SH-aLRT / UFBoot support.
  5. If CD-HIT dereplication was performed upstream, map the representative's
     cluster number back to all members.

Why patristic distance instead of raw sequence identity:
  * The tree already models multiple substitutions, site-rate heterogeneity,
    and the amino-acid replacement matrix.
  * Patristic distance is the expected number of substitutions per site; it is
    closer to true evolutionary divergence than raw identity and is consistent
    with the tree visualization.
  * The threshold remains a descriptive choice, so the script also writes a
    cutoff scan (silhouette vs t) to support a data-driven threshold choice
    rather than a conventional one.

Outputs:
  <prefix>.clusters.tsv          cluster assignment for every sequence
  <prefix>.cluster_summary.tsv   statistics for every cluster
  <prefix>.cutoff_scan.tsv       number of clusters / silhouette at each threshold
  <prefix>.cluster_report.md     short human-readable report

Usage:
  python bin/05_cluster_tree.py --tree results/04_phylo/PF00069.28/iqtree.treefile \
      --metadata results/03_domains/PF00069.28/metadata.tsv \
      --out-prefix results/05_clusters/PF00069.28
"""
from __future__ import annotations

import argparse
import csv
import math
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

try:
    from Bio import Phylo
    from scipy.cluster.hierarchy import fcluster, linkage
    from scipy.spatial.distance import squareform
except ImportError:  # pragma: no cover
    sys.exit("ERROR: Biopython, NumPy, and SciPy are required.")


def eprint(*a):
    print(*a, file=sys.stderr)


# ---------------------------------------------------------------------------
# Tree operations
# ---------------------------------------------------------------------------
def check_tree(tree) -> None:
    missing = [c.name for c in tree.get_terminals() if c.branch_length is None]
    if missing:
        raise SystemExit(
            f"ERROR: {len(missing)} tips lack branch lengths, so patristic distances "
            f"cannot be computed; for example {missing[:5]}"
        )


def patristic_matrix(tree) -> tuple[list[str], np.ndarray]:
    """Compute patristic distances between all tips in O(n^2)."""
    tips = tree.get_terminals()
    names = [t.name for t in tips]
    nodes = list(tree.find_clades())
    nid = {id(n): i for i, n in enumerate(nodes)}
    tip_index = {nid[id(t)]: i for i, t in enumerate(tips)}

    adj: list[list[tuple[int, float]]] = [[] for _ in nodes]
    for node in nodes:
        u = nid[id(node)]
        for child in node.clades:
            v = nid[id(child)]
            w = max(float(child.branch_length or 0.0), 0.0)
            adj[u].append((v, w))
            adj[v].append((u, w))

    n = len(tips)
    D = np.zeros((n, n), dtype=np.float64)
    for i, tip in enumerate(tips):
        start = nid[id(tip)]
        stack = [(start, 0.0, -1)]
        while stack:
            u, d, parent = stack.pop()
            ti = tip_index.get(u)
            if ti is not None:
                D[i, ti] = d
            for v, w in adj[u]:
                if v != parent:
                    stack.append((v, d + w, u))
    # Symmetrize to remove floating-point error
    D = (D + D.T) / 2.0
    np.fill_diagonal(D, 0.0)
    return names, D


def node_support(clade) -> tuple[str, str]:
    """Parse IQ-TREE's 'SH-aLRT/UFBoot' node label."""
    raw = ""
    if clade.name not in (None, ""):
        raw = str(clade.name)
    elif clade.confidence is not None:
        raw = str(clade.confidence)
    if not raw:
        return "NA", "NA"
    parts = [p.strip() for p in raw.replace(",", "/").split("/") if p.strip()]
    sh = parts[0] if parts else "NA"
    uf = parts[1] if len(parts) > 1 else "NA"
    return sh, uf


def clade_tipsets(tree, all_names: list[str]) -> dict[frozenset, tuple[str, str]]:
    """Map every internal-node split to (SH-aLRT, UFBoot).

    For an unrooted tree, a cluster is monophyletic if and only if some edge
    separates it from all other tips. Both the tip set of each clade and its
    complement are recorded, so an arbitrary rooting cannot cause a false
    negative.
    """
    out: dict[frozenset, tuple[str, str]] = {}
    all_tips = frozenset(all_names)
    for clade in tree.get_nonterminals(order="postorder"):
        tips = frozenset(t.name for t in clade.get_terminals())
        if len(tips) < 2:
            continue
        support = node_support(clade)
        out[tips] = support
        comp = all_tips - tips
        if len(comp) >= 2:
            # The same unrooted edge may carry support written from either side;
            # keep the first clade label that is seen
            out.setdefault(comp, support)
    return out


def silhouette(D: np.ndarray, labels: np.ndarray) -> np.ndarray:
    n = len(labels)
    uniq = np.unique(labels)
    s = np.zeros(n, dtype=np.float64)
    if len(uniq) < 2:
        return s
    for i in range(n):
        same = labels == labels[i]
        n_same = int(same.sum())
        if n_same <= 1:
            s[i] = 0.0
            continue
        a = D[i, same].sum() / (n_same - 1)
        b = math.inf
        for c in uniq:
            if c == labels[i]:
                continue
            other = labels == c
            if other.any():
                b = min(b, D[i, other].mean())
        s[i] = 0.0 if max(a, b) == 0 else (b - a) / max(a, b)
    return s


def cluster_once(D: np.ndarray, cutoff: float, method: str) -> np.ndarray:
    if len(D) == 1:
        return np.ones(1, dtype=int)
    Z = linkage(squareform(D, checks=False), method=method)
    return fcluster(Z, t=cutoff, criterion="distance")


def relabel(labels: np.ndarray, names: list[str]) -> tuple[np.ndarray, list[int]]:
    """Assign stable C01... labels by descending cluster size, then by the smallest tip name."""
    groups: dict[int, list[int]] = defaultdict(list)
    for i, lab in enumerate(labels):
        groups[int(lab)].append(i)
    order = sorted(groups.values(), key=lambda idxs: (-len(idxs), min(names[i] for i in idxs)))
    new = np.zeros(len(labels), dtype=int)
    for k, idxs in enumerate(order, 1):
        for i in idxs:
            new[i] = k
    return new, [len(x) for x in order]


def medoid(D: np.ndarray, idxs: list[int]) -> int:
    if len(idxs) == 1:
        return idxs[0]
    sub = D[np.ix_(idxs, idxs)]
    return idxs[int(np.argmin(sub.sum(axis=1)))]


def cluster_stats(D, labels, names, clades, samples, n_clusters, sil, member_counts,
                  tip_members, meta):
    rows = []
    for c in range(1, n_clusters + 1):
        idxs = [i for i, lab in enumerate(labels) if lab == c]
        others = [i for i, lab in enumerate(labels) if lab != c]
        within = D[np.ix_(idxs, idxs)]
        if len(idxs) > 1:
            iu = np.triu_indices(len(idxs), k=1)
            mean_within = float(within[iu].mean())
            max_within = float(within[iu].max())
        else:
            mean_within = max_within = 0.0
        if others:
            min_between = float(D[np.ix_(idxs, others)].min())
        else:
            min_between = float("nan")

        tipset = frozenset(names[i] for i in idxs)
        sh, uf = clades.get(tipset, ("NA", "NA"))
        mono = len(idxs) == 1 or tipset in clades
        rep_i = medoid(D, idxs)
        member_ids = []
        for i in idxs:
            member_ids.extend(tip_members.get(names[i], [names[i]]))
        sample_list = sorted({meta[m]["sample"] for m in member_ids if m in meta})
        counts = defaultdict(int)
        for m in member_ids:
            if m in meta:
                counts[meta[m]["sample"]] += 1
        if not sample_list:  # should not happen; defensive fallback
            sample_list = sorted({samples.get(names[i], "NA") for i in idxs})
            counts = defaultdict(int)
            for i in idxs:
                counts[samples.get(names[i], "NA")] += 1

        rows.append({
            "cluster": f"C{c:02d}",
            "n_seqs": len(idxs),
            "n_members": member_counts.get(c, len(idxs)),
            "n_samples": len(sample_list),
            "samples": ",".join(sample_list),
            "sample_counts": ",".join(f"{k}:{v}" for k, v in sorted(counts.items())),
            "mean_within": round(mean_within, 5),
            "max_within": round(max_within, 5),
            "min_between": "NA" if math.isnan(min_between) else round(min_between, 5),
            "mean_silhouette": round(float(np.mean([sil[i] for i in idxs])), 4),
            "monophyletic": "yes" if mono else "no",
            "SH_aLRT": sh,
            "UFBoot": uf,
            "representative": names[rep_i],
            "representative_sample": samples.get(names[rep_i], "NA"),
        })
    return rows


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tree", type=Path, required=True)
    ap.add_argument("--metadata", type=Path, required=True)
    ap.add_argument("--out-prefix", type=Path, required=True)
    ap.add_argument("--cutoff", type=float, default=float(__import__("os").environ.get("CLUSTER_CUTOFF", 0.5)))
    ap.add_argument("--linkage", choices=["average", "complete", "single"],
                    default=__import__("os").environ.get("CLUSTER_LINKAGE", "complete"))
    ap.add_argument("--min-size", type=int,
                    default=int(__import__("os").environ.get("MIN_CLUSTER_SIZE", 2)))
    ap.add_argument("--scan", default=__import__("os").environ.get("CLUSTER_SCAN", ""),
                    help="comma-separated threshold list, e.g. 0.1,0.2,0.5")
    ap.add_argument("--derep-map", type=Path, default=None)
    args = ap.parse_args()

    if not math.isfinite(args.cutoff) or args.cutoff < 0:
        raise SystemExit("ERROR: --cutoff must be a finite, non-negative number")

    out_prefix = args.out_prefix
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    meta = {}
    with open(args.metadata, "r", encoding="utf-8-sig") as fh:
        for row in csv.DictReader(fh, delimiter="\t", quoting=csv.QUOTE_NONE):
            meta[row["seq_id"]] = row
    family = next(iter(meta.values()))["family"] if meta else out_prefix.stem

    tree = Phylo.read(str(args.tree), "newick")
    check_tree(tree)
    names, D = patristic_matrix(tree)
    if len(set(names)) != len(names):
        raise SystemExit("ERROR: duplicate tip names in the tree")
    eprint(f"[{family}] tree tips={len(names)}, patristic matrix computed")

    missing_meta = [n for n in names if n not in meta]
    if missing_meta:
        raise SystemExit(
            f"ERROR: {len(missing_meta)} tree tips are missing from metadata.tsv, "
            f"for example {missing_meta[:5]}"
        )
    samples = {n: meta[n]["sample"] for n in names}
    clades = clade_tipsets(tree, names)

    labels = cluster_once(D, args.cutoff, args.linkage)
    labels, sizes = relabel(labels, names)
    n_clusters = len(sizes)
    sil = silhouette(D, labels)

    # ---------- Map dereplicated members back ---------------------------------
    tip_members: dict[str, list[str]] = {n: [n] for n in names}
    if args.derep_map and args.derep_map.exists():
        tmp: dict[str, list[str]] = defaultdict(list)
        with open(args.derep_map, "r", encoding="utf-8-sig") as fh:
            for row in csv.DictReader(fh, delimiter="\t", quoting=csv.QUOTE_NONE):
                tmp[row["representative"]].append(row["member"])
        for rep, members in tmp.items():
            if rep in tip_members:
                tip_members[rep] = members

    tip_idx = {n: i for i, n in enumerate(names)}
    label_of = {names[i]: int(labels[i]) for i in range(len(names))}
    member_counts = defaultdict(int)
    for tip in names:
        member_counts[label_of[tip]] += len(tip_members.get(tip, [tip]))

    # ---------- Cutoff scan ----------------------------------------------------
    if args.scan:
        cutoffs = [float(x) for x in args.scan.split(",") if x.strip()]
        with open(f"{out_prefix}.cutoff_scan.tsv", "w", encoding="utf-8", newline="\n") as out:
            out.write("cutoff\tn_clusters\tn_singletons\tlargest_cluster\tmean_silhouette\n")
            for t in cutoffs:
                lab = cluster_once(D, t, args.linkage)
                lab, sz = relabel(lab, names)
                sil_t = silhouette(D, lab)
                out.write(f"{t}\t{len(sz)}\t{sum(1 for s in sz if s == 1)}\t{max(sz)}\t"
                          f"{np.mean(sil_t):.4f}\n")
        eprint(f"[{family}] cutoff scan -> {out_prefix}.cutoff_scan.tsv")

    # ---------- Cluster statistics ---------------------------------------------
    stats = cluster_stats(D, labels, names, clades, samples, n_clusters, sil, member_counts,
                          tip_members, meta)
    with open(f"{out_prefix}.cluster_summary.tsv", "w", encoding="utf-8", newline="\n") as out:
        cols = ["family", "cluster", "n_seqs", "n_members", "n_samples", "samples", "sample_counts",
                "mean_within", "max_within", "min_between", "mean_silhouette",
                "monophyletic", "SH_aLRT", "UFBoot", "representative",
                "representative_sample", "passes_min_size"]
        out.write("\t".join(cols) + "\n")
        for r in stats:
            r["family"] = family
            r["passes_min_size"] = "yes" if r["n_members"] >= args.min_size else "no"
            out.write("\t".join(str(r[c]) for c in cols) + "\n")

    # ---------- Sequence assignments -------------------------------------------
    rep_of = {}
    for r in stats:
        rep_of[int(r["cluster"][1:])] = r["representative"]

    n_rows = 0
    missing_members = 0
    with open(f"{out_prefix}.clusters.tsv", "w", encoding="utf-8", newline="\n") as out:
        cols = ["family", "cluster", "cluster_size", "seq_id", "sample", "protein_id",
                "is_representative", "tree_tip", "silhouette", "domain_no", "i_evalue"]
        out.write("\t".join(cols) + "\n")
        for tip in names:
            c = label_of[tip]
            for member in tip_members.get(tip, [tip]):
                m = meta.get(member)
                if m is None:
                    missing_members += 1
                    continue
                out.write("\t".join(str(x) for x in [
                    family, f"C{c:02d}", member_counts[c], member, m["sample"], m["protein_id"],
                    "yes" if member == rep_of.get(c) else "no", tip,
                    f"{sil[tip_idx[tip]]:.4f}",
                    m.get("domain_no", "NA"), m.get("i_evalue", "NA"),
                ]) + "\n")
                n_rows += 1

    if missing_members:
        eprint(f"WARNING: {missing_members} dereplicated members are missing from metadata and were not written")

    with open(f"{out_prefix}.cluster_report.md", "w", encoding="utf-8", newline="\n") as out:
        out.write(f"# {family} tree-based clustering\n\n")
        out.write(f"- Input tree: `{args.tree}`\n")
        out.write(f"- Tree tips: {len(names)}\n")
        out.write("- Distance: patristic (expected substitutions per site)\n")
        out.write(f"- Hierarchical clustering: {args.linkage} linkage, threshold t = {args.cutoff}\n")
        out.write(f"- Clusters: {n_clusters}; mean silhouette = {np.mean(sil):.4f}\n")
        out.write(f"- Sequences written: {n_rows}\n\n")
        out.write("## Cluster statistics\n\n")
        out.write("| Cluster | Members | Tips | Samples | Mean within | Min between | Silhouette | Monophyletic | UFBoot |\n")
        out.write("|---|---:|---:|---:|---:|---:|---:|---|---|\n")
        for r in stats:
            out.write(f"| {r['cluster']} | {r['n_members']} | {r['n_seqs']} | {r['n_samples']} | "
                      f"{r['mean_within']} | "
                      f"{r['min_between']} | {r['mean_silhouette']} | {r['monophyletic']} | {r['UFBoot']} |\n")
        out.write("\n> `n_members` is the number of sequences after mapping dereplicated members "
                  "back; `n_seqs` is the number of tree tips. The threshold t is a descriptive "
                  "choice; determine it from the silhouette plateau in `.cutoff_scan.tsv` and "
                  "the biological context, rather than treating a single threshold as an "
                  "objective species boundary.\n")

    eprint(f"[{family}] {n_clusters} clusters, {n_rows} sequences -> {out_prefix}.clusters.tsv")
    return 0


if __name__ == "__main__":
    sys.exit(main())
