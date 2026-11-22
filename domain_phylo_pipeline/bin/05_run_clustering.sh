#!/usr/bin/env bash
# =============================================================================
#  05_run_clustering.sh -- patristic-distance hierarchical clustering for
#                          every family tree
# =============================================================================
set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../config/config.sh"
source "${SCRIPT_DIR}/lib/log.sh"

[[ -d "$PHYLO_DIR" ]] || die "$PHYLO_DIR is missing"
python3 - <<'PY' || die "incomplete Python dependencies (Biopython, NumPy, and SciPy are required)"
import Bio, numpy, scipy
PY

mkdir -p "$CLUSTER_DIR"

mapfile -t fam_dirs < <(find "$PHYLO_DIR" -mindepth 1 -maxdepth 1 -type d | sort)
[[ ${#fam_dirs[@]} -gt 0 ]] || die "no family directories found under $PHYLO_DIR"

declare -A SELECTED=()
if [[ -n "$FAMILY_LIST" ]]; then
    [[ -f "$FAMILY_LIST" ]] || die "FAMILY_LIST does not exist: $FAMILY_LIST"
    while IFS= read -r line; do
        line="${line%$'\r'}"
        line="${line#$'\xEF\xBB\xBF'}"
        line="${line%%#*}"
        read -r fam _rest <<< "$line" || true
        [[ -z "$fam" ]] && continue
        fam_prefix="${fam%%.*}"
        SELECTED["$fam"]=1
        SELECTED["$fam_prefix"]=1
    done < "$FAMILY_LIST"
    log "Family allowlist: ${#SELECTED[@]} keys (from $FAMILY_LIST)"
fi

todo=()
for d in "${fam_dirs[@]}"; do
    fam=$(basename "$d")
    fam_prefix="${fam%%.*}"
    if [[ ${#SELECTED[@]} -gt 0 && -z "${SELECTED[$fam]:-}" && -z "${SELECTED[$fam_prefix]:-}" ]]; then
        continue
    fi
    tree="${d}/iqtree.treefile"
    meta="${DOMAIN_DIR}/${fam}/metadata.tsv"
    [[ -s "$tree" && -s "$meta" ]] || continue
    [[ -s "${d}/status.tsv" ]] || continue
    st=$(awk -F'\t' 'NR==2{print $2}' "${d}/status.tsv")
    [[ "$st" == "ok" ]] || continue
    if [[ -s "${CLUSTER_DIR}/${fam}.clusters.tsv" && "$FORCE" -eq 0 ]]; then
        continue
    fi
    todo+=("$d")
done

log "Families pending clustering: ${#todo[@]} / ${#fam_dirs[@]} (parallelism $JOBS)"

run_one() {
    local d="$1" script_dir="$2"
    local fam; fam=$(basename "$d")
    local -a extra=()
    if [[ -s "${d}/derep_map.tsv" ]]; then
        extra=(--derep-map "${d}/derep_map.tsv")
    fi
    python3 "${script_dir}/05_cluster_tree.py" \
        --tree "${d}/iqtree.treefile" \
        --metadata "${DOMAIN_DIR}/${fam}/metadata.tsv" \
        --out-prefix "${CLUSTER_DIR}/${fam}" \
        --cutoff "$CLUSTER_CUTOFF" \
        --linkage "$CLUSTER_LINKAGE" \
        --min-size "$MIN_CLUSTER_SIZE" \
        --scan "$CLUSTER_SCAN" \
        "${extra[@]}"
}
export -f run_one
export DOMAIN_DIR CLUSTER_DIR CLUSTER_CUTOFF CLUSTER_LINKAGE MIN_CLUSTER_SIZE CLUSTER_SCAN

if [[ ${#todo[@]} -gt 0 ]]; then
    printf '%s\0' "${todo[@]}" | \
        xargs -0 -n1 -P "$JOBS" bash -c 'run_one "$2" "$1"' _ "$SCRIPT_DIR" \
        || warn "some families failed clustering"
fi

# ---------- Merge -------------------------------------------------------------
printf 'family\tcluster\tn_seqs\tn_members\tn_samples\tsamples\tsample_counts\tmean_within\tmax_within\tmin_between\tmean_silhouette\tmonophyletic\tSH_aLRT\tUFBoot\trepresentative\trepresentative_sample\tpasses_min_size\n' \
    > "${CLUSTER_DIR}/all_cluster_summary.tsv"
for f in "${CLUSTER_DIR}"/*.cluster_summary.tsv; do
    [[ -s "$f" ]] || continue
    tail -n +2 "$f" >> "${CLUSTER_DIR}/all_cluster_summary.tsv"
done

printf 'family\tcluster\tcluster_size\tseq_id\tsample\tprotein_id\tis_representative\ttree_tip\tsilhouette\tdomain_no\ti_evalue\n' \
    > "${CLUSTER_DIR}/all_clusters.tsv"
for f in "${CLUSTER_DIR}"/*.clusters.tsv; do
    [[ -s "$f" ]] || continue
    tail -n +2 "$f" >> "${CLUSTER_DIR}/all_clusters.tsv"
done

log "Clustering summary: $(($(wc -l < "${CLUSTER_DIR}/all_cluster_summary.tsv") - 1)) clusters, $(($(wc -l < "${CLUSTER_DIR}/all_clusters.tsv") - 1)) sequence assignments"
