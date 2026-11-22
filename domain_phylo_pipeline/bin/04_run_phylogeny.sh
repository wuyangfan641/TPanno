#!/usr/bin/env bash
# =============================================================================
#  04_run_phylogeny.sh -- run alignment -> trimming -> tree inference for
#                         every family directory under 03_domains
# =============================================================================
set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../config/config.sh"
source "${SCRIPT_DIR}/lib/log.sh"

[[ -d "$DOMAIN_DIR" ]] || die "$DOMAIN_DIR is missing; run 03_extract_domains.py first"
command -v mafft >/dev/null 2>&1 || die "mafft not found"
command -v trimal >/dev/null 2>&1 || die "trimal not found"
command -v "$IQTREE_BIN" >/dev/null 2>&1 || command -v iqtree2 >/dev/null 2>&1 \
    || command -v iqtree >/dev/null 2>&1 || die "iqtree2/iqtree not found"

mapfile -t fam_dirs < <(find "$DOMAIN_DIR" -mindepth 1 -maxdepth 1 -type d ! -name screen | sort)
[[ ${#fam_dirs[@]} -gt 0 ]] || die "no family directories found under $DOMAIN_DIR"

# Optional allowlist (family_list.txt / family_dirs.txt from 03b_screen_families.py)
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
    [[ -s "${d}/${fam}.domains.faa" ]] || continue
    if [[ -s "${d}/status.tsv" && "$FORCE" -eq 0 ]]; then
        st=$(awk -F'\t' 'NR==2{print $2}' "${d}/status.tsv")
        if [[ "$st" == "ok" ]]; then continue; fi
    fi
    todo+=("$d")
done

log "Total families: ${#fam_dirs[@]}; pending: ${#todo[@]} (parallelism $JOBS, $THREADS_PER_JOB threads per task)"

if [[ ${#todo[@]} -gt 0 ]]; then
    printf '%s\0' "${todo[@]}" | \
        xargs -0 -n1 -P "$JOBS" bash -c 'bash "$1/04_one_family.sh" "$2"' _ "$SCRIPT_DIR" \
        || warn "some families failed tree inference; see each family's status.tsv"
fi

# ---------- Status summary ----------------------------------------------------
ALL_STATUS="${PHYLO_DIR}/all_status.tsv"
mkdir -p "$PHYLO_DIR"
{
    printf 'family\tstatus\tn_input\tn_used\tn_aln_cols\tn_trim_cols\ttrim_fraction\tbest_model\tnotes\n'
    for d in "${fam_dirs[@]}"; do
        if [[ -s "${d}/status.tsv" ]]; then tail -n +2 "${d}/status.tsv"; fi
    done
} > "$ALL_STATUS"

n_ok=$(awk -F'\t' 'NR>1 && $2=="ok"' "$ALL_STATUS" | wc -l)
n_skip=$(awk -F'\t' 'NR>1 && $2=="skipped"' "$ALL_STATUS" | wc -l)
n_fail=$(awk -F'\t' 'NR>1 && $2=="failed"' "$ALL_STATUS" | wc -l)
log "Tree inference summary: ok=$n_ok, skipped=$n_skip, failed=$n_fail -> $ALL_STATUS"

if [[ "$n_fail" -gt 0 ]]; then
    warn "$n_fail families failed:"
    awk -F'\t' 'NR>1 && $2=="failed"{print "  - "$1": "$9}' "$ALL_STATUS" >&2
fi
