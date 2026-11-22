#!/usr/bin/env bash
# =============================================================================
#  run_pfam_scan.sh -- staged entry point for an HMM-guided full Pfam-A scan
#
#  Why the run is staged:
#    1) a full hmmscan is expensive, so extract all domain families first;
#    2) 03b_screen_families.py filters by sample coverage, copy number,
#       domain coverage, and clan redundancy;
#    3) only selected families are sent to tree inference and clustering,
#       avoiding thousands of queued trees.
#
#  Examples:
#    bash run_pfam_scan.sh
#    MIN_SAMPLES=22 MIN_SEQS=22 MAX_FAMILIES=300 bash run_pfam_scan.sh
#    # 104 cores / 800 GB recommended:
#    JOBS=24 HMMSEARCH_CPU=4 THREADS_PER_JOB=4 \
#      MIN_SAMPLES=22 MIN_SEQS=22 MAX_FAMILIES=300 bash run_pfam_scan.sh
#
#  Key environment variables (all can be overridden on the command line):
#    HMM_SCAN_DB       hmmpress'ed Pfam-A.hmm
#    HMMSEARCH_CPU     CPUs per hmmscan process (default: 8)
#    JOBS              parallel samples/families (default: 2; watch hmmscan memory)
#    THREADS_PER_JOB   threads per IQ-TREE run (default: 8)
#    MIN_SAMPLES       with 27 species: majority >= 14, conserved >= 22
#    MIN_SEQS          minimum number of sequences entering tree inference
#    MAX_FAMILIES      first-round cap, 200-500 recommended; 0 = no limit
# =============================================================================
set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

export SEARCH_MODE="hmmscan"
export HMM_SCAN_DB="${HMM_SCAN_DB:-${SCRIPT_DIR}/01_hmms/Pfam-A.hmm}"
export HMMSEARCH_CPU="${HMMSEARCH_CPU:-8}"
export JOBS="${JOBS:-2}"
export THREADS_PER_JOB="${THREADS_PER_JOB:-8}"
export MIN_SAMPLES="${MIN_SAMPLES:-22}"
export MIN_SEQS="${MIN_SEQS:-22}"
export MAX_FAMILIES="${MAX_FAMILIES:-300}"
export MIN_MEDIAN_COV="${MIN_MEDIAN_COV:-0.75}"
export MIN_MEDIAN_DOMAIN_LEN="${MIN_MEDIAN_DOMAIN_LEN:-30}"
export MAX_MEDIAN_COPIES="${MAX_MEDIAN_COPIES:-0}"
export MIN_SINGLE_COPY_FRAC="${MIN_SINGLE_COPY_FRAC:-0.0}"
export CLUSTER_LINKAGE="${CLUSTER_LINKAGE:-complete}"
export CLUSTER_CUTOFF="${CLUSTER_CUTOFF:-0.50}"
export MAX_REPORT_ROWS="${MAX_REPORT_ROWS:-200}"
export COLLAPSE_CLANS="${COLLAPSE_CLANS:-1}"

# Load the complete parameter set so 06_make_report.py can record it even
# though the child run_pipeline.sh process cannot export back to this shell.
# shellcheck disable=SC1091
source "${SCRIPT_DIR}/config/config.sh"

die() { echo "ERROR: $*" >&2; exit 1; }
log() { printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" >&2; }

[[ -s "$HMM_SCAN_DB" ]] || die "$HMM_SCAN_DB not found; run hmmpress on Pfam-A.hmm first"
[[ -s "${HMM_SCAN_DB}.h3m" ]] || die "Pfam-A is not hmmpress'ed: ${HMM_SCAN_DB}.h3m is missing"
command -v hmmscan >/dev/null 2>&1 || die "hmmscan not found"
command -v mafft >/dev/null 2>&1 || die "mafft not found"
command -v trimal >/dev/null 2>&1 || die "trimal not found"
command -v "$IQTREE_BIN" >/dev/null 2>&1 || command -v iqtree2 >/dev/null 2>&1 \
    || command -v iqtree >/dev/null 2>&1 || die "iqtree2/iqtree not found"

log "Full Pfam scan parameters: MIN_SAMPLES=$MIN_SAMPLES MIN_SEQS=$MIN_SEQS MAX_FAMILIES=$MAX_FAMILIES"
log "Resources: JOBS=$JOBS, HMMSEARCH_CPU=$HMMSEARCH_CPU, THREADS_PER_JOB=$THREADS_PER_JOB"

# ---------- 1. Input QC + sample table ---------------------------------------
# Do not overwrite an existing (possibly hand-edited) sample table.
if [[ -s "$SAMPLES_TSV" ]]; then
    log "Using the existing sample table: $SAMPLES_TSV"
    bash run_pipeline.sh --from prepare --to prepare
else
    bash run_pipeline.sh --init
fi

# ---------- 2. Full hmmscan + domain extraction ------------------------------
# Lenient thresholds are used deliberately here so that all candidates are
# written to disk. The real 22/22 screening happens in 03b; otherwise 03b
# would never see families discarded early and could not report their copy
# number or coverage.
FAMILY_LIST='' MIN_SAMPLES=1 MIN_SEQS=4 bash run_pipeline.sh --from search --to domains \
    --jobs "$JOBS" --threads "$THREADS_PER_JOB"

# ---------- 3. Family screening ----------------------------------------------
SCREEN_DIR="${DOMAIN_DIR}/screen"
if [[ -z "${CLAN_MAP:-}" && -f "${PROJECT_DIR}/01_hmms/Pfam-A.clans.tsv" ]]; then
    CLAN_MAP="${PROJECT_DIR}/01_hmms/Pfam-A.clans.tsv"
fi
clan_args=()
if [[ -n "${CLAN_MAP:-}" ]]; then
    [[ -f "$CLAN_MAP" ]] || die "CLAN_MAP does not exist: $CLAN_MAP"
    clan_args=(--clan-map "$CLAN_MAP")
fi
collapse_args=()
if [[ "${COLLAPSE_CLANS:-1}" -eq 1 ]]; then
    if [[ ${#clan_args[@]} -eq 0 ]]; then
        log "WARNING: COLLAPSE_CLANS=1 but no clan map was found; clan dereplication is skipped"
    else
        collapse_args=(--collapse-clans)
    fi
fi

log "Screening families -> $SCREEN_DIR"
python3 bin/03b_screen_families.py \
    --domain-dir "$DOMAIN_DIR" \
    --outdir "$SCREEN_DIR" \
    --min-samples "$MIN_SAMPLES" \
    --min-seqs "$MIN_SEQS" \
    --min-median-cov "$MIN_MEDIAN_COV" \
    --min-median-domain-len "$MIN_MEDIAN_DOMAIN_LEN" \
    --max-median-copies "$MAX_MEDIAN_COPIES" \
    --min-single-copy-frac "$MIN_SINGLE_COPY_FRAC" \
    --max-families "$MAX_FAMILIES" \
    "${clan_args[@]}" "${collapse_args[@]}"

# ---------- 4. Tree inference + clustering for selected families -------------
export FAMILY_LIST="${SCREEN_DIR}/family_dirs.txt"
[[ -s "$FAMILY_LIST" ]] || die "screening produced an empty family list: $FAMILY_LIST"
log "Selected $(wc -l < "$FAMILY_LIST") families; starting tree inference and clustering"
bash run_pipeline.sh --from phylo --to cluster \
    --jobs "$JOBS" --threads "$THREADS_PER_JOB"

# ---------- 5. Summary report ------------------------------------------------
python3 bin/06_make_report.py --results "$RESULTS"

log "Done. Overall report: ${REPORT_DIR}/report.md"
log "Screening details: $SCREEN_DIR/family_screen.tsv"
log "Selected families: $SCREEN_DIR/family_list.txt"
