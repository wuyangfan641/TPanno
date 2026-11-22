#!/usr/bin/env bash
# =============================================================================
#  02_run_search.sh -- domain search with HMMER
#
#  Two modes (SEARCH_MODE in config.sh):
#    hmmsearch : every .hmm in HMM_DIR against every sample
#    hmmscan   : every sample against an hmmpress'ed database (e.g. Pfam-A)
#
#  Output: results/02_search/*.domtblout
#  Note: in domtblout, the target is always the protein sequence and the query
#        is always the HMM model, so 03_extract_domains.py does not need to
#        distinguish the two modes.
# =============================================================================
set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../config/config.sh"
source "${SCRIPT_DIR}/lib/log.sh"

SAMPLES="${PREP_DIR}/samples.clean.tsv"
[[ -s "$SAMPLES" ]] || die "$SAMPLES is missing; run 00_prepare_inputs.py first"
mkdir -p "$SEARCH_DIR" "$LOG_DIR"

if [[ "$SEARCH_MODE" == "hmmscan" ]]; then
    command -v hmmscan >/dev/null 2>&1 || die "hmmscan not found; activate the conda environment"
else
    command -v hmmsearch >/dev/null 2>&1 || die "hmmsearch not found; activate the conda environment"
fi

# Check whether every model in an HMM file has GA thresholds
all_models_have_ga() {
    local f="$1" n_model n_ga
    n_model=$(grep -c '^HMMER3' "$f" 2>/dev/null || true)
    n_ga=$(grep -c '^GA' "$f" 2>/dev/null || true)
    [[ "$n_model" -gt 0 && "$n_model" -eq "$n_ga" ]]
}

run_hmmsearch() {
    local hmm="$1" sample="$2" fa="$3" out="$4"
    local -a cut_args=()
    case "$HMM_CUT_MODE" in
        ga)
            if all_models_have_ga "$hmm"; then
                cut_args=(--cut_ga)
            else
                cut_args=(-E "$HMM_EVALUE_MAX")
            fi
            ;;
        nc)     cut_args=(--cut_nc) ;;
        evalue) cut_args=(-E "$HMM_EVALUE_MAX") ;;
        *)      cut_args=(-E "$HMM_EVALUE_MAX") ;;
    esac
    hmmsearch --cpu "$HMMSEARCH_CPU" --noali "${cut_args[@]}" \
        --domtblout "$out" -o /dev/null "$hmm" "$fa" \
        || { rm -f "$out"; return 1; }
}

run_hmmscan() {
    local sample="$1" fa="$2" out="$3"
    local -a cut_args=()
    case "$HMM_CUT_MODE" in
        ga)
            if [[ "${HMM_DB_HAS_GA:-0}" -eq 1 ]]; then
                cut_args=(--cut_ga)
            else
                cut_args=(-E "$HMM_EVALUE_MAX")
            fi
            ;;
        nc)     cut_args=(--cut_nc) ;;
        evalue) cut_args=(-E "$HMM_EVALUE_MAX") ;;
        *)      cut_args=(-E "$HMM_EVALUE_MAX") ;;
    esac
    hmmscan --cpu "$HMMSEARCH_CPU" --noali "${cut_args[@]}" \
        --domtblout "$out" -o /dev/null "$HMM_SCAN_DB" "$fa" \
        || { rm -f "$out"; return 1; }
}

# Check GA availability once in the parent process. Checking it inside every
# parallel hmmscan task would re-read the entire Pfam-A.hmm file per sample.
HMM_DB_HAS_GA=0
if [[ "$SEARCH_MODE" == "hmmscan" && "$HMM_CUT_MODE" == "ga" ]]; then
    if all_models_have_ga "$HMM_SCAN_DB"; then
        HMM_DB_HAS_GA=1
    fi
fi
export HMM_DB_HAS_GA

jobbuf=$(mktemp)
trap 'rm -f "$jobbuf"' EXIT

n_jobs=0
if [[ "$SEARCH_MODE" == "hmmsearch" ]]; then
    mapfile -t hmm_files < <(find "$HMM_DIR" -maxdepth 1 -type f -name '*.hmm' | sort)
    [[ ${#hmm_files[@]} -gt 0 ]] || die "no .hmm files found in $HMM_DIR"
    log "hmmsearch mode: ${#hmm_files[@]} HMM files x $(($(wc -l < "$SAMPLES") - 1)) samples"

    while IFS=$'\t' read -r sample fa _rest; do
        sample="${sample#$'\xEF\xBB\xBF'}"
        fa="${fa%$'\r'}"
        if [[ -z "$sample" ]]; then continue; fi
        if [[ "$sample" == "sample" && "$fa" == "fasta" ]]; then continue; fi
        for hmm in "${hmm_files[@]}"; do
            hbase=$(basename "$hmm" .hmm | tr -c 'A-Za-z0-9_.-' '_')
            out="${SEARCH_DIR}/${hbase}__${sample}.domtblout"
            if [[ -s "$out" && "$FORCE" -eq 0 ]]; then
                continue
            fi
            printf '%s\0%s\0%s\0%s\0' "$hmm" "$sample" "$fa" "$out" >> "$jobbuf"
            n_jobs=$((n_jobs + 1))
        done
    done < "$SAMPLES"
else
    [[ -s "${HMM_SCAN_DB}.h3m" ]] || die "HMM database is not hmmpress'ed: ${HMM_SCAN_DB}.h3m is missing"
    log "hmmscan mode: database $(basename "$HMM_SCAN_DB") x $(($(wc -l < "$SAMPLES") - 1)) samples"

    while IFS=$'\t' read -r sample fa _rest; do
        sample="${sample#$'\xEF\xBB\xBF'}"
        fa="${fa%$'\r'}"
        if [[ -z "$sample" ]]; then continue; fi
        if [[ "$sample" == "sample" && "$fa" == "fasta" ]]; then continue; fi
        out="${SEARCH_DIR}/${sample}.domtblout"
        if [[ -s "$out" && "$FORCE" -eq 0 ]]; then
            continue
        fi
        printf '%s\0%s\0%s\0' "$sample" "$fa" "$out" >> "$jobbuf"
        n_jobs=$((n_jobs + 1))
    done < "$SAMPLES"
fi

if [[ "$n_jobs" -eq 0 ]]; then
    log "All searches are already complete; nothing to rerun"
    exit 0
fi

log "Running $n_jobs search tasks with parallelism $JOBS"

if [[ "$SEARCH_MODE" == "hmmsearch" ]]; then
    export -f all_models_have_ga run_hmmsearch
    export HMM_CUT_MODE HMM_EVALUE_MAX HMMSEARCH_CPU
    # shellcheck disable=SC2016
    xargs -0 -n4 -P "$JOBS" bash -c 'run_hmmsearch "$1" "$2" "$3" "$4"' _ < "$jobbuf" \
        || die "some hmmsearch tasks failed; check results/logs and whether the *.domtblout files are empty"
else
    export -f run_hmmscan
    export HMM_CUT_MODE HMM_EVALUE_MAX HMMSEARCH_CPU HMM_SCAN_DB HMM_DB_HAS_GA
    # shellcheck disable=SC2016
    xargs -0 -n3 -P "$JOBS" bash -c 'run_hmmscan "$1" "$2" "$3"' _ < "$jobbuf" \
        || die "some hmmscan tasks failed"
fi

log "Search complete: $(find "$SEARCH_DIR" -maxdepth 1 -name '*.domtblout' | wc -l) domtblout files"
