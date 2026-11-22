#!/usr/bin/env bash
# =============================================================================
#  run_pipeline.sh -- one-command run:
#      input QC -> HMMER -> domain extraction -> MAFFT/trimAl -> IQ-TREE
#      -> tree-based clustering -> summary report
#
#  Usage:
#    bash run_pipeline.sh --init                 # scan 00_raw/, build samples.tsv, run QC
#    bash run_pipeline.sh                        # full pipeline
#    bash run_pipeline.sh --jobs 8 --threads 4   # 8 tasks in parallel, 4 threads each
#    bash run_pipeline.sh --from phylo --to cluster
#    bash run_pipeline.sh --dry-run              # print the steps without executing them
# =============================================================================
set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_FILE="${SCRIPT_DIR}/config/config.sh"
JOBS_OVERRIDE=""
THREADS_OVERRIDE=""
FROM_STEP=""
TO_STEP=""
INIT=0
DRY_RUN=0
FORCE_FLAG=0

usage() {
    cat <<'EOF'
Usage: bash run_pipeline.sh [options]

  --init              scan 00_raw/, generate config/samples.tsv, and run input QC
  --jobs N            number of concurrent tasks (default: JOBS in config.sh)
  --threads N         threads per IQ-TREE task (default: THREADS_PER_JOB in config.sh)
  --from STEP         start at STEP: prepare | search | domains | phylo | cluster | report
  --to STEP           stop after STEP (inclusive)
  --config FILE       configuration file (default: config/config.sh)
  --force             overwrite existing results (default: resume from checkpoints)
  --dry-run           print the commands without executing them
  -h, --help          show this help
EOF
}

require_value() {
    [[ $# -ge 2 && -n "$2" ]] || { echo "ERROR: option $1 requires a value" >&2; exit 2; }
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --config)   require_value "$@"; CONFIG_FILE="$2"; shift 2 ;;
        --jobs)     require_value "$@"; JOBS_OVERRIDE="$2"; shift 2 ;;
        --threads)  require_value "$@"; THREADS_OVERRIDE="$2"; shift 2 ;;
        --from)     require_value "$@"; FROM_STEP="$2"; shift 2 ;;
        --to)       require_value "$@"; TO_STEP="$2"; shift 2 ;;
        --init)     INIT=1; shift ;;
        --force)    FORCE_FLAG=1; shift ;;
        --dry-run)  DRY_RUN=1; shift ;;
        -h|--help)  usage; exit 0 ;;
        *) echo "ERROR: unknown option: $1" >&2; usage >&2; exit 2 ;;
    esac
done

[[ -f "$CONFIG_FILE" ]] || { echo "ERROR: configuration file not found: $CONFIG_FILE" >&2; exit 2; }
# shellcheck disable=SC1090
source "$CONFIG_FILE"
source "${SCRIPT_DIR}/bin/lib/log.sh"

if [[ -n "$JOBS_OVERRIDE" ]]; then JOBS="$JOBS_OVERRIDE"; fi
if [[ -n "$THREADS_OVERRIDE" ]]; then THREADS_PER_JOB="$THREADS_OVERRIDE"; fi
if [[ "$FORCE_FLAG" -eq 1 ]]; then FORCE=1; fi
export JOBS THREADS_PER_JOB FORCE

mkdir -p "$RESULTS" "$LOG_DIR" "$REPORT_DIR"
LOG_FILE="${LOG_DIR}/run_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG_FILE") 2>&1

log "Project directory: $PROJECT_DIR"
log "Configuration file: $CONFIG_FILE"
log "Parallelism: JOBS=$JOBS, THREADS_PER_JOB=$THREADS_PER_JOB, SEARCH_MODE=$SEARCH_MODE"

# ---------- Dependency checks ------------------------------------------------
need() { command -v "$1" >/dev/null 2>&1 || die "command not found: $1 (activate the matching conda environment)"; }
check_python_modules() {
    local modules="$1"
    python3 - "$modules" <<'PY' || die "missing Python modules: ${modules} (see env/environment.yml)"
import importlib
import sys
for name in sys.argv[1].split(","):
    importlib.import_module(name.strip())
PY
}

STEPS=(prepare search domains phylo cluster report)
if [[ "$INIT" -eq 1 ]]; then
    STEPS=(prepare)
fi

# Select the step interval
start=0; end=${#STEPS[@]}
if [[ -n "$FROM_STEP" ]]; then
    found=0
    for i in "${!STEPS[@]}"; do
        if [[ "${STEPS[$i]}" == "$FROM_STEP" ]]; then
            start=$i; found=1; break
        fi
    done
    [[ "$found" -eq 1 ]] || die "unknown --from step: $FROM_STEP (valid: ${STEPS[*]})"
fi
if [[ -n "$TO_STEP" ]]; then
    found=0
    for i in "${!STEPS[@]}"; do
        if [[ "${STEPS[$i]}" == "$TO_STEP" ]]; then
            end=$((i + 1)); found=1; break
        fi
    done
    [[ "$found" -eq 1 ]] || die "unknown --to step: $TO_STEP (valid: ${STEPS[*]})"
fi
(( start < end )) || die "--from must come before --to"

# Check only the dependencies needed by the selected steps
for ((i=start; i<end; i++)); do
    case "${STEPS[$i]}" in
        prepare) need python3; check_python_modules Bio ;;
        search)
            if [[ "$SEARCH_MODE" == "hmmscan" ]]; then need hmmscan; else need hmmsearch; fi ;;
        domains) need python3; check_python_modules Bio ;;
        phylo)   need python3; need mafft; need trimal
                 command -v "$IQTREE_BIN" >/dev/null 2>&1 || command -v iqtree2 >/dev/null 2>&1 \
                     || command -v iqtree >/dev/null 2>&1 || die "iqtree2/iqtree not found" ;;
        cluster) need python3; check_python_modules Bio,numpy,scipy ;;
        report)  need python3 ;;
    esac
done

# ---------- Software versions ------------------------------------------------
VERSION_FILE="${REPORT_DIR}/software_versions.txt"
{
    echo "# Generated: $(date '+%Y-%m-%d %H:%M:%S %Z')"
    echo "# Host: $(uname -a)"
    echo "python3: $(python3 --version 2>&1)"
    if command -v hmmsearch >/dev/null 2>&1; then
        echo "hmmsearch: $(hmmsearch -h 2>&1 | head -n 2 | tr '\n' ' ')"
    fi
    if command -v hmmscan >/dev/null 2>&1; then
        echo "hmmscan: $(hmmscan -h 2>&1 | head -n 2 | tr '\n' ' ')"
    fi
    echo "mafft: $(mafft --version 2>&1)"
    echo "trimal: $(trimal --version 2>&1 | head -n 1)"
    if command -v iqtree2 >/dev/null 2>&1; then
        echo "iqtree2: $(iqtree2 --version 2>&1 | head -n 1)"
    elif command -v iqtree >/dev/null 2>&1; then
        echo "iqtree: $(iqtree --version 2>&1 | head -n 1)"
    fi
    if command -v cd-hit >/dev/null 2>&1; then
        echo "cd-hit: $(cd-hit 2>&1 | head -n 2 | tr '\n' ' ')"
    fi
} > "$VERSION_FILE" 2>&1 || true

# ---------- Execution --------------------------------------------------------
run_step() {
    local name="$1"; shift
    local t0=$SECONDS
    log "================ STEP: $name ================"
    if [[ "$DRY_RUN" -eq 1 ]]; then
        log "[dry-run] $*"
        return 0
    fi
    "$@"
    log "================ STEP: $name done ($((SECONDS - t0))s) ================"
}

PREP_CMD=(python3 "${SCRIPT_DIR}/bin/00_prepare_inputs.py"
          --raw-dir "$RAW_DIR" --outdir "$PREP_DIR")
if [[ "$INIT" -eq 1 ]]; then
    PREP_CMD+=(--write-samples "$SAMPLES_TSV")
elif [[ -f "$SAMPLES_TSV" ]]; then
    PREP_CMD+=(--samples "$SAMPLES_TSV")
fi
if [[ "$FORCE" -eq 1 ]]; then PREP_CMD+=(--force); fi

DOMAIN_CMD=(python3 "${SCRIPT_DIR}/bin/03_extract_domains.py"
            --search-dir "$SEARCH_DIR"
            --samples "${PREP_DIR}/samples.clean.tsv"
            --outdir "$DOMAIN_DIR")
if [[ -n "$FAMILY_LIST" ]]; then DOMAIN_CMD+=(--family-list "$FAMILY_LIST"); fi

for ((i=start; i<end; i++)); do
    case "${STEPS[$i]}" in
        prepare) run_step prepare "${PREP_CMD[@]}" ;;
        search)  run_step search  bash "${SCRIPT_DIR}/bin/02_run_search.sh" ;;
        domains) run_step domains "${DOMAIN_CMD[@]}" ;;
        phylo)   run_step phylo   bash "${SCRIPT_DIR}/bin/04_run_phylogeny.sh" ;;
        cluster) run_step cluster bash "${SCRIPT_DIR}/bin/05_run_clustering.sh" ;;
        report)  run_step report  python3 "${SCRIPT_DIR}/bin/06_make_report.py" --results "$RESULTS" ;;
    esac
done

if [[ "$INIT" -eq 1 ]]; then
    cat <<EOF

[init] Sample table: ${SAMPLES_TSV}
       Check that the sample names and paths are correct (27 rows) and edit
       them manually if necessary. When ready, run: bash run_pipeline.sh
EOF
else
    log "All steps completed. Main report: ${REPORT_DIR}/report.md"
    log "Log: ${LOG_FILE}"
fi
