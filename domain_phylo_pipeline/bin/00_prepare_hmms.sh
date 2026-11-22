#!/usr/bin/env bash
# =============================================================================
#  00_prepare_hmms.sh -- extract selected HMMs from Pfam-A.hmm by name
#
#  Example:
#    cat > pfam_ids.txt <<EOF
#    PF00069.28   # Pkinase
#    PF07714.17   # PK_Tyr_Ser-Thr
#    EOF
#    bash bin/00_prepare_hmms.sh --pfam /db/Pfam-A.hmm --names pfam_ids.txt \
#         --out 01_hmms/selected.hmm --press
#
#  Note: hmmsearch mode is intended for a selected set of families (tens to
#        hundreds of models). To scan the full Pfam library, use
#        SEARCH_MODE=hmmscan with an hmmpress'ed Pfam-A.hmm instead.
# =============================================================================
set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../config/config.sh"
source "${SCRIPT_DIR}/lib/log.sh"

PFAM=""
NAMES=""
OUT=""
PRESS=0

usage() {
    sed -n '2,16p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
}

require_value() {
    [[ $# -ge 2 && -n "$2" ]] || die "option $1 requires a value"
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --pfam)   require_value "$@"; PFAM="$2"; shift 2 ;;
        --names)  require_value "$@"; NAMES="$2"; shift 2 ;;
        --out)    require_value "$@"; OUT="$2"; shift 2 ;;
        --press)  PRESS=1; shift ;;
        -h|--help) usage; exit 0 ;;
        *) die "unknown option: $1" ;;
    esac
done
[[ -n "$PFAM" && -f "$PFAM" ]] || die "use --pfam to specify Pfam-A.hmm"
[[ -n "$NAMES" && -f "$NAMES" ]] || die "use --names to specify the model-name list"
[[ -n "$OUT" ]] || die "use --out to specify the output .hmm file"

command -v hmmfetch >/dev/null 2>&1 || die "hmmfetch not found"
mkdir -p "$(dirname "$OUT")"

tmp_names=$(mktemp)
trap 'rm -f "$tmp_names"' EXIT
grep -v '^[[:space:]]*#' "$NAMES" | sed '1s/^\xEF\xBB\xBF//' | awk 'NF{print $1}' > "$tmp_names"
n=$(wc -l < "$tmp_names")
[[ "$n" -gt 0 ]] || die "no valid model names found in $NAMES"

log "Extracting $n models from $PFAM -> $OUT"
hmmfetch -f "$PFAM" "$tmp_names" > "$OUT"
[[ -s "$OUT" ]] || die "hmmfetch produced an empty file (do the model names exactly match the NAME fields in Pfam?)"

if [[ "$PRESS" -eq 1 ]]; then
    hmmpress -f "$OUT"
    log "hmmpress completed: ${OUT}.h3m"
fi
log "Done. Point HMM_DIR in config.sh to the directory containing this file."
