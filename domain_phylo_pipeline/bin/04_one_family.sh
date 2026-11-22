#!/usr/bin/env bash
# =============================================================================
#  04_one_family.sh -- one family: MAFFT/hmmalign -> trimAl -> IQ-TREE
#  Called in parallel by 04_run_phylogeny.sh; can also be run directly:
#      bash bin/04_one_family.sh results/03_domains/PF00069.28
# =============================================================================
set -uo pipefail

FAM_DIR="${1:?usage: 04_one_family.sh <family_dir>}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../config/config.sh"
source "${SCRIPT_DIR}/lib/log.sh"

[[ -d "$FAM_DIR" ]] || { echo "ERROR: family directory not found: $FAM_DIR" >&2; exit 1; }
FAM_DIR="$(cd "$FAM_DIR" && pwd)"
FAM="$(basename "$FAM_DIR")"
FAA="${FAM_DIR}/${FAM}.domains.faa"
STATUS="${FAM_DIR}/status.tsv"
ALN="${FAM_DIR}/${FAM}.aln.fa"
TRIM="${FAM_DIR}/${FAM}.trim.fa"
PREFIX="${FAM_DIR}/iqtree"

mkdir -p "$FAM_DIR"

write_status() {
    local status="$1" n_in="$2" n_used="$3" aln_cols="$4" trim_cols="$5" model="$6" notes="$7"
    local trim_frac="NA"
    if [[ "$aln_cols" =~ ^[0-9]+$ && "$aln_cols" -gt 0 && "$trim_cols" =~ ^[0-9]+$ ]]; then
        trim_frac=$(awk -v a="$trim_cols" -v b="$aln_cols" 'BEGIN{printf "%.3f", a/b}')
    fi
    {
        printf 'family\tstatus\tn_input\tn_used\tn_aln_cols\tn_trim_cols\ttrim_fraction\tbest_model\tnotes\n'
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "$FAM" "$status" "$n_in" "$n_used" "$aln_cols" "$trim_cols" "$trim_frac" "$model" "$notes"
    } > "$STATUS"
}

fail() {
    local msg="$1"
    warn "[$FAM] $msg"
    write_status "failed" "${N_IN:-NA}" "${N_USED:-NA}" "NA" "NA" "NA" "$msg"
    exit 1
}

skip() {
    local msg="$1"
    log "[$FAM] skip: $msg"
    write_status "skipped" "${N_IN:-NA}" "${N_USED:-NA}" "NA" "NA" "NA" "$msg"
    exit 0
}

[[ -s "$FAA" ]] || fail "$FAA is missing"
if [[ -s "$STATUS" && "$FORCE" -eq 0 ]]; then
    st=$(awk -F'\t' 'NR==2{print $2}' "$STATUS")
    if [[ "$st" == "ok" ]]; then exit 0; fi
fi

N_IN=$(awk '/^>/{n++} END{print n+0}' "$FAA")
N_USED="$N_IN"
if [[ "$N_IN" -lt "$MIN_SEQS" ]]; then
    skip "sequence count $N_IN < MIN_SEQS=$MIN_SEQS"
fi

# ---------- Optional: CD-HIT dereplication (very large families only) --------
ALN_IN="$FAA"
if [[ "$DEREPLICATE" -eq 1 && "$DEREP_MAX_SEQS" -gt 0 && "$N_IN" -gt "$DEREP_MAX_SEQS" ]]; then
    if command -v cd-hit >/dev/null 2>&1; then
        log "[$FAM] n=$N_IN > $DEREP_MAX_SEQS; running CD-HIT dereplication (-c $DEREP_ID)"
        NR="${FAM_DIR}/${FAM}.nr.faa"
        if ! cd-hit -i "$FAA" -o "$NR" -c "$DEREP_ID" -n 5 -T 1 -M 0 -d 0 \
                > "${FAM_DIR}/cdhit.log" 2>&1; then
            fail "cd-hit failed"
        fi
        if ! python3 "${SCRIPT_DIR}/lib/cdhit_map.py" \
                --clstr "${NR}.clstr" --out "${FAM_DIR}/derep_map.tsv" \
                2>> "${FAM_DIR}/cdhit.log"; then
            fail "failed to parse the cd-hit .clstr file"
        fi
        ALN_IN="$NR"
        N_USED=$(awk '/^>/{n++} END{print n+0}' "$NR")
    else
        warn "[$FAM] cd-hit is not installed; skipping dereplication (n=$N_IN, this may be slow)"
    fi
fi

# ---------- Multiple sequence alignment --------------------------------------
case "$ALIGNER" in
    hmmalign)
        if [[ -z "$HMM_SOURCE" ]]; then
            warn "[$FAM] ALIGNER=hmmalign but HMM_SOURCE is empty; falling back to MAFFT"
            ALIGNER="mafft"
        elif command -v hmmfetch >/dev/null 2>&1 \
                && hmmfetch "$HMM_SOURCE" "$FAM" > "${FAM_DIR}/${FAM}.hmm" 2>/dev/null \
                && [[ -s "${FAM_DIR}/${FAM}.hmm" ]]; then
            log "[$FAM] hmmalign (HMM coordinate frame)"
            if ! hmmalign --trim --outformat afa -o "$ALN" "${FAM_DIR}/${FAM}.hmm" "$ALN_IN" \
                    > "${FAM_DIR}/hmmalign.log" 2>&1; then
                warn "[$FAM] hmmalign failed; falling back to MAFFT"
                ALIGNER="mafft"
            else
                sed -i '/^>/! s/~/-/g' "$ALN"
            fi
        else
            warn "[$FAM] could not fetch the HMM for $FAM from HMM_SOURCE; falling back to MAFFT"
            ALIGNER="mafft"
        fi
        ;;
    mafft) ;;
    *) fail "unknown ALIGNER=$ALIGNER (expected mafft or hmmalign)" ;;
esac

if [[ "$ALIGNER" == "mafft" ]]; then
    command -v mafft >/dev/null 2>&1 || fail "mafft not found"
    mafft_args=()
    case "$MAFFT_MODE" in
        linsi) mafft_args=(--localpair --maxiterate "$MAFFT_MAXITERATE") ;;
        einsi) mafft_args=(--genafpair --maxiterate "$MAFFT_MAXITERATE") ;;
        ginsi) mafft_args=(--globalpair --maxiterate "$MAFFT_MAXITERATE") ;;
        auto)
            if [[ "$N_USED" -le 200 ]]; then
                mafft_args=(--localpair --maxiterate "$MAFFT_MAXITERATE")
            else
                mafft_args=(--auto)
            fi
            ;;
        *) fail "unknown MAFFT_MODE=$MAFFT_MODE" ;;
    esac
    # shellcheck disable=SC2206
    extra=(${MAFFT_EXTRA:-})
    log "[$FAM] MAFFT ${mafft_args[*]} (n=$N_USED)"
    if ! mafft --thread 1 --anysymbol "${mafft_args[@]}" "${extra[@]}" "$ALN_IN" > "$ALN" \
            2> "${FAM_DIR}/mafft.log"; then
        fail "MAFFT failed"
    fi
fi

[[ -s "$ALN" ]] || fail "alignment is empty"

# ---------- trimAl ------------------------------------------------------------
command -v trimal >/dev/null 2>&1 || fail "trimal not found"
case "$TRIMAL_MODE" in
    gappyout|automated1|strict)
        trim_args=("-$TRIMAL_MODE")
        ;;
    *)
        # shellcheck disable=SC2206
        trim_args=($TRIMAL_MODE)
        ;;
esac
log "[$FAM] trimAl ${trim_args[*]}"
if ! trimal -in "$ALN" -out "$TRIM" -fasta "${trim_args[@]}" \
        > "${FAM_DIR}/trimal.log" 2>&1; then
    fail "trimAl failed"
fi
[[ -s "$TRIM" ]] || fail "trimAl output is empty"

ALN_COLS=$(awk 'BEGIN{m=0} /^>/{if(n>m)m=n; n=0; next} {n+=length($0)} END{if(n>m)m=n; print m}' "$ALN")
TRIM_COLS=$(awk 'BEGIN{m=0} /^>/{if(n>m)m=n; n=0; next} {n+=length($0)} END{if(n>m)m=n; print m}' "$TRIM")
N_TRIM=$(awk '/^>/{n++} END{print n+0}' "$TRIM")

log "[$FAM] alignment $ALN_COLS columns -> $TRIM_COLS columns after trimming, $N_TRIM sequences"
if [[ "$TRIM_COLS" -lt "$TRIMAL_MIN_COL" ]]; then
    skip "only $TRIM_COLS columns after trimming < TRIMAL_MIN_COL=$TRIMAL_MIN_COL"
fi
if [[ "$N_TRIM" -lt "$MIN_SEQS" ]]; then
    skip "only $N_TRIM sequences after trimming < MIN_SEQS=$MIN_SEQS"
fi
TRIM_FRAC=$(awk -v a="$TRIM_COLS" -v b="$ALN_COLS" 'BEGIN{printf "%.3f", a/b}')
if awk -v f="$TRIM_FRAC" -v t="$MIN_COL_FRACTION" 'BEGIN{exit !(f < t)}'; then
    warn "[$FAM] only ${TRIM_FRAC} of the alignment columns were retained (< $MIN_COL_FRACTION); the tree may be unstable"
fi

# ---------- IQ-TREE -----------------------------------------------------------
IQ_BIN="$IQTREE_BIN"
if ! command -v "$IQ_BIN" >/dev/null 2>&1; then
    if command -v iqtree2 >/dev/null 2>&1; then IQ_BIN="iqtree2"
    elif command -v iqtree >/dev/null 2>&1; then IQ_BIN="iqtree"
    else fail "iqtree2/iqtree not found"; fi
fi

iq_args=(-s "$TRIM" -st AA -m "$IQTREE_MODEL" -B "$UFBOOT" -bnni -alrt "$ALRT"
         -T "$THREADS_PER_JOB" -seed "$IQTREE_SEED" --prefix "$PREFIX" --quiet)
[[ -n "$IQTREE_MSET" ]] && iq_args+=(-mset "$IQTREE_MSET")
# shellcheck disable=SC2206
iq_extra=(${IQTREE_EXTRA:-})
iq_args+=("${iq_extra[@]}")

if [[ -n "$OUTGROUP" ]]; then
    if grep -Fqx ">${OUTGROUP}" "$TRIM"; then
        iq_args+=(-o "$OUTGROUP")
    else
        warn "[$FAM] OUTGROUP=$OUTGROUP is not present in the alignment; keeping the tree unrooted"
    fi
fi

log "[$FAM] IQ-TREE inference ($IQ_BIN, model=$IQTREE_MODEL, UFBoot=$UFBOOT)"
if ! "$IQ_BIN" "${iq_args[@]}" > "${FAM_DIR}/iqtree.stdout.log" 2>&1; then
    fail "IQ-TREE failed; see iqtree.stdout.log"
fi
[[ -s "${PREFIX}.treefile" ]] || fail "IQ-TREE did not produce a .treefile"

BEST_MODEL=$(grep -m1 'Best-fit model' "${PREFIX}.iqtree" 2>/dev/null \
    | sed 's/.*: //; s/[[:space:]].*//' | tr -d '\r' || true)
[[ -n "$BEST_MODEL" ]] || BEST_MODEL="$IQTREE_MODEL"

write_status "ok" "$N_IN" "$N_USED" "$ALN_COLS" "$TRIM_COLS" "$BEST_MODEL" \
    "aligner=$ALIGNER;trimal=$TRIMAL_MODE"
if [[ "$KEEP_INTERMEDIATE" -eq 0 ]]; then
    rm -f "$ALN" "$TRIM"
fi
log "[$FAM] done: $BEST_MODEL"
