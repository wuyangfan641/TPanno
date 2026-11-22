#!/usr/bin/env bash
# =============================================================================
#  gca_faa.sh -- batch-download protein FASTA files by GCA accession
#
#  Usage: bash scripts/gca_faa.sh <input.tsv>
#
#  Input format: one record per line, tab-separated:
#      <prefix>    <GCA_accession>
#  A header row is optional and is skipped if present.
#
#  Output: proteins/<prefix>.faa
# =============================================================================
set -uo pipefail

if [[ $# -ne 1 ]]; then
    echo "Usage: $0 <input.tsv>" >&2
    echo "Input format: <prefix><TAB><GCA_accession>, one record per line" >&2
    exit 1
fi

input_file="$1"
if [[ ! -f "$input_file" ]]; then
    echo "ERROR: input file not found: $input_file" >&2
    exit 1
fi

for cmd in wget unzip; do
    command -v "$cmd" >/dev/null 2>&1 || { echo "ERROR: $cmd is required" >&2; exit 1; }
done

output_dir="proteins"
mkdir -p "$output_dir"

n_ok=0
n_skip=0
n_fail=0
line_no=0

while IFS=$'\t' read -r prefix accession _rest || [[ -n "${prefix:-}" ]]; do
    line_no=$((line_no + 1))
    prefix="${prefix%$'\r'}"
    accession="${accession%$'\r'}"
    # Strip a UTF-8 BOM from the first line (Windows editors may add one).
    if [[ "$line_no" -eq 1 ]]; then
        prefix="${prefix#$'\xEF\xBB\xBF'}"
    fi
    prefix="$(echo "$prefix" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')"
    accession="$(echo "$accession" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')"

    # Skip blank lines and comments
    [[ -z "$prefix" ]] && continue
    [[ "$prefix" == \#* ]] && continue

    # Skip an optional header row (but not a sample literally named "sample")
    if [[ "${prefix,,}" =~ ^(prefix|species|sample)$ ]] && \
       [[ "${accession,,}" =~ ^(assembly|accession|gca|fasta)$ ]]; then
        continue
    fi

    if [[ ! "$accession" =~ ^GC[AF]_[0-9]+\.[0-9]+$ ]]; then
        echo "WARNING: line $line_no: invalid GCA accession '$accession'; skipping '$prefix'" >&2
        n_fail=$((n_fail + 1))
        continue
    fi

    safe_name="$(echo "$prefix" | tr ' ' '_' | tr -cd 'a-zA-Z0-9_.-')"
    if [[ -z "$safe_name" ]]; then
        echo "WARNING: line $line_no: prefix '$prefix' has no usable characters; skipping" >&2
        n_fail=$((n_fail + 1))
        continue
    fi

    output_faa="${output_dir}/${safe_name}.faa"
    if [[ -f "$output_faa" ]]; then
        echo "  [skip] $output_faa already exists"
        n_skip=$((n_skip + 1))
        continue
    fi

    echo "  [get] $prefix ($accession) -> $output_faa"
    download_url="https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/${accession}/download?include_annotation_type=PROT_FASTA"

    temp_dir="$(mktemp -d)"
    zip_file="${temp_dir}/${safe_name}.zip"
    ok=0

    if wget --tries=3 --timeout=60 --header="Accept: application/zip" -qO "$zip_file" "$download_url"; then
        if unzip -q "$zip_file" -d "$temp_dir"; then
            faa_file="$(find "$temp_dir/ncbi_dataset/data" -type f -name "*.faa" -print -quit 2>/dev/null)"
            if [[ -n "$faa_file" && -f "$faa_file" ]]; then
                mv "$faa_file" "$output_faa"
                echo "         saved protein FASTA"
                ok=1
            else
                echo "ERROR: no .faa file found in the downloaded archive for $prefix ($accession)" >&2
            fi
        else
            echo "ERROR: could not unzip the archive for $prefix ($accession)" >&2
        fi
    else
        echo "ERROR: download failed for $prefix ($accession)" >&2
    fi

    rm -rf "$temp_dir"
    if [[ "$ok" -eq 1 ]]; then
        n_ok=$((n_ok + 1))
    else
        n_fail=$((n_fail + 1))
        rm -f "$output_faa"
    fi
done < "$input_file"

echo "Done: $n_ok downloaded, $n_skip skipped, $n_fail failed."
[[ "$n_fail" -eq 0 ]]
