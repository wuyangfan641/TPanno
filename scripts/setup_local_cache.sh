#!/usr/bin/env bash
# =============================================================================
#  setup_local_cache.sh -- populate the EGAPx local_cache with rsync
#
#  Run from the repository root:
#      bash scripts/setup_local_cache.sh
#
#  The cache is written to ./local_cache and can be passed to EGAPx with
#  -lc ./local_cache.
# =============================================================================
set -euo pipefail

command -v rsync >/dev/null 2>&1 || { echo "ERROR: rsync is required" >&2; exit 1; }

mkdir -p local_cache
cd local_cache || { echo "ERROR: cannot enter local_cache" >&2; exit 1; }

mkdir -p gnomon ortholog_references target_proteins taxonomy \
         reference_sets misc cmsearch busco_downloads/lineages

echo "Syncing gnomon..."
rsync --copy-links --recursive --times --verbose \
    rsync://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/EGAP/support_data/gnomon/2 gnomon/

echo "Syncing ortholog_references..."
rsync --copy-links --recursive --times --verbose \
    rsync://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/EGAP/support_data/ortholog_references/3 ortholog_references/

echo "Syncing target_proteins..."
rsync --copy-links --recursive --times --verbose \
    rsync://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/EGAP/support_data/target_proteins/2 target_proteins/

echo "Syncing taxonomy..."
rsync --copy-links --recursive --times --verbose \
    rsync://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/EGAP/support_data/taxonomy/2 taxonomy/

echo "Syncing reference_sets..."
rsync --copy-links --recursive --times --verbose \
    rsync://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/EGAP/support_data/reference_sets/3 reference_sets/

echo "Syncing misc..."
rsync --copy-links --recursive --times --verbose \
    rsync://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/EGAP/support_data/misc/2 misc/

echo "Syncing cmsearch..."
rsync --copy-links --recursive --times --verbose \
    rsync://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/EGAP/support_data/cmsearch/1 cmsearch/

echo "All local_cache data synchronized."
