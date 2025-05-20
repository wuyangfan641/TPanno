#!/bin/bash
set -euo pipefail

# 创建主目录并进入
mkdir -p local_cache
cd local_cache || { echo "无法进入local_cache目录"; exit 1; }

# 创建子目录
mkdir -p gnomon ortholog_references target_proteins taxonomy reference_sets misc

# 同步gnomon数据
echo "同步gnomon数据..."
rsync --copy-links --recursive --times --verbose rsync://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/EGAP/support_data/gnomon/2 gnomon/

# 同步ortholog_references数据
echo "同步ortholog_references数据..."
rsync --copy-links --recursive --times --verbose rsync://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/EGAP/support_data/ortholog_references/2 ortholog_references/

# 同步target_proteins数据
echo "同步target_proteins数据..."
rsync --copy-links --recursive --times --verbose rsync://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/EGAP/support_data/target_proteins/2 target_proteins/

# 同步taxonomy数据
echo "同步taxonomy数据..."
rsync --copy-links --recursive --times --verbose rsync://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/EGAP/support_data/taxonomy/1 taxonomy/

# 同步reference_sets数据
echo "同步reference_sets数据..."
rsync --copy-links --recursive --times --verbose rsync://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/EGAP/support_data/reference_sets/2 reference_sets/

# 同步misc数据
echo "同步misc数据..."
rsync --copy-links --recursive --times --verbose rsync://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/EGAP/support_data/misc/2 misc/

echo "所有数据同步完成！"
