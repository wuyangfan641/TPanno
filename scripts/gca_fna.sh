#!/bin/bash

# 脚本: gca.sh
# 作用: 通过GCA编号批量下载基因组文件
# 用法: ./gca.sh <输入文件>
# 输入文件格式: 每行包含物种名和GCA编号，用制表符分隔 (例如 vv.tsv)

# 检查输入参数
if [[ $# -ne 1 ]]; then
    echo "用法: $0 <输入文件>"
    echo "输入文件格式：每行包含物种名和GCA编号，用制表符分隔"
    exit 1
fi

input_file="$1"

# 检查输入文件是否存在
if [[ ! -f "$input_file" ]]; then
    echo "错误：输入文件 '$input_file' 不存在"
    exit 1
fi

# 创建输出目录
output_dir="genomes"
mkdir -p "$output_dir"

# 读取并处理每一行
while IFS=$'\t' read -r species gca || [[ -n "$species" ]]; do
    # 跳过空行和注释行
    [[ -z "$species" || "$species" == \#* ]] && continue

    # 清理物种名和GCA编号中的前后空格和回车符
    species=$(echo "$species" | sed 's/^[ \t]*//;s/[ \t\r]*$//')
    gca=$(echo "$gca" | sed 's/^[ \t]*//;s/[ \t\r]*$//')

    # 如果清理后GCA为空，则跳过
    if [[ -z "$gca" ]]; then
        continue
    fi

    # 将物种名转换为安全的文件名 (替换空格和特殊字符)
    safe_species_name=$(echo "$species" | tr ' ' '_' | tr -cd 'a-zA-Z0-9_.-')

    echo "正在处理: $species ($gca)"

    output_fna="$output_dir/${safe_species_name}.fna"

    # 检查目标文件是否已存在
    if [ -f "$output_fna" ]; then
        echo "  文件 '$output_fna' 已存在，跳过。"
        continue
    fi

    # 使用NCBI Datasets API v2构建下载链接
    download_url="https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/${gca}/download?include_annotation_type=GENOME_FASTA"
    
    # 创建临时目录用于下载和解压
    temp_dir=$(mktemp -d)
    zip_file="${temp_dir}/${safe_species_name}.zip"

    echo "  正在下载蛋白质组..."
    if wget --header="Accept: application/zip" -qO "$zip_file" "$download_url"; then
        echo "  下载完成，正在解压..."
        unzip -q "$zip_file" -d "$temp_dir"

        # 在解压文件中查找 .fna 文件
        fna_file=$(find "$temp_dir/ncbi_dataset/data" -type f -name "*.fna" 2>/dev/null)

        if [[ -f "$fna_file" ]]; then
            # 移动并重命名文件
            mv "$fna_file" "$output_fna"
            echo "  成功: 基因组已保存为 '$output_fna'"
        else
            echo "  错误: 在下载的压缩包中未找到 .faa 文件。"
        fi
    else
        echo "  错误: 下载失败 '$species' ($gca)."
    fi

    # 清理临时文件
    rm -rf "$temp_dir"

done < "$input_file"

echo "所有任务已完成。" 
