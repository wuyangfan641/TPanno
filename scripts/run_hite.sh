#!/bin/bash

# 检查输入参数
if [ $# -ne 1 ]; then
    echo "用法: $0 <输入目录>"
    exit 1
fi

input_dir="$1"

# 验证输入目录是否存在
if [ ! -d "$input_dir" ]; then
    echo "错误：输入目录 $input_dir 不存在。"
    exit 1
fi

# 创建主输出目录（如果不存在）
mkdir -p "../hite_output"

# 遍历输入目录中的所有文件
for input_file in "$input_dir"/*; do
    if [ -f "$input_file" ]; then
        echo "正在处理文件: $input_file"
        
        # 提取文件名前缀（去除扩展名）
        filename=$(basename -- "$input_file")
        prefix="${filename%.*}"
        
        # 为每个文件创建独立的输出子目录
        output_subdir="../hite_output/$prefix"
        mkdir -p "$output_subdir"
        
        # 运行HiTE分析文件
        python main.py --genome "$input_file" --outdir "$output_subdir" --thread 12 --annotate 1 --plant 0
        
        # 重命名输出文件
        output_file="$output_subdir/genome.rename.fa"
        if [ -f "$output_file" ]; then
            mv "$output_file" "$output_subdir/${prefix}.fa.masked"
            echo "已重命名：$output_subdir/${prefix}.fa.masked"
        else
            echo "警告：未找到输出文件 $output_file"
        fi
    fi
done

echo "处理完成！"
