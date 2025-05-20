#!/bin/bash

# 检查是否输入了SRR编号
if [ $# -ne 1 ]; then
    echo "错误：请提供一个SRA编号作为参数"
    echo "用法：$0 SRRxxxx"
    exit 1
fi

SRR="$1"

# 创建输出目录
mkdir -p sradir

# 第一步：使用prefetch下载SRA数据
echo "正在下载 $SRR..."
prefetch "$SRR"

# 检查prefetch是否成功
if [ $? -ne 0 ]; then
    echo "错误：SRA数据下载失败"
    exit 1
fi

# 第二步：使用fasterq-dump转换格式
echo "正在转换 $SRR 为FASTA..."
fasterq-dump --skip-technical \
             --threads 6 \
             --split-files \
             --seq-defline ">\$ac.\$si.\$ri" \
             --fasta \
             -O sradir/ \
             "./$SRR"

# 检查fasterq-dump是否成功
if [ $? -ne 0 ]; then
    echo "错误：FASTA转换失败"
    exit 1
fi

echo "操作成功完成！结果保存在 sradir/ 目录"