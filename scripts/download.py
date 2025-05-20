#!/usr/bin/env python3
import ftplib
import re
import os
import sys
import subprocess
from urllib.parse import urlparse
from urllib.request import urlretrieve

def acc2path(acc):
    """将GCA编号转换为NCBI FTP基础路径"""
    parts = acc.split('_')
    if len(parts) < 2:
        raise ValueError(f"无效的GCA编号: {acc}")
    acc_part = parts[1].split('.')[0]
    return f"/genomes/all/GCA/{acc_part[0:3]}/{acc_part[3:6]}/{acc_part[6:9]}"

def extract_version(subdir):
    """从子目录名提取版本号"""
    match = re.search(r'\.(\d+)(?:_|$)', subdir)
    return int(match.group(1)) if match else 0

def get_latest_fna_url(ftp, gca):
    """获取指定GCA编号的最新版.fna.gz文件URL"""
    base_path = acc2path(gca)
    try:
        ftp.cwd(base_path)
    except ftplib.error_perm:
        print(f"错误: 路径不存在 {base_path} (GCA: {gca})")
        return None

    # 获取子目录列表并选择最新版本
    dirs = []
    ftp.retrlines('LIST', dirs.append)
    subdirs = [line.split()[-1] for line in dirs if line.startswith('d')]
    if not subdirs:
        print(f"错误: {base_path} 下无子目录 (GCA: {gca})")
        return None

    # 按版本号排序
    version_dirs = [(extract_version(d), d) for d in subdirs]
    version_dirs.sort(reverse=True, key=lambda x: x[0])
    latest_dir = version_dirs[0][1]
    full_path = f"{base_path}/{latest_dir}"

    # 获取文件列表
    ftp.cwd(latest_dir)
    files = []
    ftp.retrlines('LIST', files.append)
    filenames = [line.split()[-1] for line in files if not line.startswith('d')]

    # 筛选符合条件的.fna.gz文件
    for fname in filenames:
        if fname.endswith('.fna.gz') and '_cds_' not in fname and '_rna_' not in fname:
            return f"ftp://ftp.ncbi.nlm.nih.gov{full_path}/{fname}"
    return None


def download_file(url, species):
    """使用wget下载文件并重命名"""
    # 清理物种名中的特殊字符
    species_clean = species.replace(' ', '_').replace('/', '_').strip()
    parsed = urlparse(url)
    filename = os.path.basename(parsed.path)
    ext = '.fna.gz' if filename.endswith('.fna.gz') else ''
    output_name = f"{species_clean}{ext}"

    # 构造wget命令
    wget_cmd = [
        'wget',
        '-c',
        '--show-progress',  # 显示进度条
        '-O', output_name,  # 指定输出文件名
        url
    ]

    try:
        print(f"下载 {species} -> {output_name}")
        result = subprocess.run(
            wget_cmd,
            check=True,  # 检查命令执行状态
            stderr=subprocess.STDOUT  # 捕获错误输出
        )
        print(f"成功: {output_name}")
        return True
    except subprocess.CalledProcessError as e:
        print(f"下载失败 (状态码 {e.returncode}): {e.output.decode()}")
        if os.path.exists(output_name):
            os.remove(output_name)
        return False
    except Exception as e:
        print(f"未知错误: {str(e)}")
        if os.path.exists(output_name):
            os.remove(output_name)
        return False

def main():
    if len(sys.argv) != 2:
        print("用法: ./genome.py <genome.tsv>")
        sys.exit(1)

    tsv_file = sys.argv[1]
    if not os.path.exists(tsv_file):
        print(f"错误: 文件 {tsv_file} 不存在")
        sys.exit(1)

    # 读取genome.tsv文件
    gcadata = []
    with open(tsv_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) < 2:
                print(f"警告: 跳过无效行: {line}")
                continue
            gcadata.append((parts[0], parts[1]))  # (GCA, species)

    # 连接FTP
    ftp = ftplib.FTP('ftp.ncbi.nlm.nih.gov')
    ftp.login()

    # 获取所有GCA的下载链接
    download_list = []
    for gca, species in gcadata:
        url = get_latest_fna_url(ftp, gca)
        if url:
            download_list.append((species, url))
        else:
            print(f"警告: 未找到 {gca} 的FTP路径")

    ftp.quit()

    # 下载文件
    success = 0
    for species, url in download_list:
        if download_file(url, species):
            success += 1

    print(f"完成! 成功下载 {success}/{len(download_list)} 个文件")

if __name__ == '__main__':
    main()
