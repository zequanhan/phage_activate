import sys
import os
import argparse
import re
from Bio import SeqIO
import subprocess

# 获取当前脚本的目录
current_dir = os.path.dirname(os.path.abspath(__file__))

# 将要导入的目录设为与脚本同一层的目录
sys.path.append(os.path.join(current_dir, 'total_step'))
from generate_result import *
from total_step_integrate_tfbs_and_promoter import GenomeAnalyzer

# 解析命令行参数
parser = argparse.ArgumentParser(description='Analyze genome using specified paths.')
parser.add_argument('-fasta_path', type=str, required=True, help='Path to the FASTA file')
parser.add_argument('-output_dir', type=str, required=True, help='Path to the output directory')
args = parser.parse_args()

# 获取 accession 名字
accession = os.path.basename(args.fasta_path).split('.')[0]

# 检查 output_dir 是否存在，如果不存在则创建
output_dir = os.path.join(args.output_dir, accession)
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# 指定自定义数据库路径
db_path = os.path.join(current_dir, 'blast_of_virus', 'prophage_virus.db')

# 使用 Prokka 注释 FASTA 文件，并将输出文件保存在指定目录
prokka_command = [
    'prokka',
    '--outdir', output_dir,
    '--prefix', accession,
    '--force',
    '--kingdom', 'Viruses',
    '--proteins', db_path,  # 使用自定义蛋白数据库
    args.fasta_path
]

# 运行 Prokka 并检查输出
try:
    env = os.environ.copy()
    env['PATH'] = '/home/hanzequan/miniconda3/envs/DPProm/bin:' + env['PATH']
    subprocess.run(prokka_command, check=True, env=env)
except subprocess.CalledProcessError as e:
    print(f"Error: Prokka command failed with return code {e.returncode}")
    sys.exit(1)

# 获取 Prokka 生成的 GBK 文件路径
gbk_file = os.path.join(output_dir, f"{accession}.gbk")
if not os.path.exists(gbk_file):
    print(f"Error: Prokka did not generate the expected GBK file at {gbk_file}")
    sys.exit(1)
print('GBK file generated:', gbk_file)

# 修正 GBK 文件中的日期格式
corrected_gbk_file = os.path.join(output_dir, f"{accession}_corrected.gbk")
with open(gbk_file, "r", encoding="utf-8") as infile, open(corrected_gbk_file, "w", encoding="utf-8") as outfile:
    for line in infile:
        if line.startswith("LOCUS"):
            line = re.sub(r'(\d{2})-(\d{1,2})月-(\d{4})', r'\1-AUG-\3', line)
        outfile.write(line)

print('Corrected GBK file saved:', corrected_gbk_file)

# 实例化 GenomeAnalyzer 并进行分析
analyzer = GenomeAnalyzer(corrected_gbk_file, args.fasta_path, output_dir)
meme_results, df_promoters, pwm_df = analyzer.analyze_genome()

# 检查 meme_results 是否为空
if meme_results is None:
    print("Error: No results found from genome analysis. Please check your input files and Prokka output.")
    sys.exit(1)

# 调用 analyze_tfbs_modification 进行进一步分析
analyze_tfbs_modification(meme_results, df_promoters, pwm_df, genbank_path=corrected_gbk_file, output_path=output_dir)

print("Genome analysis completed successfully. Results are saved in:", output_dir)

