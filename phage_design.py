import sys
import os
import argparse
import tempfile
from Bio import SeqIO

current_dir = os.path.dirname(os.path.abspath(__file__))

# 将要导入的目录设为与脚本同一层的目录
sys.path.append(os.path.join(current_dir, 'total_step'))
from generate_result import *
from total_step_integrate_tfbs_and_promoter import GenomeAnalyzer

parser = argparse.ArgumentParser(description='Analyze genome using specified paths.')
parser.add_argument('-gbk_path', type=str, required=True, help='Path to the gbk file')
parser.add_argument('-output_dir', type=str, required=True, help='Path to the output directory')

args = parser.parse_args()

# 获取脚本当前所在的目录
current_dir = os.path.dirname(os.path.abspath(__file__))

# 获取 accession 名字
accession = os.path.basename(args.gbk_path).split('.')[0]

# 检查 output_dir 是否存在，如果不存在则创建
output_dir = os.path.join(args.output_dir, accession)
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# 生成临时FASTA文件
with tempfile.NamedTemporaryFile(delete=False, suffix=".fasta") as temp_fasta:
    fasta_path = temp_fasta.name
    with open(fasta_path, "w") as fasta_file:
        for record in SeqIO.parse(args.gbk_path, "genbank"):
            SeqIO.write(record, fasta_file, "fasta")

# 实例化 GenomeAnalyzer 并进行分析
analyzer = GenomeAnalyzer(args.gbk_path, fasta_path, output_dir)
meme_results, df_promoters, pwm_df = analyzer.analyze_genome()

# 调用 analyze_tfbs_modification 进行进一步分析
analyze_tfbs_modification(meme_results, df_promoters, pwm_df, genbank_path=args.gbk_path, output_path=output_dir)

# 删除临时FASTA文件
os.remove(fasta_path)

