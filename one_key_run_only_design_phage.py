import sys
import os

# 获取脚本当前所在的目录
current_dir = os.path.dirname(os.path.abspath(__file__))

# 将要导入的目录设为与脚本同一层的目录
sys.path.append(os.path.join(current_dir, 'total_step'))

# 导入相关包
from generate_result import *
from total_step_integrate_tfbs_and_promoter import GenomeAnalyzer



accession = 'KX897981'
gbk_path = "/home/hanzequan/saher_file/test_data/"
fasta_path = '/home/hanzequan/saher_file/test_data/'
output_dir = f'/home/hanzequan/saher_file/test_data/{accession}'

# 添加 accession 并添加相应的扩展名
gbk_file = f"{gbk_path}{accession}.gbk"
fasta_file = f"{fasta_path}{accession}.fasta"
# output_path = f"{output_dir}{accession}"
analyzer = GenomeAnalyzer(gbk_file, fasta_file, output_dir)
meme_results, df_promoters, pwm_df=analyzer.analyze_genome()
analyze_tfbs_modification(meme_results, df_promoters, pwm_df,genbank_path=gbk_file,output_path=output_dir)
