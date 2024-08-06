import subprocess
import tempfile
import pandas as pd
import numpy as np
import logomaker
import matplotlib.pyplot as plt
from Bio.Seq import Seq
import os
import pandas as pd
from Bio import motifs
import re
from Bio import SeqIO
from Bio.Blast import NCBIXML
import os
import subprocess
import subprocess
import tempfile
import pandas as pd
import logomaker
import matplotlib.pyplot as plt
from Bio.Seq import Seq
import os
import pandas as pd
from Bio import motifs
import re


def find_directories_with_string(root_directory, target_string, file_extension=".txt"):
    """
    遍历指定目录及其所有子目录，查找包含特定字符串的文件，并返回这些文件所在的目录。
    """
    directories = set()
    for root, dirs, files in os.walk(root_directory):
        for file in files:
            if file.endswith(file_extension):
                file_path = os.path.join(root, file)
                with open(file_path, 'r', encoding='utf-8') as f:
                    if target_string in f.read():
                        directories.add(root)
    return directories


def build_motif_matrices(directory, sequence_count_occurrences):
    """
    在指定目录及其所有子目录中查找所有的'meme.xml'文件，并构建motif矩阵。
    返回该目录及其所有子目录中找到的所有motif数据的DataFrame。
    """
    all_motifs_data = pd.DataFrame()
    found_xml = False  # 标记是否找到至少一个meme.xml文件

    for root, dirs, files in os.walk(directory):
        for file in files:
            if file == 'meme.xml':
                found_xml = True
                xml_file = os.path.join(root, file)
                txt_file = os.path.join(os.path.dirname(root), os.path.basename(root) + '.txt')

                try:
                    with open(txt_file, 'r') as txt:
                        for line in txt.readlines():
                            line = line.strip()
                            if line:
                                sequence_ids = line.split(',')
                                sequence_count = len(sequence_ids)
                                if sequence_count in sequence_count_occurrences:
                                    sequence_count_occurrences[sequence_count] += 1
                                else:
                                    sequence_count_occurrences[sequence_count] = 0
                except FileNotFoundError:
                    print(f"Warning: Corresponding text file not found: {txt_file}")

                with open(xml_file) as f:
                    meme_record = motifs.parse(f, "meme")

                motif_index = 1
                for motif in meme_record:
                    sequences = [Seq(str(instance)) for instance in motif.instances]
                    motifs_data = []
                    for instance in motif.instances:
                        sequence_name = instance.sequence_name
                        id = motif.id
                        consensus = motif.name
                        e_value = motif.evalue
                        num_occurrences = len(motif.instances)
                        suffix = f".{sequence_count_occurrences[sequence_count]}" if sequence_count_occurrences[
                                                                                         sequence_count] > 0 else ""
                        motif_data = {
                            "Number": f"{sequence_count}{suffix}",
                            "Layer": id,
                            "Strand": instance.strand,
                            "Start": instance.start,
                            "p-value": instance.pvalue,
                            "e-value": e_value,
                            "Sequence": str(instance),
                            "Motif": consensus
                        }
                        motifs_data.append(motif_data)
                    if motifs_data:
                        motifs_df = pd.DataFrame(motifs_data)
                        all_motifs_data = pd.concat([all_motifs_data, motifs_df], ignore_index=True)
                    motif_index += 1

    if not found_xml:
        raise FileNotFoundError(f"No MEME output file found in directory {directory} or its subdirectories.")

    return all_motifs_data


#使用needle算法合并序列并找到最关键的motif

def run_needle(seq1, seq2):  # 使用needle算法合并序列
    seq1_rc = str(Seq(seq1).reverse_complement())
    results = {}   
    sequences_to_test = [(seq1, "original"), (seq1_rc, "reverse complement")]

    # 获取脚本当前所在的目录
    current_dir = os.path.dirname(os.path.abspath(__file__))
    # 构建 needle 的相对路径
    needle_path = os.path.join(current_dir, '..', 'usr', 'bin', 'needle')

    for sequence, label in sequences_to_test:
        with tempfile.NamedTemporaryFile('w+', delete=False) as temp_seq1, \
                tempfile.NamedTemporaryFile('w+', delete=False) as temp_seq2:
            temp_seq1.write(f">seq1\n{sequence}\n")
            temp_seq2.write(f">seq2\n{seq2}\n")
            temp_seq1.flush()
            temp_seq2.flush()
            
            output_path = tempfile.NamedTemporaryFile('w+', delete=False).name
            command = [ 
                needle_path,  # 使用相对路径
                '-asequence', temp_seq1.name,
                '-bsequence', temp_seq2.name,
                '-gapopen', '10',
                '-gapextend', '0.5',
                '-outfile', output_path,
                '-auto', 'yes'
            ]

            # 打印命令以进行调试
            print("Running command: ", " ".join(command))

            # 执行命令
            os.system(" ".join(command))

            # 读取结果文件内容并存储在results字典中
            with open(output_path, 'r') as result_file:
                results[label] = result_file.read()

            try:
                subprocess.run(command, check=True)
                details = extract_alignment_details(output_path)
                if details:
                    results[label] = details
            except subprocess.CalledProcessError as e:
                print(f"Error running needle with {label}: {e}")
            except FileNotFoundError as e:
                print(f"Error: Output file not found with {label}: {e}")

    if results:
        max_label = max(results, key=lambda x: results[x].get('Score', float('-inf')))
        return pd.DataFrame([results[max_label]])
    else:
        return pd.DataFrame()


def extract_alignment_details(file_path):
    details = {}
    try:
        with open(file_path, 'r') as file:
            content = file.read()
            length_match = re.search(r'# Length:\s*(\d+)', content)
            identity_match = re.search(r'# Identity:\s*(\d+)/(\d+)\s*\(\s*(\d+\.\d+)%\)', content)
            similarity_match = re.search(r'# Similarity:\s*(\d+)/(\d+)\s*\(\s*(\d+\.\d+)%\)', content)
            gaps_match = re.search(r'# Gaps:\s*(\d+)/(\d+)\s*\(\s*(\d+\.\d+)%\)', content)
            score_match = re.search(r'# Score:\s*(\d+\.\d+)', content)

            if length_match:
                details['Length'] = int(length_match.group(1))
            if identity_match:
                details['Identity'] = identity_match.groups()  # Save as tuple
            if similarity_match:
                details['Similarity'] = similarity_match.groups()  # Save as tuple
            if gaps_match:
                details['Gaps'] = gaps_match.groups()  # Save as tuple
            if score_match:
                details['Score'] = float(score_match.group(1))
    except Exception as e:
        print(f"Failed to read or parse the results file: {e}")
    return details


def run_comparisons_on_motifs(df):
    first_motifs = {}
    unique_numbers = df['Number'].unique()
    # 收集每个 Number 的每个 Layer 的第一个 Motif
    for number in unique_numbers:
        for layer in ['motif_1', 'motif_2', 'motif_3']:
            if not df[(df['Number'] == number) & (df['Layer'] == layer)].empty:
                motif = df[(df['Number'] == number) & (df['Layer'] == layer)]['Motif'].iloc[0]
                first_motifs[motif] = 0  # 使用序列作为键，初始化状态为0
    results_list = []

    # 将字典的键（序列）和其状态用于比对
    for motif1 in first_motifs.keys():
        for motif2 in first_motifs.keys():
            if motif1 != motif2:  # 确保不与自身比对
                result_df = run_needle(motif1, motif2)

                if not result_df.empty:
                    max_identity = float(result_df['Identity'][0][-1])
                    if max_identity > 70:
                        # 如果 Identity 大于 70，增加对应序列的状态
                        first_motifs[motif1] += 1
                    # 收集比对结果
                    for _, row in result_df.iterrows():
                        row['Original Motif'] = motif1
                        row['Target Motif'] = motif2
                        row['State Change'] = first_motifs[motif1]
                        results_list.append(row.to_dict())
    return pd.DataFrame(results_list)


def merge_sequences_based_on_identity(final_df):
    # 从 final_df 中提取序列和其状态
    motif_counts = final_df['Original Motif'].value_counts().to_dict()
    motifs = list(motif_counts.keys())

    # 初始化序列状态
    motif_states = {motif: count for motif, count in motif_counts.items()}

    # 相互比对序列
    for i, motif1 in enumerate(motifs):
        for motif2 in motifs[i + 1:]:
            if motif1 != motif2 and motif1 in motif_states and motif2 in motif_states:
                result_df = run_needle(motif1, motif2)
                if not result_df.empty:
                    max_identity = float(result_df['Identity'][0][-1])
                    if max_identity > 70:  #最大相似度
                        # 确定哪个序列具有更高的状态，进行合并
                        larger_motif = motif1 if motif_states[motif1] >= motif_states[motif2] else motif2
                        smaller_motif = motif2 if motif1 == larger_motif else motif1

                        # 更新状态并合并
                        motif_states[larger_motif] += motif_states.pop(smaller_motif, 0)

    # 创建最终的 DataFrame 显示序列和其最终状态
    final_results = [{'Motif': motif, 'Final State': state} for motif, state in motif_states.items()]
    return pd.DataFrame(final_results)


# 假设 run_needle 函数返回的是一个包含比对结果详情的字典

def create_motif_visualization(accession):
    root_directory = '/home/hanzequan/test_bectiral/operator_recongize/all_tree'
    directories = find_directories_with_string(root_directory, accession)
    # 汇总所有目录中的motif数据
    all_motifs_df = pd.DataFrame()
    sequence_count_occurrences = {}  # 初始化计数字典
    for directory in directories:
        try:
            meme_df = build_motif_matrices(directory, sequence_count_occurrences)  # 传递计数字典
            all_motifs_df = pd.concat([all_motifs_df, meme_df], ignore_index=True)
        except FileNotFoundError as e:
            print(e)

    comparison_results_df = run_comparisons_on_motifs(all_motifs_df)
    final_df = comparison_results_df[
        comparison_results_df['State Change'] == comparison_results_df['State Change'].max()]
    best_motif = all_motifs_df[all_motifs_df['Motif'] == merge_sequences_based_on_identity(final_df)[
        merge_sequences_based_on_identity(final_df).index == 0]['Motif'][0]]
    sequences = best_motif['Sequence'].to_list()

    # 寻找最短的序列长度
    min_length = min(len(seq) for seq in sequences)

    # 截断所有序列到最短长度
    sequences_truncated = [seq[:min_length] for seq in sequences]
    sequences_truncated = [Seq(seq) for seq in sequences_truncated]

    # 使用BioPython创建motif
    m = motifs.create(sequences_truncated)
    pwm = m.counts.normalize(pseudocounts={'A': 0, 'C': 0, 'G': 0, 'T': 0})
    pwm_df = pd.DataFrame({nucleotide: pwm[nucleotide] for nucleotide in 'ACGT'})
    ic = logomaker.transform_matrix(pwm_df, from_type='probability', to_type='information')

    # 生成并展示序列标志
    fig, ax = plt.subplots(figsize=(10, 3))
    logo = logomaker.Logo(ic, ax=ax)
    ax.set_ylabel('bits', fontsize=12)
    ax.set_ylim(bottom=0)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.yaxis.tick_left()
    plt.show()
def extract_sequence(gbk_path, protein_names):
    """
    从GenBank文件中提取第一个特定蛋白质名称的序列，并在找到第一个匹配后停止搜索。
    """
    record = SeqIO.read(gbk_path, "genbank")
    for feature in record.features:
        if feature.type == "CDS" and "product" in feature.qualifiers:
            product_name = feature.qualifiers["product"][0].lower()
            if any(protein_name.lower() in product_name for protein_name in protein_names):
                print('#'*30)
                print('find repressor protein_name:',protein_name)
                return feature.qualifiers["translation"][0]  # 返回第一个匹配的序列
    return None  # 如果没有找到匹配项，则返回None


def write_sequences_to_fasta(sequences, output_path):
    """
    将序列写入FASTA文件。
    """
    with open(output_path, "w") as file:
        for idx, seq in enumerate(sequences):
            file.write(f">protein_{idx}\n{seq}\n")

def run_blastp(query_fasta, db_path, output_path):
    """
    使用blastp进行蛋白质比对。
    """
    command = [
        'blastp',
        '-query', query_fasta,
        '-db', db_path,
        '-out', output_path,
        '-outfmt', "6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore"
    ]
    subprocess.run(command, check=True)


def parse_blast_results(blast_output):
    """
    解析BLAST比对结果，并返回特定信息。
    """
    column_names = ["query_accver", "subject_accver", "percentage_identity", "alignment_length", "mismatches",
                    "gap_opens", "query_start", "query_end", "subject_start", "subject_end", "evalue", "bit_score"]
    try:
        df = pd.read_csv(blast_output, sep="\t", names=column_names)
        if not df.empty:
            # 取得最佳匹配（第一行）
            best_hit = df.iloc[0]
            # 格式化subject_accver以去除'.'后的内容
            formatted_acc = best_hit['subject_accver'].split('.')[0]
            # 提取蛋白名和e值
            protein_name = best_hit['subject_accver']  # 假设subject_accver包含了蛋白名称
            evalue = best_hit['evalue']
            return {
                'Accession': formatted_acc,
                'E-value': evalue,
                'Protein Name': protein_name
            }
        else:
            return None
    except Exception as e:
        print(f"Error reading BLAST output: {e}")
        return None

def main(gbk_path):
    protein_names = [
        'repressor', 'transcriptional regulator', 'immunity repressor', 
        'transcriptional repressor', 'Cro/CI family transcriptional regulator', 
        'Hxr', 'CI protein', 'CII-like transcriptional activator'
    ]
    query_fasta = "/home/public_new/dowlond_phage/query_proteins.fasta"
    db_path = "/home/public_new/dowlond_phage/all_phage_tree/blast_db"
    blast_output = "/home/public_new/dowlond_phage/blast_results.tsv"
    
    # 解析GenBank文件，提取序列
    sequence = extract_sequence(gbk_path, protein_names)
    if sequence:
        write_sequences_to_fasta([sequence], query_fasta)  # 包装序列为列表
        # 运行BLAST比对
        run_blastp(query_fasta, db_path, blast_output)
        # 分析BLAST结果
        result = parse_blast_results(blast_output)
        if result:
            accession = result['Accession']
            print(f"Accession: {accession}, E-value: {result['E-value']}, Protein Name: {result['Protein Name']}")
            create_motif_visualization(accession)  # 调用函数构建图片并传递 gbk_path
        else:
            print("No BLAST hits found or failed to parse BLAST output.")

# 示例函数调用
#if __name__ == "__main__":
#    gbk_path = "/home/public_new/dowlond_phage/phage_gbk/NC_001884.gbk"
#    main(gbk_path)


