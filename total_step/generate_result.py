import subprocess
from Bio import motifs
from Bio.Seq import Seq
import pandas as pd
import os
import matplotlib.pyplot as plt
import logomaker
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
def get_pwms():
    pwm_35 = {
        'A': [0.0000, 0.0784, 0.0362, 0.4894, 0.3605, 0.4208],
        'C': [0.1109, 0.0656, 0.0747, 0.2851, 0.3605, 0.0769],
        'G': [0.1267, 0.0181, 0.6192, 0.1041, 0.0000, 0.2225],
        'T': [0.7624, 0.8379, 0.2700, 0.1214, 0.2790, 0.2798]
    }
    pwm_10 = {
        'A': [0.0097, 1.0000, 0.3363, 0.5335, 0.4963, 0.0781],
        'C': [0.0618, 0.0000, 0.1190, 0.1094, 0.2299, 0.0268],
        'G': [0.1042, 0.0000, 0.0856, 0.1317, 0.1399, 0.0000],
        'T': [0.8244, 0.0000, 0.4591, 0.2254, 0.1339, 0.8951]
    }
    return pwm_35, pwm_10

def round_score(score):
    return round(score, 2)
def calculate_pwm_score(sequence, pwm):
    score = 0
    for i, base in enumerate(sequence):
        score += pwm[base][i]
    return score

def scan_sequence_for_regions_and_create_dataframe(full_sequence, window_size_35=6, window_size_10=6, gap_range=(14, 20), weight_35=1.5):
    """Scan the full DNA sequence for -35 and -10 regions, and create a DataFrame for both forward and reverse strands."""
    # 调用get_pwms获取PWM矩阵
    pwm_35, pwm_10 = get_pwms()

    # 正向链扫描
    best_result_forward = find_best_promoter_region(full_sequence, pwm_35, pwm_10, window_size_35, window_size_10, gap_range, strand='+', weight_35=weight_35)

    # 反向链扫描 (反向互补序列)
    reverse_complement_seq = str(Seq(full_sequence).reverse_complement())
    best_result_reverse = find_best_promoter_region(reverse_complement_seq, pwm_35, pwm_10, window_size_35, window_size_10, gap_range, strand='-', weight_35=weight_35)

    # 比较正向链和反向链的得分，选择最高的结果
    if best_result_forward["total score"] >= best_result_reverse["total score"]:
        best_result = best_result_forward
    else:
        best_result = best_result_reverse

    df = pd.DataFrame([best_result])
    return df

def find_best_promoter_region(sequence, pwm_35, pwm_10, window_size_35, window_size_10, gap_range, strand, weight_35):
    best_35 = None
    best_35_score = float('-inf')
    best_35_start = None
    best_result = None

    for i in range(len(sequence) - window_size_35 + 1):
        window_sequence_35 = sequence[i:i + window_size_35]
        score_35 = calculate_pwm_score(window_sequence_35, pwm_35)

        valid_10_found = False
        for gap in range(gap_range[0], gap_range[1] + 1):
            start_10 = i + window_size_35 + gap
            if start_10 + window_size_10 <= len(sequence):
                window_sequence_10 = sequence[start_10:start_10 + window_size_10]
                score_10 = calculate_pwm_score(window_sequence_10, pwm_10)
                
                # 使用权重调整-35区域的得分
                total_score = round_score(score_35 * weight_35 + score_10)

                if total_score > best_35_score:
                    best_35_score = total_score
                    best_result = {
                        "-35 sequence": window_sequence_35,
                        "-35 score": round_score(score_35),
                        "start position -35": i,
                        "end position -35": i + window_size_35,
                        "-10 sequence": window_sequence_10,
                        "-10 score": round_score(score_10),
                        "start position -10": start_10,
                        "end position -10": start_10 + window_size_10,
                        "total score": total_score,
                        "promoter strand": strand  # 添加链方向
                    }

    return best_result if best_result else {
        "-35 sequence": None,
        "-35 score": None,
        "start position -35": None,
        "end position -35": None,
        "-10 sequence": None,
        "-10 score": None,
        "start position -10": None,
        "end position -10": None,
        "total score": float('-inf'),
        "promoter strand": strand
    }
####重新定义find_matching_regions函数，返回链的数据####
def find_overlapping_region_relative(predicted_start, predicted_end, real_start, real_end):
    # 计算重叠区域的起始和结束位置（相对于基因组的绝对位置）
    overlap_start_absolute = max(predicted_start, real_start)
    overlap_end_absolute = min(predicted_end, real_end)

    # 计算重叠区域的相对位置
    # 对于预测序列来说，重叠区域的起始位置相对于预测序列的开始
    overlap_start_relative_to_predicted = overlap_start_absolute - predicted_start
    # 重叠区域的结束位置相对于预测序列的开始
    overlap_end_relative_to_predicted = overlap_end_absolute - predicted_start

    return overlap_start_relative_to_predicted, overlap_end_relative_to_predicted

def find_matching_regions_with_relative_info(predicted_seq, predicted_start, predicted_end, real_seq, real_start, real_end, strand):
    # 计算相对于预测序列的重叠区域的相对起始和结束位置
    overlap_start_relative, overlap_end_relative = find_overlapping_region_relative(predicted_start, predicted_end, real_start, real_end)

    if strand == '-':
        real_seq = complement(real_seq)  # 确保complement函数是正确的

    # 提取重叠区域的序列，使用相对位置信息
    overlapping_seq = predicted_seq[overlap_start_relative:overlap_end_relative + 1]  # +1因为Python切片不包括结束索引

    # 返回重叠区域的相对位置信息和序列
    return [(overlap_start_relative, overlap_end_relative, strand)]

def find_matching_regions(predicted_seq, real_seq, strand):#算法有问题，以修改为find_matching_regions_with_relative_info函数
    if strand == '-':
        real_seq = complement(real_seq)
    return [(i, i + len(real_seq), strand) for i in range(len(predicted_seq) - len(real_seq) + 1) if predicted_seq[i:i + len(real_seq)] == real_seq]

def complement(seq):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join([complement[base] for base in seq][::-1])#不拿到反向互补序列，而是直接反向

######################已经定义的函数######################
def calculate_pwm_minus(meme_results):
    complementary_sequences = [str(Seq(seq).reverse_complement()) for seq in meme_results['sequence']]
    m = motifs.create(complementary_sequences)
    pwm_minus_matrix = m.counts.normalize(pseudocounts={'A':0.6, 'C':0.6, 'G':0.6, 'T':0.6})
    pwm_minus_df = pd.DataFrame(pwm_minus_matrix)  # Convert to DataFrame
    return pwm_minus_df


def modify_tfbs_and_calculate_score(sequences, tfbs_positions, pwm_plus, pwm_minus, promoter_regions):
    """
    修改序列中的TFBS并计算转录因子结合的分数。
    """
    modified_sequences = []
    binding_scores = []
    
    for seq, (tfbs_start, tfbs_end, strand), promoter in zip(sequences, tfbs_positions, promoter_regions):
        start_35, end_35, start_10, end_10 = promoter
        pwm = pwm_minus if strand == '-' else pwm_plus
        lowest_scoring_bases = pwm.idxmin()
        modified_seq = list(seq)
        score = 0
        n = 0  # TFBS区域内的相对位置
        
        for i in range(tfbs_start, tfbs_end):
            original_base = seq[i]
            new_base = lowest_scoring_bases[n] if not (start_35 <= i < end_35 or start_10 <= i < end_10) else original_base
            modified_seq[i] = new_base

            # 计算改变后的碱基或原碱基的PWM分数
            score += pwm.loc[new_base, n]

            n += 1  # 更新TFBS区域内的相对位置
            if n >= len(pwm.columns):
                n = 0  # 重置迭代变量

        modified_sequences.append(''.join(modified_seq))
        binding_scores.append(score)

    return modified_sequences, binding_scores, promoter_regions



def generate_marked_sequence(original_seq, modified_seq):
    """
    生成标记修改的序列字符串。
    original_seq: 原始序列字符串。
    modified_seq: 修改后的序列字符串。
    """
    top_layer = original_seq
    middle_layer = []
    bottom_layer = modified_seq

    for orig_base, mod_base in zip(top_layer, bottom_layer):
        if orig_base != mod_base:
            middle_layer.append('|')  # 标记修改
        else:
            middle_layer.append(' ')  # 未修改

    # 将三层合并为一个字符串
    marked_seq = '\n'.join([''.join(top_layer), ''.join(middle_layer), ''.join(bottom_layer)])
    return marked_seq
#def apply_modifications_and_save_as_fasta(modifications, genbank_path, fasta_output_path):
    # 从GenBank文件读取基因组序列
   # record = SeqIO.read(genbank_path, "genbank")
    #genome_seq = str(record.seq)

    # 应用所有的修改
    #for (start, end), modified_seq in modifications.items():
     #   genome_seq = genome_seq[:start] + modified_seq + genome_seq[end:]

    # 创建新的SeqRecord
   # modified_record = SeqRecord(Seq(genome_seq), id=record.id, description="Modified Genome")

    # 保存修改后的基因组序列到FASTA文件
    #SeqIO.write(modified_record, fasta_output_path, "fasta")
def apply_modifications_and_save_as_fasta(modifications, genbank_path, fasta_output_path):
    """应用修改并保存为 FASTA 文件"""
    genome_record = SeqIO.read(genbank_path, "genbank")
    genome_seq = str(genome_record.seq)

    for (start, end), modified_seq in modifications.items():
        start_adj = start  # 如果需要，进行调整
        expected_length = end - start_adj  # 计算期望长度
        
        # 检查修改前后的长度是否匹配
        if len(modified_seq) != expected_length:
            raise ValueError(f"Modification length mismatch: {start}-{end} expected length {expected_length}, got {len(modified_seq)}")
        
        # 应用修改
        genome_seq = genome_seq[:start_adj] + modified_seq + genome_seq[end:]

    with open(fasta_output_path, "w") as output_file:
        output_file.write(f">{genome_record.id}\n")
        output_file.write("\n".join([genome_seq[i:i+60] for i in range(0, len(genome_seq), 60)]))

    print(f"Modified genome saved to {fasta_output_path}")
def get_sequence_id_from_genbank(genbank_path):
    record = SeqIO.read(genbank_path, "genbank")
    return record.id
def analyze_tfbs_modification(meme_results, df_promoters, pwm_df, genbank_path, output_path, window_size_35=6, window_size_10=6, gap_range=(14, 20), weight_35=1.5, extension=10):
    pwm_minus_df = calculate_pwm_minus(meme_results)
    sequence_id = get_sequence_id_from_genbank(genbank_path)
    
    # 读取基因组序列
    genome_record = SeqIO.read(genbank_path, "genbank")
    genome_seq = str(genome_record.seq)
    
    used_predicted_indices = set()
    used_real_indices = set()
    promoter_positions = []
    modifications = {}
    change_score_results = []  # 收集所有要打印的信息
    promoter_info_list = []  # 用于收集promoter和重叠的TFBS信息

    promoter_count = 0  # 用于跟踪启动子的个数

    # 更新df_promoters以包括延长后的区域
    df_promoters = df_promoters.copy()  # 创建副本以避免修改原始数据
    for idx, row in df_promoters.iterrows():
        extended_start, extended_end, extended_seq = extend_promoter_region(
            row['start'], row['end'], extension, genome_seq
        )
        df_promoters.at[idx, 'start'] = extended_start
        df_promoters.at[idx, 'end'] = extended_end
        df_promoters.at[idx, 'sequence'] = extended_seq

    for real_index, real in meme_results.iterrows():
        if real_index in used_real_indices:
            continue
            
        for predicted_index, predicted in df_promoters.iterrows():
            if predicted_index in used_predicted_indices:
                continue
            
            if (predicted['start'] <= real['end']) and (predicted['end'] >= real['start']):
                predicted_seq = predicted['sequence']
                real_seq = real['sequence']
                strand = real['strand']

                matching_regions = find_matching_regions_with_relative_info(predicted_seq, predicted['start'], predicted['end'], real_seq, real['start'], real['end'], strand)
                
                if matching_regions:
                    used_predicted_indices.add(predicted_index)
                    used_real_indices.add(real_index)

                    promoter_start = predicted['start']
                    promoter_end = predicted['end']
                    promoter_positions.append((promoter_start, promoter_end))

                    tfbs_positions = [(region[0], region[1], strand) for region in matching_regions]

                    # 使用整合的函数扫描并找到最佳的-35和-10区域
                    top_result = scan_sequence_for_regions_and_create_dataframe(predicted_seq, window_size_35, window_size_10, gap_range, weight_35)

                    if not top_result.empty:
                        top_result = top_result.iloc[0]
                        best_promoter_region = (top_result['start position -35'], top_result['start position -35'] + window_size_35,
                                                top_result['start position -10'], top_result['start position -10'] + window_size_10)
                        promoter_strand = top_result['promoter strand']  # 获取promoter strand信息
                        seq_35 = top_result['-35 sequence']  # 获取-35序列
                        seq_10 = top_result['-10 sequence']  # 获取-10序列
                    
                    modified_seqs, scores, promoter_regions = modify_tfbs_and_calculate_score(
                        [predicted_seq],
                        tfbs_positions,
                        pwm_df.T, pwm_minus_df.T,
                        [best_promoter_region]
                    )
                    key = (promoter_start, promoter_end)
                    modifications[key] = modified_seqs[0]
                    
                    marked_seq = generate_marked_sequence(predicted_seq, modified_seqs[0])                    
                    
                    promoter_count += 1  # 对每个找到的启动子递增计数器
                    promoter_header = f"Promoter {promoter_count}: Start = {promoter_start}, End = {promoter_end}"
                    change_score_results.append(promoter_header)
                    change_score_results.append(f"Modified Sequence with Scores:\n{marked_seq}")
                    change_score_results.append(f"Binding Scores: {scores[0]}")

                    # 查找与promoter重叠的meme_results中的TFBS
                    overlapping_tfbs = meme_results[(meme_results['start'] <= promoter_end) & (meme_results['end'] >= promoter_start)]

                    # 将promoter的-10和-35信息以及重叠的TFBS信息收集到列表中
                    for _, tfbs in overlapping_tfbs.iterrows():
                        # 计算 TFBS 相对于 Promoter 的重叠部分位置
                        relative_overlap_start = max(0, tfbs['start'] - promoter_start)
                        relative_overlap_end = min(promoter_end - promoter_start, tfbs['end'] - promoter_start)
                        overlap_seq = predicted_seq[relative_overlap_start:relative_overlap_end + 1]

                        promoter_info_list.append({
                            'Promoter': f'Promoter {promoter_count}',
                            '-35 Start': best_promoter_region[0],
                            '-35 End': best_promoter_region[1],
                            '-10 Start': best_promoter_region[2],
                            '-10 End': best_promoter_region[3],
                            '-35 Sequence': seq_35,
                            '-10 Sequence': seq_10,
                            'TFBS Sequence': tfbs['sequence'],
                            'TFBS Start': tfbs['start'],
                            'TFBS End': tfbs['end'],
                            'Relative Overlap Start': relative_overlap_start,
                            'Relative Overlap End': relative_overlap_end,
                            'Relative Overlap Sequence': overlap_seq,
                            'Promoter Strand': promoter_strand,  # 包含promoter strand信息
                            'Score': tfbs['score'],
                            'Sequence': predicted_seq
                        })

    # 将promoter信息转换为DataFrame
    promoter_info_df = pd.DataFrame(promoter_info_list)

    # 保存promoter信息为CSV文件
    promoter_info_csv_path = os.path.join(output_path, f"{sequence_id}_promoter_info.csv")
    promoter_info_df.to_csv(promoter_info_csv_path, index=False)
    print(f"Promoter information saved to {promoter_info_csv_path}")

    fasta_output_path = os.path.join(output_path, f"{sequence_id}_modified_genome.fasta")
    apply_modifications_and_save_as_fasta(modifications, genbank_path, fasta_output_path)

    change_score_result_path = os.path.join(output_path, f"{sequence_id}_change_score_result.txt")
    with open(change_score_result_path, "w") as file:
        file.write("\n".join(change_score_results))
    print(f"Change score results saved to {change_score_result_path}")

    return modifications, promoter_info_df
