import subprocess
from Bio import motifs
from Bio.Seq import Seq
import pandas as pd
import os
import matplotlib.pyplot as plt
import logomaker
import matplotlib.lines as mlines
import matplotlib.patches as patches
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

<<<<<<< HEAD
def scan_sequence_for_regions_and_create_dataframe(full_sequence, window_size_35=6, window_size_10=6,
                                                   gap_range=(14, 20), weight_35=1.5, strand='+'):
    """扫描全序列的 -35 和 -10 区域，创建包含正向和反向链结果的 DataFrame。"""
    # 获取 PWM 矩阵
    pwm_35, pwm_10 = get_pwms()
    seq_length = len(full_sequence)

    # 根据给定的 strand，决定初始扫描顺序
    if strand == '+':
        # 首先扫描正向链（原始序列）
        best_result_first = find_best_promoter_region(
            full_sequence, pwm_35, pwm_10, window_size_35, window_size_10, gap_range,
            strand='+', weight_35=weight_35, is_reverse=False
        )
        # 然后扫描反向链（反向互补序列）
        reverse_complement_seq = str(Seq(full_sequence).reverse_complement())
        best_result_second = find_best_promoter_region(
            reverse_complement_seq, pwm_35, pwm_10, window_size_35, window_size_10, gap_range,
            strand='-', weight_35=weight_35, is_reverse=True, seq_length=seq_length
        )
    elif strand == '-':
        # 首先扫描反向链（反向互补序列）
        reverse_complement_seq = str(Seq(full_sequence).reverse_complement())
        best_result_first = find_best_promoter_region(
            reverse_complement_seq, pwm_35, pwm_10, window_size_35, window_size_10, gap_range,
            strand='-', weight_35=weight_35, is_reverse=True, seq_length=seq_length
        )
        # 然后扫描正向链（原始序列）
        best_result_second = find_best_promoter_region(
            full_sequence, pwm_35, pwm_10, window_size_35, window_size_10, gap_range,
            strand='+', weight_35=weight_35, is_reverse=False
        )
    else:
        # 如果 strand 无效，默认使用 '+'
        best_result_first = find_best_promoter_region(
            full_sequence, pwm_35, pwm_10, window_size_35, window_size_10, gap_range,
            strand='+', weight_35=weight_35, is_reverse=False
        )
        reverse_complement_seq = str(Seq(full_sequence).reverse_complement())
        best_result_second = find_best_promoter_region(
            reverse_complement_seq, pwm_35, pwm_10, window_size_35, window_size_10, gap_range,
            strand='-', weight_35=weight_35, is_reverse=True, seq_length=seq_length
        )

    # 比较得分，选择更高的结果
    if best_result_first["total score"] >= best_result_second["total score"]:
        best_result = best_result_first
    else:
        best_result = best_result_second
        # 更新 strand
        strand = best_result['promoter strand']

    df = pd.DataFrame([best_result])
    return df, strand  # 返回更新后的 strand

def find_best_promoter_region(sequence, pwm_35, pwm_10, window_size_35, window_size_10,
                              gap_range, strand, weight_35, is_reverse=False, seq_length=None):
    best_result = None
    best_total_score = float('-inf')
=======
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
>>>>>>> 03027939e00730d23a5515fb8f7fcf63b107702a

    for i in range(len(sequence) - window_size_35 + 1):
        window_sequence_35 = sequence[i:i + window_size_35]
        score_35 = calculate_pwm_score(window_sequence_35, pwm_35)

<<<<<<< HEAD
=======
        valid_10_found = False
>>>>>>> 03027939e00730d23a5515fb8f7fcf63b107702a
        for gap in range(gap_range[0], gap_range[1] + 1):
            start_10 = i + window_size_35 + gap
            if start_10 + window_size_10 <= len(sequence):
                window_sequence_10 = sequence[start_10:start_10 + window_size_10]
                score_10 = calculate_pwm_score(window_sequence_10, pwm_10)
<<<<<<< HEAD

                # 使用权重调整 -35 区域的得分
                total_score = score_35 * weight_35 + score_10

                if total_score > best_total_score:
                    best_total_score = total_score
                    start_position_35 = i
                    end_position_35 = i + window_size_35
                    start_position_10 = start_10
                    end_position_10 = start_10 + window_size_10

                    # 如果是在反向互补序列中，需要将位置映射回原始序列
                    if is_reverse and seq_length is not None:
                        start_position_35_mapped = seq_length - end_position_35
                        end_position_35_mapped = seq_length - start_position_35
                        start_position_10_mapped = seq_length - end_position_10
                        end_position_10_mapped = seq_length - start_position_10
                    else:
                        start_position_35_mapped = start_position_35
                        end_position_35_mapped = end_position_35
                        start_position_10_mapped = start_position_10
                        end_position_10_mapped = end_position_10

                    best_result = {
                        "-35 sequence": window_sequence_35,
                        "-35 score": round_score(score_35),
                        "start position -35": start_position_35_mapped,
                        "end position -35": end_position_35_mapped,
                        "-10 sequence": window_sequence_10,
                        "-10 score": round_score(score_10),
                        "start position -10": start_position_10_mapped,
                        "end position -10": end_position_10_mapped,
                        "total score": round_score(total_score),
                        "promoter strand": strand  # 添加链方向
                    }

    # 如果没有找到任何结果，返回空结果
    if best_result is None:
        best_result = {
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

    return best_result
=======
                
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
>>>>>>> 03027939e00730d23a5515fb8f7fcf63b107702a
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
def plot_sequence_with_annotations(predicted_seq, modified_seq, best_promoter_region, tfbs_regions,
                                   predicted_index, real_index, output_dir):
    """
    绘制序列，突出显示 -35、-10、TFBS 区域。

    - 对 -35 和 -10 区域，只标记方框，不添加连接线。
    - 对 TFBS 区域，标记方框，并添加连接线。
    - 在下方的序列中，仅展示 TFBS 区域的碱基变化，其他位置不显示。

    参数：
    - predicted_seq: 原始序列
    - modified_seq: 改造后的序列
    - best_promoter_region: (-35_start, -35_end, -10_start, -10_end)
    - tfbs_regions: [(tfbs_start, tfbs_end, strand), ...]
    - predicted_index: 预测序列的索引
    - real_index: 实际序列的索引
    - output_dir: 输出目录
    """
    seq_length = len(predicted_seq)
    
    # 调整碱基之间的间距
    x_scale = 0.5  # 可以调整此值以增加或减少间距
    x_positions = [i * x_scale for i in range(seq_length)]
    fig_width = max(10, x_positions[-1] + 2)
    if fig_width > 20:
        fig_width = 20  # 限制图像的最大宽度为20英寸

    fig_height = 4  # 调整图像高度

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))  # 设置图像尺寸

    # 绘制上方的原始序列
    for i, nucleotide in enumerate(predicted_seq):
        ax.text(x_positions[i], 1.5, nucleotide, ha='center', va='center', fontsize=10, color='black')

    # 在下方的序列中，仅展示 TFBS 区域的碱基变化，其他位置不显示
    for i in range(seq_length):
        in_tfbs = False
        for tfbs_start, tfbs_end, strand in tfbs_regions:
            if tfbs_start <= i <= tfbs_end:
                in_tfbs = True
                break
        if in_tfbs:
            nucleotide = modified_seq[i]
            ax.text(x_positions[i], 0.5, nucleotide, ha='center', va='center', fontsize=10, color='black')
        else:
            # 不显示非 TFBS 区域的碱基
            pass

    # 添加连接线，仅在 TFBS 区域
    line_positions = set()
    # 收集 TFBS 区域的位置信息
    for tfbs_start, tfbs_end, strand in tfbs_regions:
        for pos in range(tfbs_start, tfbs_end + 1):
            if 0 <= pos < seq_length:
                line_positions.add(pos)

    # 添加连接线
    for i in line_positions:
        if predicted_seq[i] != modified_seq[i]:
            line_color = 'red'  # 修改的位置用红色线连接
            line_width = 1.5
        else:
            line_color = 'gray'  # 未修改的位置用灰色线连接
            line_width = 0.5
        ax.add_line(mlines.Line2D([x_positions[i], x_positions[i]], [1.4, 0.6], color=line_color, linewidth=line_width))

    # 添加颜色方框标记区域
    # -35 区域，蓝色方框
    start_35, end_35 = best_promoter_region[0], best_promoter_region[1]
    rect_35_top = patches.Rectangle((x_positions[start_35] - x_scale/2, 1.2),
                                    x_scale * (end_35 - start_35), 0.6,
                                    linewidth=1, edgecolor='blue', facecolor='none')
    ax.add_patch(rect_35_top)

    # -10 区域，橙色方框
    start_10, end_10 = best_promoter_region[2], best_promoter_region[3]
    rect_10_top = patches.Rectangle((x_positions[start_10] - x_scale/2, 1.2),
                                    x_scale * (end_10 - start_10), 0.6,
                                    linewidth=1, edgecolor='orange', facecolor='none')
    ax.add_patch(rect_10_top)

    # TFBS 区域，红色方框，上下方序列都标记
    for tfbs_start, tfbs_end, strand in tfbs_regions:
        # 上方序列的方框
        rect_tfbs_top = patches.Rectangle((x_positions[tfbs_start] - x_scale/2, 1.2),
                                          x_scale * (tfbs_end - tfbs_start + 1), 0.6,
                                          linewidth=1, edgecolor='red', facecolor='none')
        ax.add_patch(rect_tfbs_top)
        # 下方序列的方框
        rect_tfbs_bottom = patches.Rectangle((x_positions[tfbs_start] - x_scale/2, 0.2),
                                             x_scale * (tfbs_end - tfbs_start + 1), 0.6,
                                             linewidth=1, edgecolor='red', facecolor='none')
        ax.add_patch(rect_tfbs_bottom)

    # 添加 Promoter 名称和索引信息，居中显示
    ax.text(x_positions[seq_length // 2], 2.0,
            f'Predicted Index: {predicted_index}, Real Index: {real_index}',
            ha='center', va='center', fontsize=12, color='black')

    # 调整图例
    legend_elements = [
        patches.Patch(edgecolor='blue', facecolor='none', label='-35 Region'),
        patches.Patch(edgecolor='orange', facecolor='none', label='-10 Region'),
        patches.Patch(edgecolor='red', facecolor='none', label='TFBS Region'),
        mlines.Line2D([], [], color='red', linewidth=1.5, label='Modified Base'),
        mlines.Line2D([], [], color='gray', linewidth=0.5, label='Unchanged Base')
    ]
    ax.legend(handles=legend_elements, loc='upper center', bbox_to_anchor=(0.5, -0.1),
              ncol=3, fontsize=10, frameon=False)

    # 调整图形设置
    ax.set_xlim(-1, x_positions[-1] + 1)
    ax.set_ylim(0, 2.5)
    ax.axis('off')

    # 调整图形的边距，去除多余空白
    plt.subplots_adjust(top=0.85, bottom=0.15)

    # 保存并关闭图像，设置 dpi=300
    image_filename = f"sequence_match_{predicted_index}_{real_index}.png"
    save_path = os.path.join(output_dir, image_filename)
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"图像已保存到 {save_path}")
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
<<<<<<< HEAD

def extend_promoter_region(promoter_start, promoter_end, extension, genome_seq):
    """延长 Promoter 区域并提取延长后的序列"""
    extended_start = max(0, promoter_start - extension)
    extended_end = min(promoter_end + extension, len(genome_seq))
    extended_seq = genome_seq[extended_start:extended_end]
    return extended_start, extended_end, extended_seq
def plot_sequences_from_dataframe(df, modifications_dict, output_dir):
    """
    从 DataFrame 中读取数据，并绘制序列图像。

    参数：
    - df: 包含 promoter 信息的 DataFrame
    - modifications_dict: 一个字典，键为 (promoter_start, promoter_end)，值为改造后的序列
    - output_dir: 输出图像的目录
    """
    # 创建输出目录（如果不存在）
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # 按照 'Promoter' 列进行分组
    grouped = df.groupby('Promoter')
    for promoter_name, group in grouped:
        # 获取序列和位置信息
        sequence = group['Sequence'].iloc[0]
        seq_length = len(sequence)
        seq_indices = range(seq_length)
        
        # 获取 -35 和 -10 区域的位置
        start_35 = int(group['-35 Start'].iloc[0])
        end_35 = int(group['-35 End'].iloc[0])
        start_10 = int(group['-10 Start'].iloc[0])
        end_10 = int(group['-10 End'].iloc[0])
        
        # 获取 TFBS 区域的位置
        tfbs_starts = group['Relative Overlap Start'].astype(int).tolist()
        tfbs_ends = group['Relative Overlap End'].astype(int).tolist()
        
        # 获取改造后的序列
        promoter_start = int(group['TFBS Start'].min())
        promoter_end = int(group['TFBS End'].max())
        key = (promoter_start, promoter_end)
        modified_sequence = modifications_dict.get(key, sequence)
        
        # 创建颜色列表
        nucleotide_colors = ['black'] * seq_length
        
        # 标记 -35 区域为蓝色
        for i in range(start_35, end_35):
            if 0 <= i < seq_length:
                nucleotide_colors[i] = 'blue'
        
        # 标记 -10 区域为黄色
        for i in range(start_10, end_10):
            if 0 <= i < seq_length:
                nucleotide_colors[i] = 'yellow'
        
        # 标记 TFBS 区域为红色
        for start, end in zip(tfbs_starts, tfbs_ends):
            for i in range(start, end + 1):
                if 0 <= i < seq_length:
                    nucleotide_colors[i] = 'red'
        
        # 标记改造的碱基为绿色
        for i, (orig_base, mod_base) in enumerate(zip(sequence, modified_sequence)):
            if orig_base != mod_base:
                nucleotide_colors[i] = 'green'
        
        # 绘制序列
        fig, ax = plt.subplots(figsize=(10, 3))
        for i, (nucleotide, color) in enumerate(zip(sequence, nucleotide_colors)):
            ax.text(i, 1, nucleotide, ha='center', va='center', fontsize=10, color=color)
        
        # 绘制改造后的序列
        for i, nucleotide in enumerate(modified_sequence):
            ax.text(i, 0, nucleotide, ha='center', va='center', fontsize=10, color=nucleotide_colors[i])
        
        # 添加 Promoter 名称
        ax.text(0, 1.5, promoter_name, ha='left', va='center', fontsize=12, color='black')
        
        # 设置轴
        ax.set_xlim(-1, seq_length + 1)
        ax.set_ylim(-1, 2)
        ax.axis('off')
        plt.tight_layout()
        
        # 保存图像
        image_filename = f"{promoter_name.replace(' ', '_')}.png"
        save_path = os.path.join(output_dir, image_filename)
        plt.savefig(save_path)
        plt.close()
        print(f"{promoter_name} 的图像已保存到 {save_path}")


def analyze_tfbs_modification(meme_results, df_promoters, pwm_df, genbank_path, output_path, window_size_35=6, window_size_10=6, gap_range=(14, 20), weight_35=1.5, extension=10):
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    
    pwm_minus_df = calculate_pwm_minus(meme_results)
    sequence_id = get_sequence_id_from_genbank(genbank_path)
    
    # 获取 pwm_35 和 pwm_10
    pwm_35, pwm_10 = get_pwms()
    
=======
def analyze_tfbs_modification(meme_results, df_promoters, pwm_df, genbank_path, output_path, window_size_35=6, window_size_10=6, gap_range=(14, 20), weight_35=1.5, extension=10):
    pwm_minus_df = calculate_pwm_minus(meme_results)
    sequence_id = get_sequence_id_from_genbank(genbank_path)
    
>>>>>>> 03027939e00730d23a5515fb8f7fcf63b107702a
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
<<<<<<< HEAD
                strand = real['strand']  # 获取初始的 strand
                print('predicted_seq, strand:', predicted_seq, strand)
                matching_regions = find_matching_regions_with_relative_info(
                    predicted_seq, predicted['start'], predicted['end'],
                    real_seq, real['start'], real['end'], strand
                )
=======
                strand = real['strand']

                matching_regions = find_matching_regions_with_relative_info(predicted_seq, predicted['start'], predicted['end'], real_seq, real['start'], real['end'], strand)
>>>>>>> 03027939e00730d23a5515fb8f7fcf63b107702a
                
                if matching_regions:
                    used_predicted_indices.add(predicted_index)
                    used_real_indices.add(real_index)
<<<<<<< HEAD
    
                    promoter_start = predicted['start']
                    promoter_end = predicted['end']
                    promoter_positions.append((promoter_start, promoter_end))
    
                    tfbs_positions = [(region[0], region[1], strand) for region in matching_regions]
    
                    # 使用整合的函数扫描并找到最佳的 -35 和 -10 区域
                    top_result, updated_strand = scan_sequence_for_regions_and_create_dataframe(
                        predicted_seq, window_size_35, window_size_10, gap_range, weight_35, strand=strand
                    )
                    # 更新 strand
                    strand = updated_strand

                    if not top_result.empty:
                        top_result = top_result.iloc[0]
                        best_promoter_region = (int(top_result['start position -35']), int(top_result['end position -35']),
                                                int(top_result['start position -10']), int(top_result['end position -10']))
=======

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
>>>>>>> 03027939e00730d23a5515fb8f7fcf63b107702a
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
                    
<<<<<<< HEAD
                    # 生成标记的序列，显示修改位置
=======
>>>>>>> 03027939e00730d23a5515fb8f7fcf63b107702a
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

<<<<<<< HEAD
                    # **在此处调用绘图函数，传入 marked_seq**
                    plot_sequence_with_annotations(
                        predicted_seq=predicted_seq,
                        modified_seq=modified_seqs[0],
                        best_promoter_region=best_promoter_region,
                        tfbs_regions=tfbs_positions,
                        predicted_index=predicted_index,
                        real_index=real_index,
                     #   pwm_35=pwm_35,
                     #   pwm_10=pwm_10,
                        output_dir=output_path,
                       # marked_seq=marked_seq  # 传入标记的序列
                    )

=======
>>>>>>> 03027939e00730d23a5515fb8f7fcf63b107702a
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
<<<<<<< HEAD


=======
>>>>>>> 03027939e00730d23a5515fb8f7fcf63b107702a
