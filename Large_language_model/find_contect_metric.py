import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO


##  按照基因上的顺序寻找合适的蛋白
# 定义函数：从GenBank文件中提取目标蛋白质的核苷酸序列、氨基酸序列和所有基因间序列
# def extract_protein_and_intergenic_sequences(gbk_file, temp_dir, protein_names, pfam_hmm, exclude_keywords=None):
#     """
#     从GenBank文件中提取第一个符合指定条件的蛋白质的核苷酸序列（对于负链已反向互补）和氨基酸序列，
#     以及所有基因间序列，并保存到指定的文件中。

#     参数：
#     - gbk_file (str): 输入的GenBank文件路径。
#     - temp_dir (str): 用于存放临时文件的目录。
#     - protein_names (list): 需要提取的蛋白质名称列表。
#     - pfam_hmm (str): Pfam HMM数据库的路径。
#     - exclude_keywords (list, optional): 排除的关键词列表。

#     返回：
#     - nucleotide_seq (str): 目标蛋白质的核苷酸序列（对于负链，已反向互补）。
#     - protein_seq (str): 目标蛋白质的氨基酸序列。
#     - intergenic_sequences (list of str): 所有基因间的核苷酸序列列表。
#     - domain_nt_positions (tuple): HTH结构域在核苷酸序列中的起始和结束位置（已扩展20 bp）。
#     - strand (int): 目标基因的链方向（+1 或 -1）。
#     """
#     from Bio.Seq import Seq
#     # 设置默认的排除关键词
#     if exclude_keywords is None:
#         exclude_keywords = ['Cro', 'anti', 'cro', 'CRO']

#     # 检查输入的GenBank文件是否存在
#     if not os.path.exists(gbk_file):
#         print(f"文件 {gbk_file} 不存在。")
#         return None, None, None, None, None

#     # 解析GenBank文件
#     records = list(SeqIO.parse(gbk_file, "genbank"))
#     if not records:
#         print(f"无法解析GenBank文件 {gbk_file}。")
#         return None, None, None, None, None

#     record = records[0]  # 假设只有一个记录

#     # 收集所有CDS特征
#     cds_features = [feature for feature in record.features if feature.type == "CDS" and "product" in feature.qualifiers]
#     if not cds_features:
#         print("GenBank文件中没有CDS特征。")
#         return None, None, None, None, None

#     # 对CDS特征按起始位置排序
#     cds_features_sorted = sorted(cds_features, key=lambda f: f.location.start)
#     # 寻找第一个符合条件的蛋白质
#     matching_cds = None
#     for feature in cds_features_sorted:
#         product = feature.qualifiers["product"][0]
#         print('product',product)
#         # 排除包含不需要的关键词的蛋白质（不区分大小写）
#         if any(keyword.lower() in product.lower() for keyword in exclude_keywords):
#             continue

# 检查蛋白质名称是否在指定列表中（不区分大小写）
#         if any(protein_name.lower() in product.lower() for protein_name in protein_names):
#             matching_cds = feature
#             break  # 找到第一个匹配的CDS，退出循环

#     if not matching_cds:
#         print("未找到符合条件的蛋白质。")
#         return None, None, None, None, None

#     # 获取基因位置信息
#     location = matching_cds.location
#     gene_start = int(location.start)  # Biopython使用0-based索引
#     gene_end = int(location.end)
#     strand = location.strand  # +1 或 -1

#     # 提取匹配CDS的核苷酸序列
#     try:
#         # 提取核苷酸序列（保持原始方向）
#         nucleotide_seq = matching_cds.extract(record.seq)

#         # 获取翻译后的氨基酸序列
#         protein_seq = matching_cds.qualifiers.get("translation", [""])[0]
#         if not protein_seq:
#             # 根据链方向翻译蛋白质序列
#             if strand == 1:
#                 protein_seq = nucleotide_seq.translate(to_stop=True)
#             elif strand == -1:
#                 protein_seq = nucleotide_seq.reverse_complement().translate(to_stop=True)
#             else:
#                 print("未知的链方向。")
#                 return None, None, None, None, None
#     except Exception as e:
#         print(f"提取序列时出错: {e}")
#         return None, None, None, None, None

#     if not nucleotide_seq or not protein_seq:
#         print("未能提取到蛋白质的核苷酸序列或氨基酸序列。")
#         return None, None, None, None, None

#     # 准备FASTA条目
#     strand_symbol = '+' if strand == 1 else '-'

#     # 对负链的序列进行反向互补操作，然后保存到FASTA文件
#     fasta_sequence_protein_nucleotide = str(nucleotide_seq)
#     if strand == -1:
#         fasta_sequence_protein_nucleotide = str(Seq(fasta_sequence_protein_nucleotide).reverse_complement())

#     fasta_header_protein_nucleotide = f">{record.id}_{matching_cds.qualifiers['product'][0].replace(' ', '_')}_nucleotide|Location:{gene_start+1}-{gene_end}({strand_symbol})"

#     # 保存蛋白质的核苷酸序列到FASTA文件
#     output_fasta = os.path.join(temp_dir, "extracted_sequences.fasta")
#     fasta_entries = [f"{fasta_header_protein_nucleotide}\n{fasta_sequence_protein_nucleotide}\n"]

#     # 收集基因间序列（保持原始方向，不进行反向互补）
#     intergenic_sequences = []
#     intergenic_counter = 1
#     for i in range(len(cds_features_sorted) - 1):
#         current_cds = cds_features_sorted[i]
#         next_cds = cds_features_sorted[i + 1]

#         # 定义基因间区域
#         intergenic_start = current_cds.location.end
#         intergenic_end = next_cds.location.start

#         if intergenic_end > intergenic_start:
#             intergenic_seq = record.seq[intergenic_start:intergenic_end]
#             intergenic_header = f">{record.id}_Intergenic_{intergenic_counter}|Location:{intergenic_start+1}-{intergenic_end}(+)"
#             intergenic_sequence = str(intergenic_seq)
#             fasta_entries.append(f"{intergenic_header}\n{intergenic_sequence}\n")
#             intergenic_sequences.append(intergenic_sequence)
#             intergenic_counter += 1

#     # 将所有序列保存到FASTA文件
#     with open(output_fasta, "w") as fasta_out:
#         fasta_out.writelines(fasta_entries)
#     print(f"序列已保存到 {output_fasta}")

#     # 运行hmmscan，寻找HTH结构域
#     domain_nt_positions = find_helix_turn_helix_domain(
#         nucleotide_seq=nucleotide_seq,
#         protein_seq=protein_seq,
#         pfam_hmm=pfam_hmm,
#         temp_dir=temp_dir,
#         strand=strand
#     )

#     if domain_nt_positions is None:
#         print("未能找到HTH结构域，程序终止。")
#         return None, None, None, None, None

#     return str(nucleotide_seq), str(protein_seq), intergenic_sequences, domain_nt_positions, strand

# 按照列表的顺序寻找合适的蛋白
def extract_protein_and_intergenic_sequences(gbk_file, temp_dir, protein_names, pfam_hmm, exclude_keywords=None):
    """
    从GenBank文件中提取符合指定条件的蛋白质的核苷酸序列（对于负链已反向互补）和氨基酸序列，
    以及所有基因间序列，并保存到指定的文件中。

    参数：
    - gbk_file (str): 输入的GenBank文件路径。
    - temp_dir (str): 用于存放临时文件的目录。
    - protein_names (list): 需要提取的蛋白质名称列表。
    - pfam_hmm (str): Pfam HMM数据库的路径。
    - exclude_keywords (list, optional): 排除的关键词列表。

    返回：
    - nucleotide_seq (str): 目标蛋白质的核苷酸序列（对于负链，已反向互补）。
    - protein_seq (str): 目标蛋白质的氨基酸序列。
    - intergenic_sequences (list of str): 所有基因间的核苷酸序列列表。
    - domain_nt_positions (tuple): HTH结构域在核苷酸序列中的起始和结束位置（已扩展20 bp）。
    - strand (int): 目标基因的链方向（+1 或 -1）。
    """
    from Bio.Seq import Seq
    # 设置默认的排除关键词
    if exclude_keywords is None:
        exclude_keywords = ['Cro', 'anti', 'cro', 'CRO']

    # 检查输入的GenBank文件是否存在
    if not os.path.exists(gbk_file):
        print(f"文件 {gbk_file} 不存在。")
        return None, None, None, None, None

    # 解析GenBank文件
    records = list(SeqIO.parse(gbk_file, "genbank"))
    if not records:
        print(f"无法解析GenBank文件 {gbk_file}。")
        return None, None, None, None, None

    record = records[0]  # 假设只有一个记录

    # 收集所有CDS特征
    cds_features = [feature for feature in record.features if feature.type == "CDS" and "product" in feature.qualifiers]
    if not cds_features:
        print("GenBank文件中没有CDS特征。")
        return None, None, None, None, None

    # 对CDS特征按起始位置排序
    cds_features_sorted = sorted(cds_features, key=lambda f: f.location.start)
    # 寻找符合条件的蛋白质，按照 protein_names 的顺序查找
    matching_cds = None
    for protein_name in protein_names:
        for feature in cds_features_sorted:
            product = feature.qualifiers["product"][0]
            # 排除包含不需要的关键词的蛋白质（不区分大小写）
            if any(keyword.lower() in product.lower() for keyword in exclude_keywords):
                continue

            # 检查蛋白质名称是否与当前的 protein_name 匹配（不区分大小写）
            if protein_name.lower() in product.lower():
                matching_cds = feature
                break  # 找到匹配的CDS，退出内循环
        if matching_cds:
            break  # 找到匹配的CDS，退出外循环

    if not matching_cds:
        print("未找到符合条件的蛋白质。")
        return None, None, None, None, None

    # 获取基因位置信息
    location = matching_cds.location
    gene_start = int(location.start)
    gene_end = int(location.end)
    strand = location.strand

    # 提取匹配CDS的核苷酸序列
    try:
        nucleotide_seq = matching_cds.extract(record.seq)

        # 获取翻译后的氨基酸序列
        protein_seq = matching_cds.qualifiers.get("translation", [""])[0]
        if not protein_seq:
            # 根据链方向翻译蛋白质序列
            if strand == 1:
                protein_seq = nucleotide_seq.translate(to_stop=True)
            elif strand == -1:
                protein_seq = nucleotide_seq.reverse_complement().translate(to_stop=True)
            else:
                print("未知的链方向。")
                return None, None, None, None, None
    except Exception as e:
        print(f"提取序列时出错: {e}")
        return None, None, None, None, None

    if not nucleotide_seq or not protein_seq:
        print("未能提取到蛋白质的核苷酸序列或氨基酸序列。")
        return None, None, None, None, None

    # 准备FASTA条目
    strand_symbol = '+' if strand == 1 else '-'
    fasta_sequence_protein_nucleotide = str(nucleotide_seq)
    if strand == -1:
        fasta_sequence_protein_nucleotide = str(Seq(fasta_sequence_protein_nucleotide).reverse_complement())

    fasta_header_protein_nucleotide = f">{record.id}_{matching_cds.qualifiers['product'][0].replace(' ', '_')}_nucleotide|Location:{gene_start + 1}-{gene_end}({strand_symbol})"

    # 保存蛋白质的核苷酸序列到FASTA文件
    output_fasta = os.path.join(temp_dir, "extracted_sequences.fasta")
    fasta_entries = [f"{fasta_header_protein_nucleotide}\n{fasta_sequence_protein_nucleotide}\n"]

    # 收集基因间序列
    intergenic_sequences = []
    intergenic_counter = 1
    for i in range(len(cds_features_sorted) - 1):
        current_cds = cds_features_sorted[i]
        next_cds = cds_features_sorted[i + 1]
        intergenic_start = current_cds.location.end
        intergenic_end = next_cds.location.start

        if intergenic_end > intergenic_start:
            intergenic_seq = record.seq[intergenic_start:intergenic_end]
            intergenic_header = f">{record.id}_Intergenic_{intergenic_counter}|Location:{intergenic_start + 1}-{intergenic_end}(+)"
            intergenic_sequence = str(intergenic_seq)
            fasta_entries.append(f"{intergenic_header}\n{intergenic_sequence}\n")
            intergenic_sequences.append(intergenic_sequence)
            intergenic_counter += 1

    # 将所有序列保存到FASTA文件
    with open(output_fasta, "w") as fasta_out:
        fasta_out.writelines(fasta_entries)
    print(f"序列已保存到 {output_fasta}")

    # 运行hmmscan，寻找HTH结构域
    domain_nt_positions = find_helix_turn_helix_domain(
        nucleotide_seq=nucleotide_seq,
        protein_seq=protein_seq,
        pfam_hmm=pfam_hmm,
        temp_dir=temp_dir,
        strand=strand
    )

    if domain_nt_positions is None:
        print("未能找到HTH结构域，程序终止。")
        return None, None, None, None, None

    return str(nucleotide_seq), str(protein_seq), intergenic_sequences, domain_nt_positions, strand


# 定义函数：使用hmmscan寻找HTH结构域
def find_helix_turn_helix_domain(nucleotide_seq, protein_seq, pfam_hmm, temp_dir, strand):
    """
    使用hmmscan在蛋白质序列中寻找多个抑制蛋白相关的结构域，并根据基因方向计算其核苷酸位置。

    参数：
    - nucleotide_seq (str): 目标蛋白质的核苷酸序列（对于负链，已反向互补）。
    - protein_seq (str): 蛋白质的氨基酸序列。
    - pfam_hmm (str): Pfam HMM数据库的路径。
    - temp_dir (str): 用于存放临时文件的目录。
    - strand (int): 基因的链方向（+1 或 -1）。

    返回：
    - domain_nt_positions (tuple): HTH结构域在核苷酸序列中的起始和结束位置（已扩展20 bp）。
    """
    import subprocess
    from Bio.Seq import Seq

    # 定义目标结构域描述的集合
    target_domain_descriptions = [
        'Helix-turn-helix',
        'Bacterial regulatory proteins',
        'Repressor',
        'Transcriptional repressor',
        'HTH domain',
        'Regulatory protein',
        'DNA-binding domain',
        'Arc-like DNA binding domain'
    ]

    # 将蛋白质序列保存到临时FASTA文件
    protein_fasta = os.path.join(temp_dir, "protein_seq.fasta")
    with open(protein_fasta, "w") as f:
        f.write(f">protein_of_interest\n{protein_seq}\n")

    # 运行hmmscan
    hmmscan_output = os.path.join(temp_dir, "hmmscan_output.txt")
    command = [
        "hmmscan",
        "--domtblout", hmmscan_output,
        pfam_hmm,
        protein_fasta
    ]

    try:
        print("正在运行hmmscan...")
        subprocess.run(command, check=True)
        print("hmmscan运行完成。")
    except subprocess.CalledProcessError as e:
        print(f"运行hmmscan时出错: {e}")
        return None

    # 解析hmmscan输出，寻找目标结构域
    target_domains = []
    with open(hmmscan_output, "r") as infile:
        for line in infile:
            if line.startswith("#"):
                continue
            columns = line.strip().split()
            if len(columns) < 23:
                continue
            description = ' '.join(columns[22:])
            # 检查描述中是否包含任何目标结构域描述
            if any(target_desc.lower() in description.lower() for target_desc in target_domain_descriptions):
                domain_score = float(columns[13])
                ali_from = int(columns[17])
                ali_to = int(columns[18])
                target_domains.append({
                    'description': description,
                    'domain_score': domain_score,
                    'ali_from': ali_from,
                    'ali_to': ali_to
                })

    if not target_domains:
        print("在hmmscan结果中未找到目标结构域。")
        return None

    # 选择得分最高的结构域
    best_domain = max(target_domains, key=lambda x: x['domain_score'])
    print(f"找到的最佳结构域: {best_domain['description']}, 得分: {best_domain['domain_score']}")
    print(f"氨基酸位置: {best_domain['ali_from']}-{best_domain['ali_to']}")

    gene_length = len(nucleotide_seq)

    # 计算结构域在核苷酸序列中的位置
    if strand == 1:
        # 正链，从前向后
        nt_start = (best_domain['ali_from'] - 1) * 3
        nt_end = best_domain['ali_to'] * 3 - 1
    elif strand == -1:
        # 负链，从后向前，需要反向计算
        nt_end = gene_length - (best_domain['ali_from'] - 1) * 3 - 1
        nt_start = gene_length - best_domain['ali_to'] * 3
    else:
        print("未知的链方向。")
        return None

    # 扩展20 bp，并确保不越界
    domain_nt_start_in_seq = max(0, nt_start - 20)
    domain_nt_end_in_seq = min(gene_length - 1, nt_end + 20)

    print(f"结构域在核苷酸序列中的位置（扩展20 bp）: {domain_nt_start_in_seq}-{domain_nt_end_in_seq}")

    return (domain_nt_start_in_seq, domain_nt_end_in_seq)


# 定义输入和输出路径
gbk_file = "/home/hanzequan/test_set_megadna/Bxb1.gbk"
pfam_hmm = "/home/hanzequan/test_set_megadna/Pfam/Pfam-A.hmm"
temp_dir = "/home/hanzequan/test_set_megadna/temp_dir"
output_dir = "/home/hanzequan/test_set_megadna/batch_figure"
output_fasta = "/home/hanzequan/test_set_megadna/temp_dir/extracted_sequences.fasta"

# 创建临时目录
os.makedirs(temp_dir, exist_ok=True)
os.makedirs(output_dir, exist_ok=True)

# 定义感兴趣的蛋白质名称列表和排除关键词
protein_names = protein_names
exclude_keywords = ["Cro", "anti", "cro", "CRO"]
# 提取目标蛋白质的序列和HTH结构域的位置
nucleotide_seq, protein_seq, intergenic_sequences, domain_nt_positions, strand = extract_protein_and_intergenic_sequences(
    gbk_file=gbk_file,
    temp_dir=temp_dir,
    protein_names=protein_names,
    pfam_hmm=pfam_hmm,
    exclude_keywords=exclude_keywords
)
if not nucleotide_seq or not protein_seq:
    print("未能提取到序列，程序终止。")
    exit(1)

# 初始化 Alphabet 实例
alphabet = Alphabet()

# 读取 extracted_sequences.fasta 文件中的所有序列
sequences = []
sequence_sources = []  # 新增一个列表，记录序列的来源
for record in SeqIO.parse(output_fasta, 'fasta'):
    sequences.append(str(record.seq))
    if 'Intergenic' in record.id:
        sequence_sources.append('intergenic')
    else:
        sequence_sources.append('protein')

# 检查是否至少有一个序列
if not nucleotide_seq or not protein_seq:
    print("未能提取到序列，程序终止。")
    exit(1)

# 初始化 Alphabet 实例
alphabet = Alphabet()

# 读取 extracted_sequences.fasta 文件中的所有序列
sequences = []
sequence_sources = []  # 新增一个列表，记录序列的来源
for record in SeqIO.parse(output_fasta, 'fasta'):
    sequences.append(str(record.seq))
    if 'Intergenic' in record.id:
        sequence_sources.append('intergenic')
    else:
        sequence_sources.append('protein')

# 检查是否至少有一个序列
if not sequences:
    print(f"在 '{output_fasta}' 中未找到任何序列。")
    exit(1)
else:
    # 获取第一个序列及其长度（抑制蛋白序列）
    first_sequence = sequences[0]
    first_sequence_length = len(first_sequence)

    # 初始化已处理的序列索引集合（不包括 first_sequence）
    used_indices = set()

    # 初始化轮次计数
    round_count = 1

    # 获取所有序列的索引列表（从1开始，因为0是抑制蛋白序列）
    sequence_indices = list(range(1, len(sequences)))

    # 循环处理，直到所有其他序列都被使用或剩余序列长度不足
    while True:
        print(f"\n开始第 {round_count} 轮处理...")

        # 初始化 final_sequence，始终以 first_sequence 开头
        final_sequence = first_sequence
        total_length = len(final_sequence)
        used_indices_in_round = []

        # 计算剩余未使用的序列总长度
        remaining_length = sum(len(sequences[idx]) for idx in sequence_indices if idx not in used_indices)

        # 如果剩余序列长度不足，则停止处理
        if remaining_length < first_sequence_length / 2:
            print("剩余序列长度小于抑制蛋白序列长度的一半，停止处理。")
            break

        # 定义目标和最大长度
        TARGET_LENGTH = 1700
        MAX_LENGTH = 2000

        # 逐个添加未使用的序列，直到达到目标长度或无法添加更多
        for idx in sequence_indices:
            if idx in used_indices:
                continue  # 跳过已使用的序列
            next_sequence = sequences[idx]
            next_length = len(next_sequence)

            if total_length + next_length <= TARGET_LENGTH:
                # 添加序列，不超过目标长度
                final_sequence += next_sequence
                total_length += next_length
                used_indices_in_round.append(idx)
                used_indices.add(idx)
            elif total_length + next_length <= MAX_LENGTH:
                # 添加序列，超过目标但不超过最大长度
                final_sequence += next_sequence
                total_length += next_length
                used_indices_in_round.append(idx)
                used_indices.add(idx)
                break  # 达到或接近目标，停止添加
            else:
                # 添加该序列会超过最大长度，跳过
                continue

        # 如果除了 first_sequence，没有其他序列被添加，则跳过该轮
        if not used_indices_in_round:
            print("没有可添加的未使用序列，停止处理。")
            break

        # 打印当前组合的序列长度和索引
        print(f"组合序列长度为 {total_length} bp，包含序列索引：[0] + {used_indices_in_round}")

        # 开始运行您的代码
        # 输入的 DNA 序列
        sequence = final_sequence

        # 记录开始时间
        start_time = time.time()
        print(sequence)

        # 计算雅可比矩阵
        jac = get_categorical_jacobian(sequence, model, device, alphabet)

        # 计算雅可比矩阵的平均值
        jac_mean = np.mean(jac)

        # 更新 Jacobian 矩阵的指定值
        for i in range(len(sequence)):
            for j in range(4):
                sequence2 = list(sequence)
                sequence2[i] = alphabet.nt[1 + j]  # 替换第 i 个碱基
                sequence2 = ''.join(sequence2)
                codon_start = (i // 3) * 3
                codon_end = codon_start + 3
                codon = sequence2[codon_start:codon_end]  # 提取对应的三联密码子

                # 判断是否为特定的三联密码子
                if codon in ['TTA', 'CTA', 'TCA']:
                    jac[i, j, :, :] = jac_mean  # 更新指定位置的 Jacobian
                    jac[:, :, i, j] = jac_mean

        # 对 Jacobian 矩阵中心化并对称化
        for i in range(4):
            jac -= jac.mean(i, keepdims=True)
        jac = (jac + jac.transpose(2, 3, 0, 1)) / 2

        # 计算接触矩阵
        contacts = get_contacts(jac)

        # 绘制并保存图像
        # 创建 contacts 的副本
        contacts2 = contacts.copy()

        # 将下三角的元素设为零
        lower_triangle_indices = np.tril_indices(len(contacts), -1)  # -1 排除对角线
        contacts2[lower_triangle_indices] = 0

        # 创建图像
        plt.figure(figsize=(15, 15))
        plt.imshow(contacts2, cmap='Greys')
        plt.colorbar()
        plt.clim([0, 0.01])

        # 设置横坐标和纵坐标的刻度每100一次
        plt.xticks(range(0, len(contacts), 100), rotation=90)
        plt.yticks(range(0, len(contacts), 100))

        # 在图中添加水平和垂直线，标记蛋白质序列的区域
        protein_nt_start_in_seq = 0  # 抑制蛋白序列的起始位置为 0
        protein_nt_end_in_seq = len(nucleotide_seq) - 1  # 抑制蛋白序列的结束位置
        plt.axhline(y=protein_nt_start_in_seq, color='green', linewidth=1, label='Protein Region Start')
        plt.axhline(y=protein_nt_end_in_seq, color='green', linewidth=1, label='Protein Region End')
        plt.axvline(x=protein_nt_start_in_seq, color='green', linewidth=1)
        plt.axvline(x=protein_nt_end_in_seq, color='green', linewidth=1)

        # 在图中添加水平线，表示 HTH 结构域的位置
        domain_nt_start_in_seq, domain_nt_end_in_seq = domain_nt_positions
        hth_start_global = protein_nt_start_in_seq + domain_nt_start_in_seq
        hth_end_global = protein_nt_start_in_seq + domain_nt_end_in_seq
        plt.axhline(y=hth_start_global, color='blue', linewidth=1, label='HTH Domain Start')
        plt.axhline(y=hth_end_global, color='blue', linewidth=1, label='HTH Domain End')

        # 计算 HTH 模体与其他序列的互作，排除蛋白质序列区域
        sequence_length = len(sequence)
        protein_indices = set(range(protein_nt_start_in_seq, protein_nt_end_in_seq + 1))
        hth_indices = range(hth_start_global, hth_end_global + 1)

        # 初始化列表存储互作信息
        hth_interactions = []

        for i in hth_indices:
            for j in range(sequence_length):
                # 排除 j 在蛋白质序列区域内的情况
                if j not in protein_indices:
                    interaction_strength = contacts[i, j]
                    if interaction_strength > 0.05:  # 设定阈值
                        hth_interactions.append({
                            'hth_position': i,
                            'other_position': j,
                            'interaction_strength': interaction_strength
                        })

        # 按照互作强度排序
        interactions_list_sorted = sorted(hth_interactions, key=lambda x: x['interaction_strength'], reverse=True)

        # 在图像上叠加 HTH 模体的互作位置
        x_positions = [interaction['other_position'] for interaction in interactions_list_sorted]
        y_positions = [interaction['hth_position'] for interaction in interactions_list_sorted]
        interaction_strengths = [interaction['interaction_strength'] for interaction in interactions_list_sorted]

        # 根据互作强度调整点的大小（可选）
        sizes = [max(1, (s - 0.05) * 100) for s in interaction_strengths]  # 调整系数以适应图像

        plt.scatter(x_positions, y_positions, color='red', s=sizes, label='HTH Interactions')

        # 添加图例
        plt.legend(loc='upper right')

        # 保存图像
        figure_filename = f'contacts_round_{round_count}_with_hth_interactions.png'
        figure_path = os.path.join(output_dir, figure_filename)
        plt.savefig(figure_path)
        plt.close()
        print(f"带有 HTH 互作位置的图像已保存到 {figure_path}")

        # 输出时间
        end_time = time.time()
        elapsed_time = end_time - start_time
        print(f"第 {round_count} 轮处理完成，耗时: {elapsed_time:.2f} 秒")

        # 检查是否所有其他序列都已使用
        if len(used_indices) == len(sequences) - 1:
            print("所有序列均已处理完成。")
            break

        # 增加轮次计数
        round_count += 1

# 清理临时文件（可选）
try:
    os.remove(os.path.join(temp_dir, "protein_seq.fasta"))
    os.remove(os.path.join(temp_dir, "hmmscan_output.txt"))
except OSError as e:
    print(f"Error removing temporary files: {e}")
