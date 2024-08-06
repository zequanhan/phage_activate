import re
import pandas as pd
import numpy as np
from Bio.Seq import Seq
from Bio import motifs
import logomaker
from scipy.stats import pearsonr
from scipy.spatial.distance import cosine

def build_pwm(sequences):
    max_length = max(len(seq) for seq in sequences)
    counts = {base: np.zeros(max_length) for base in 'ACGT'}
    
    for seq in sequences:
        padded_seq = seq + '-' * (max_length - len(seq))
        for idx, char in enumerate(padded_seq):
            if char in counts:
                counts[char][idx] += 1
    
    total_counts = np.zeros(max_length)
    for base in 'ACGT':
        total_counts += counts[base]
    
    pwm = {base: np.divide(counts[base], total_counts, out=np.zeros_like(counts[base]), where=total_counts!=0) for base in 'ACGT'}
    
    return pwm

def process_sequences(file_path, target_accession):
    with open(file_path, 'r') as file:
        content = file.read()

    cleaned_content = content.replace('..', ',')
    entries = cleaned_content.strip().split('>')[1:]

    header_pattern = re.compile(r'RefSeq:\s*([^,]+),.*?\((\d+)\s*[,.]\s*(\d+)\)')
    data = []

    for entry in entries:
        lines = entry.strip().split('\n')
        header = lines[0]
        sequence = ''.join(lines[1:]).strip()
        match = header_pattern.search(header)
        if match:
            accession = match.group(1).strip()
            start = match.group(2)
            end = match.group(3)
            data.append({
                'accession': accession,
                'start': start,
                'end': end,
                'sequence': sequence
            })

    df = pd.DataFrame(data)
    if target_accession not in df['accession'].values:
        print(f"Accession {target_accession} not found in the dataset.")
        return None
    
    sequences = df[df['accession'] == target_accession]['sequence'].tolist()
    min_length = min(len(seq) for seq in sequences)
    sequences_truncated = [seq[:min_length] for seq in sequences]
    sequences_truncated = [Seq(seq) for seq in sequences_truncated]
    m = motifs.create(sequences_truncated)
    pwm = m.counts.normalize(pseudocounts={'A': 0, 'C': 0, 'G': 0, 'T': 0})
    pwm_df_true = pd.DataFrame({nucleotide: pwm[nucleotide] for nucleotide in 'ACGT'})

    return pwm_df_true

def predict_pwm(root_directory, target_accession, sequence_count_occurrences):
    directories = find_directories_with_string(root_directory, target_accession)
    all_motifs_df = pd.DataFrame()

    for directory in directories:
        try:
            meme_df = build_motif_matrices(directory, sequence_count_occurrences)
            all_motifs_df = pd.concat([all_motifs_df, meme_df], ignore_index=True)
        except FileNotFoundError as e:
            print(e)
    
    if all_motifs_df.empty:
        print(f"No motif data found for accession {target_accession}.")
        return None

    comparison_results_df = run_comparisons_on_motifs(all_motifs_df)
    final_df = comparison_results_df[comparison_results_df['State Change'] == comparison_results_df['State Change'].max()]
    best_motif = all_motifs_df[all_motifs_df['Motif'] == merge_sequences_based_on_identity(final_df)[merge_sequences_based_on_identity(final_df).index == 0]['Motif'][0]]
    sequences = best_motif['Sequence'].tolist()

    min_length = min(len(seq) for seq in sequences)
    sequences_truncated = [seq[:min_length] for seq in sequences]
    sequences_truncated = [Seq(seq) for seq in sequences_truncated]
    m = motifs.create(sequences_truncated)
    pwm = m.counts.normalize(pseudocounts={'A': 0, 'C': 0, 'G': 0, 'T': 0})
    pwm_df_pred = pd.DataFrame({nucleotide: pwm[nucleotide] for nucleotide in 'ACGT'})

    return pwm_df_pred

def kl_divergence(p, q):
    """计算KL散度，假设p和q是同样大小的numpy数组，且q中没有零值。"""
    return np.sum(np.where(p != 0, p * np.log(p / q), 0))

def compare_pwms(pwm_df_true, pwm_df_pred):
    # 确保两个矩阵在行数上长度一致
    max_rows = max(pwm_df_true.shape[0], pwm_df_pred.shape[0])
    
    # 在较短的矩阵中添加新的行，填充值为0.25
    if pwm_df_true.shape[0] < max_rows:
        additional_rows = pd.DataFrame(0.25, index=range(max_rows - pwm_df_true.shape[0]), columns=pwm_df_true.columns)
        pwm_df_true = pd.concat([pwm_df_true, additional_rows], ignore_index=True)
        
    if pwm_df_pred.shape[0] < max_rows:
        additional_rows = pd.DataFrame(0.25, index=range(max_rows - pwm_df_pred.shape[0]), columns=pwm_df_pred.columns)
        pwm_df_pred = pd.concat([pwm_df_pred, additional_rows], ignore_index=True)

    # 转换为numpy数组
    pwm_true = pwm_df_true.values
    pwm_pred = pwm_df_pred.values

    # 展平矩阵以计算皮尔森相关系数和余弦相似度
    pwm_true_flat = pwm_true.flatten()
    pwm_pred_flat = pwm_pred.flatten()

    # 确保数组长度一致
    min_length = min(len(pwm_true_flat), len(pwm_pred_flat))
    pwm_true_flat = pwm_true_flat[:min_length]
    pwm_pred_flat = pwm_pred_flat[:min_length]

    # 计算皮尔森相关系数
    correlation, _ = pearsonr(pwm_true_flat, pwm_pred_flat)

    # 计算余弦相似度
    cos_sim = 1 - cosine(pwm_true_flat, pwm_pred_flat)

    # 添加小数以避免除以零错误
    epsilon = 1e-10
    pwm_true += epsilon
    pwm_pred += epsilon

    # 计算总的KL散度
    total_kl_div = sum([kl_divergence(pwm_true[i, :], pwm_pred[i, :]) for i in range(max_rows)])

    # 返回分数矩阵
    scores = {
        'Pearson Correlation': correlation,
        'Cosine Similarity': cos_sim,
        'Total KL Divergence': total_kl_div
    }
    return scores


def analyze_accession(file_path, root_directory, target_accession):
    sequence_count_occurrences = {}
    pwm_df_true = process_sequences(file_path, target_accession)
    if pwm_df_true is None:
        return None
    
    pwm_df_pred = predict_pwm(root_directory, target_accession, sequence_count_occurrences)
    if pwm_df_pred is None:
        return None
    
    scores = compare_pwms(pwm_df_true, pwm_df_pred)
    return target_accession, scores

def main(file_path, root_directory, accessions):
    results = []
    for accession in accessions:
        result = analyze_accession(file_path, root_directory, accession)
        if result:
            target_accession, scores = result
            results.append([target_accession, scores['Pearson Correlation'], scores['Cosine Similarity'], scores['Total KL Divergence']])
        else:
            print(f"Analysis for accession {accession} could not be completed.")
    
    results_df = pd.DataFrame(results, columns=['Accession', 'Pearson Correlation', 'Cosine Similarity', 'Total KL Divergence'])
    return results_df

# Example usage
file_path = '/home/hanzequan/test_bectiral/operator_recongize/operator_phisit.csv'
root_directory = '/home/hanzequan/test_bectiral/operator_recongize/all_tree'
accessions = df['accession'].unique()

results_df = main(file_path, root_directory, accessions)
print(results_df)

