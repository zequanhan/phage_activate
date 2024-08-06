import numpy as np

def read_fasta(filePath):
    seqs = []
    seq = ''
    with open(filePath) as f:
        for each in f:
            if each.startswith('>'):
                seqs.append(seq)
                seq = ''
            else:
                seq += each.replace('\n', '')
        seqs.append(seq)
      
        seqs.pop(0)
    # seqs = "".join(seqs)
    return seqs



def CGcontent(seq):
    C,G = 0,0
    length = len(seq)
    for temp in seq:
        if temp == 'C':
            C += 1
        elif temp == 'G':
            G += 1
    CG_value = (C+G)/length
    return CG_value


def seq_split(seq):
    seq_split = []
    for i in range(len(seq)):
        temp = seq[i: i+2]
        if len(temp) == 2:
            seq_split.append(temp)
    return seq_split

def self_complementory(seq):
    length = len(seq)//2
    rev_seq = ''
    pre_seq = seq[0: length]
    pos_seq = seq[length:]
    complement = {'C':'G', 'G':'C', 'T':'A', 'A':'T'}
    for temp in pre_seq:
        rev_seq += complement[temp]
    if rev_seq == pos_seq:
        return True
    else:
        return False


def free_energy(seq):
    sseq = seq_split(seq)
    fe = {'AA': -1.00, 'TT':-1.00, 'AT':-0.88, 'TA':-0.58, 'CA':-1.45, 'TG':-1.45, 'GT':-1.44, 'AC':-1.44,
          'CT':-1.28, 'AG':-1.28, 'GA':-1.30, 'TC':-1.30, 'CG':-2.17, 'GC':-2.24, 'GG':-1.84, 'CC':-1.84}
    fe_value = 0
    length = len(sseq)
    # se_GC = 0.98
    # se_AT = 1.03
    se_GC2 = 1.96
    se_AT2 = 2.06
    se_GC_AT = 2.01
 
    sym = 0.43
    if self_complementory(seq):
        fe_value += sym
  
    if seq[0] == 'G' or seq[0] == 'C':
        if seq[length] == 'G' or seq[length] == 'C':
            fe_value += se_GC2
        elif seq[length] == 'A' or seq[length] == 'T':
            fe_value += se_GC_AT
    elif seq[0] == 'A' or seq[0] == 'T':
        if seq[length] == 'A' or seq[length] == 'T':
            fe_value += se_AT2
        elif seq[length] == 'G' or seq[length] == 'C':
            fe_value += se_GC_AT

    for temp in sseq:
        for key in fe:
            if temp == key:
                fe_value += fe[key]

    return fe_value

def atgc_ratio(seq):
    A, T, G, C = 0, 0, 0, 0
    for temp in seq:
        if temp == 'A':
            A += 1
        elif temp == 'T':
            T += 1
        elif temp == 'G':
            G += 1
        elif temp == 'C':
            C += 1
    return (A+T)/(G+C)


def Z_curve(seq):
    A, T, G, C = 0, 0, 0, 0
    for temp in seq:
        if temp == 'A':
            A += 1
        elif temp == 'T':
            T += 1
        elif temp == 'G':
            G += 1
        elif temp == 'C':
            C += 1
    return (A+G)-(C+T), (A+C)-(G+T), (A+T)-(G+C)

def cumulativeSkew(seq):
    A, T, G, C = 0, 0, 0, 0
    for temp in seq:
        if temp == 'A':
            A += 1
        elif temp == 'T':
            T += 1
        elif temp == 'G':
            G += 1
        elif temp == 'C':
            C += 1

    if C + G == 0 and A + T ==0:
        return 10, 10
    if C + G == 0:
        return 10, (A - T) / (A + T)
    if A + T == 0:
        return (G - C) / (C + G), 10

    return (G - C) / (C + G), (A - T) / (A + T)

def com_seq_ATCG(seqs):
    f = []
    for seq in seqs:
        f.append(atgc_ratio(seq))
    return np.array(f)

def com_seq_free_energy(seqs):
    f = []
    for seq in seqs:
        f.append(free_energy(seq))
    return np.array(f)

def com_seq_CGcontent(seqs):
    CG = []
    for seq in seqs:
        CG.append(CGcontent(seq))
    return np.array(CG)

def com_seq_Z_curve(seqs):
    feature = []
    for seq in seqs:
        feature.append(Z_curve(seq))
    return np.array(feature)

def com_seq_cumulativeSkew(seqs):
    feature = []
    for seq in seqs:
        feature.append(cumulativeSkew(seq))
    return np.array(feature)

def com_seq_feature(seqs):
    feature = np.zeros((len(seqs), 7))
    # feature[:, 0] = com_seq_ATCG(seqs)
    feature[:, 0] = com_seq_CGcontent(seqs)
    feature[:, 1:4] = com_seq_Z_curve(seqs)
    feature[:, 4:6] = com_seq_cumulativeSkew(seqs)
    feature[:, 6] = com_seq_free_energy(seqs)
    return feature

