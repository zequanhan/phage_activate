import numpy as np


# input ['ACGT'] return [[1 2 3 4]]
def number_encoder(seqs, max_len):
    id = 'ACGT'
    data_e = []
    for i in range(len(seqs)):
        length = len(seqs[i])
        elemt, st = [], seqs[i]
        for j in st:
            index = id.index(j) + 1
            elemt.append(index)
        if length < max_len:
            elemt += [0] * (max_len - length)
        data_e.append(elemt)
    return np.array(data_e)


