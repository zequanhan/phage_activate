import numpy as np
# import matplotlib
# matplotlib.rcParams['font.sans-serif']=['SimHei']
# matplotlib.rcParams['axes.unicode_minus']=False

# read fasta file
def getData(filePath, flag):
    seqs = []
    header = []
    seq = ''
    with open(filePath) as f:
        for each in f:
            if each.startswith('>'):
                header.append(each.replace('\n', ''))
                seqs.append(seq)
                seq = ''
            else:
                seq += each.replace('\n', '')
        seqs.append(seq)
     
        seqs.pop(0)
    if flag:
        return seqs, header
    else:
        return seqs


def posandneg_catch(data):
    # preprocessing label and data
    l = len(data)
    chongfu = 0
    for i in range(l):
        ll = len(data)
        idx = []
        each = data[i]
        j = i + 1
        bo = False
        while j < ll:
            if (data[j] == each):
                
                # label[i] += label[j]
                idx.append(j)
                bo = True
            j += 1
        t = [i] + idx
        if bo:
            #print(t)
            chongfu += 1
           # print(data[t[0]])
            #print(data[t[1]])
        data = np.delete(data, idx, axis=0)
        # label = np.delete(label, idx, axis=0)

        if i == len(data)-1:
            break
    
    # return data, chongfu
    return chongfu

def posorneg_catch(data):
    chongfu = 0
    id = []
    for i in range(len(data)):
        for j in range(i+1, len(data)):
            if (data[i] == data[j]):
                if j not in id:
                    id.append(j)
                    chongfu += 1
                    # print(i, j)
    #print(id)
    data = np.delete(data, id, axis=0)
    return data, chongfu


# get longest length
def maxLength(seqs):
    max_len = 0
    for i in range(len(seqs)):
        st = seqs[i]
        if (len(st) > max_len): max_len = len(st)
    return max_len


# get shortest length
def minLength(seqs):
    min_len = len(seqs[0])
    for i in range(len(seqs)):
        st = seqs[i]
        if (len(st) < min_len): min_len = len(st)
    return min_len


def createLabel(seqs, flag):
    label = []
    length = len(seqs)
    if flag:
        label = [1] * length
    else:
        label = [0] * length
    return np.array(label)








