import os
import numpy as np
from dataprocess import getData


def get_complement(header):
    left_site, right_site = [], []
    for temp in header:
        start_index = temp.find('(')+1
        mid_index = temp.find('..')
        end_index = temp.find(')')
        left_site.append(int(temp[start_index:mid_index]))
        right_site.append(int(temp[mid_index+2:end_index]))
    return left_site, right_site


def get_score(header):
    scores = []
    for temp in header:
        start_index = temp.find('[')+1
        end_index = temp.find(']')
        scores.append(temp[start_index: end_index])
    # print(score)
    return scores


def write_merge(filepath, seqs, left_site, right_site, scores):
    f = open(filepath, 'w')
    for i in range(len(seqs)):
        f.write('>promoter' + str(i) + ' complement (' + str(left_site[i]) + '..' + str(right_site[i]) + ')' + ' score = ' + scores[i] +  '\n')
        f.write(seqs[i] + '\n')
    f.close()
    # f = open(filepath, 'w')
    # f.write('>promoter' + ' complement (' + str(left_site) + '..' + str(right_site) + ')' + ' score = ' + scores + '\n')
    # f.write(seqs + '\n')
    # f.close()

def del_seqs(merge_seqs, left_site, right_site, socres):
    id = []
    for i in range(len(merge_seqs)):
        for j in range(i+1, len(merge_seqs)):
            if merge_seqs[i] == merge_seqs[j]:
                id.append(i)
    #print(id)
    merge_seqs= np.delete(merge_seqs, id, axis=0)
    left_site = np.delete(left_site, id, axis=0)
    right_site = np.delete(right_site, id, axis=0)
    socres = np.delete(socres, id, axis=0)
    return merge_seqs, left_site, right_site, socres

# i = 587, j = 612
def split_seq(seqs, left_site, right_site, i, j):
    new_seqs, ls, rs = [], [], []
    temp_seq = seqs[i]
    #print('temp_seq ', temp_seq)
    l, r = i, 0
    if i == j-1:
        new_seqs.append(seqs[i])
        ls.append(i)
        rs.append(i)
        return new_seqs, ls, rs

    for x in range(i, j):
        y = x + 1
        #print('seqs[y] ', seqs[y])
     
        if right_site[x] - left_site[y] + 1 >= 0 and right_site[y]-left_site[y] + 1 >= 0:
            for k in range(right_site[x] - left_site[y] + 1, right_site[y]-left_site[y] + 1):
                temp_seq = temp_seq + seqs[y][k]
           # print('y = ', y)
            r = y
 
        else:
            new_seqs.append(temp_seq)
            ls.append(l)
            rs.append(r)
            temp_seq = seqs[x+1]
            l = r + 1
            r = y
            if y+1 < len(seqs):
                split_seq(seqs, left_site, right_site, x+1, y+1)
            else:
            #    print('that is the sequence')
                new_seqs.append(temp_seq)
                ls.append(l)
                rs.append(r)

    new_seqs.append(temp_seq)
    ls.append(l)
    rs.append(r)
   # print(new_seqs)
   # print('ls:', ls)
   # print('rs:', rs)
    return new_seqs, ls, rs

def is_find(left_index, right_index, left_site, right_site, merge_seqs, seqs):
    index = []
  
    for i in range(len(merge_seqs)):
        
        for j in range(len(seqs)):
            if merge_seqs[i] == seqs[j] and left_site[left_index[i]] == left_site[j] and right_site[right_index[i]] == right_site[j]:
                index.append(j)
    return index
    

def merge(filepath):

    seqs, header = getData(filepath, True)
    left_site, right_site = get_complement(header)
    scores = get_score(header)
    Merge_seqs, Left_site, Right_site, Scores = [], [], [], []
    i = 0
    while i < len(seqs):
       # print('i=', i)
        j = i + 1
        while j < len(seqs):
            if len(seqs[i]) == len(seqs[j]):
                j += 1
        #        print('1 j=', j)
            else:
         #       print('2 j = ', j)
                break

        merge_seqs, left_index, right_index = split_seq(seqs, left_site, right_site, i, j-1)
        index = is_find(left_index, right_index, left_site, right_site, merge_seqs, seqs)
        if index != []:

            for ii in index:
                Merge_seqs.append(seqs[ii])
                Left_site.append(left_site[ii])
                Right_site.append((right_site[ii]))
                Scores.append(scores[ii])
        i = j
   # print(Merge_seqs, Left_site, Right_site, Scores)
    return Merge_seqs, Left_site, Right_site, Scores


# def max_score(seqs, left_site, right_site, scores):
#     if len(seqs) == 0:
#         return None, None, None, None
#     elif len(seqs) == 1:
#         return seqs[0], left_site[0], right_site[0], scores[0]
#     elif len(seqs) >= 2:
#         index = 0
#         max_score = scores[0]
#         for i in range(len(seqs)):
#             if max_score <= scores[i]:
#                 max_score = scores[i]
#                 index = i
#         return seqs[index], left_site[index], right_site[index], scores[index]


def max_score(scores, index, i):
 
    if scores[index] < scores[i]:
        return i
    else:
        return index

def max_lenseq(seqs, left_site, right_site, scores):
    if len(seqs) == 0:
        return None, None, None, None
    elif len(seqs) == 1:
        return seqs[0], left_site[0], right_site[0], scores[0]
    elif len(seqs) >= 2:
        index = 0
        max_length = len(scores[0])
        for i in range(len(seqs)):
            if max_length < len(scores[i]):
                max_length = len(scores[i])
                index = i
            elif max_length == len(scores[i]):
                index = max_score(scores, index, i)
        return seqs[index], left_site[index], right_site[index], scores[index]
def merge_seqs(resultpath, aftermergepath):
    filelist = os.listdir(resultpath)
    filelist.sort(key=lambda x: int(x[x.find('print') + 5: x.find('.')]))
    for i in filelist:
        filepath1 = resultpath + '/' + i
  
        merge_seqs, Left_site, Right_site, scores = merge(filepath1)
        merge_seqs, Left_site, Right_site, scores = del_seqs(merge_seqs, Left_site, Right_site, scores)

        if merge_seqs.size >0:
            filepath2 = aftermergepath + '/' + i
            write_merge(filepath2, merge_seqs, Left_site, Right_site, scores)
# def merge_seqs(resultpath, aftermergepath):
#     filelist = os.listdir(resultpath)
#     filelist.sort(key=lambda x: int(x[x.find('print') + 5: x.find('.')]))
#     for i in filelist:
#         filepath1 = resultpath + '/' + i
  
#         merge_seqs, Left_site, Right_site, scores = merge(filepath1)
#         merge_seqs, Left_site, Right_site, scores = del_seqs(merge_seqs, Left_site, Right_site, scores)


#         #if merge_seqs != []:2023.12.27修改,因为矩阵不能与列表对比
#         if merge_seqs.size>0:
#             # write_seq, write_left, write_right, write_score = max_lenseq(merge_seqs, Left_site, Right_site, scores)

#             filepath2 = aftermergepath + '/' + i
#             #print(filepath2)
#             # if write_seq != None:
#             write_merge(filepath2, merge_seqs, Left_site, Right_site, scores)
