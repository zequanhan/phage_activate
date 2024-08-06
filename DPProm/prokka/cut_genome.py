import random
# cut_genome-->(feature_genome_main)-->predict_independ-->merge-->independ_score

#def gbk_reader(gbk_path):
#    position = []
#    with open(gbk_path, 'r', encoding='utf-8',errors='ignore') as f:
#        for line in f:
#            if ' CDS ' in line:
#                line = line.replace('\t', '').replace('\n', '').replace(' ', '').replace('CDS', '').replace(
 #                   'complement(', '').replace(')', '')
 #               head, _, tail = line.partition('..')
 #               print(head, tail)
 #              position.append([int(head), int(tail)])
 #   return position
def gbk_reader(gbk_path):
    position = []
    with open(gbk_path, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            if ' CDS ' in line:
                # 移除多余的字符
                line = line.replace('\t', '').replace('\n', '').replace(' ', '').replace('CDS', '').replace('complement(', '').replace(')', '')
                
                # 检查是否包含 '..'，表明是一个位置范围
                if '..' in line:
                    head, _, tail = line.partition('..')
                    
                    # 仅当头尾都是数字时添加位置
                    if head.isdigit() and tail.isdigit():
                        position.append([int(head), int(tail)])
                    else:
                        print(f"Skipped complex position: {line}")
    return position

def get_data(filePath):
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
    return seqs



def non_coding_area(genome, position):
    length = len(genome)
    new_position = []
    for i in range(len(position) - 1):
        if position[i][1] < position[i + 1][0]:
            lenth_genome = 50
            p = [position[i][1] -lenth_genome, position[i + 1][0] +lenth_genome]#修改编码区向基因内延伸长度
           
            #p = [position[i][1], position[i + 1][0]]
            #print('提取基因间间隔区:',p)
            new_position.append(p)
    if position[0][0] > 1:
        new_position.insert(0, [1, position[0][0] - 1])
    if position[-1][1] < length:
        new_position.insert(-1, [position[-1][1] + 1, length])
    return new_position



def sort(position):

    for i in range(len(position)-1):
        for j in range(i+1, len(position)):
            if position[i][0] > position[j][0]:
                temp_left, temp_right = position[j][0], position[j][1]
                position[j][0], position[j][1] = position[i][0], position[i][1]
                position[i][0], position[i][1] = temp_left, temp_right

    return position
    

def cut_genome(genome, position):
    new_position = non_coding_area(genome, position)
    #print('before sorting')
    #print(new_position)
 
    new_position = sort(new_position)
    #print('after sorting')
    #print(new_position)
    #print('genome length', len(genome))
    #print(new_position[-1])

    if new_position[-1][1] == len(genome):
        pass#修改能够取到基因后的非编码区
        #new_position.pop()
    #print(new_position)
    nca = []
    for p in new_position:
        seq = genome[p[0] - 1: p[1]]
        nca.append(seq)
    return nca, new_position
    
    
    
# def cut_genome(genome, position):
    # new_position = non_coding_area(genome, position)
    # nca = []
    # for p in new_position:
        # seq = genome[p[0] - 1: p[1]]
        # nca.append(seq)
    # return nca, new_position



def getKmers(sequence, position, w):
    all_seqs, all_poss = [], []
    #print('step:', 1)
    for x in range(0, len(sequence)-w+1, 1):
        seqs = sequence[x:x+w].upper()
        pos = [position[0]+x, position[0]+w+x-1]
        all_seqs.append(seqs)
        all_poss.append(pos)

    return all_seqs, all_poss



def getseq(nca, new_position, independ_test_seqs_path, min_length=21, max_length=99):

    for i in range(len(nca)):
        s, p = [], []
       # print("the %d sequence, the length is %d" % (i, len(nca[i])))
        for window_size in range(min_length, max_length+1):

            if len(nca[i]) >= window_size:
                #print('window size:', window_size)
                seqs, poss = getKmers(nca[i], new_position[i], window_size)
                s = s + seqs
                p = p + poss

        if s != [] and p != []:
            write(independ_test_seqs_path, s, p, i)


def write(independ_test_seqs_path, seqs, position, num):
    #print("write to file")
    f = open(independ_test_seqs_path + str(num) + '.fasta', 'a')
    for i in range(len(seqs)):
        f.write('>promoter ' + str(i) + ' complement(' + str(position[i][0]) + '..' + str(position[i][1]) + ')' + '\n')
        f.write(seqs[i])
        f.write('\n')
    f.close()

def isATCG(seqs):
    flag = 1
    for seq in seqs:
        if seq not in ['A', 'C', 'G', 'T']:
            flag = 0
    return flag

def replace(seqs):
    new_seqs = []
    for seq in seqs:
        if isATCG(seq) == 0:
            #print('Character substitution in progress')
            seq = seq.upper()
            new_seq = ''
            for s in seq:
                new_s = s.replace('W', random.choice('AT')).replace('S', random.choice('CG')).replace('R', random.choice('AG')) \
                    .replace('Y', random.choice('CT')).replace('K', random.choice('GT')).replace('M',random.choice('AC')) \
                    .replace('B', random.choice('CGT')).replace('D', random.choice('AGT')) \
                    .replace('H', random.choice('ACT')).replace('V', random.choice('ACG')) \
                    .replace('N', random.choice('ACGT'))
                new_seq += new_s
            new_seqs.append(new_seq)
        else:
            # print('No character substitution required')
            new_seqs.append(seq)
    return new_seqs

def cut_genome_seqs(gbk_path, genome_path, independ_test_seqs_path):
    position = gbk_reader(gbk_path)
    genomes = get_data(genome_path)
    nca, nca_position = cut_genome(genomes[0], position)
    new_nca = replace(nca)
    # print(nca)
    # print(len(nca))
    # print(nca_position)
    getseq(new_nca, nca_position, independ_test_seqs_path)


