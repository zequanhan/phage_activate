import os
from dataprocess import getData


def read_after_merge_result(filepath):
    fileList = os.listdir(filepath)
    # fileList[0] = 'Host'
    host_filepath = filepath + '/' + fileList[0] + '/'
    phage_filepath = filepath + '/' + fileList[1] + '/'
    hostList = os.listdir(host_filepath)
    phageList = os.listdir(phage_filepath)

    host_seqs, host_headers = [], []
    phage_seqs, phage_headers = [], []
    for i in hostList:
        hostfile = host_filepath + i
        seqs, headers = getData(hostfile, True)
        host_seqs += seqs
        host_headers += headers
    for j in phageList:
        phagefile = phage_filepath + j
        seqs, headers = getData(phagefile, True)
        phage_seqs += seqs
        phage_headers += headers

    #print(len(host_seqs+phage_seqs), len(host_headers+phage_headers))
    return host_seqs+phage_seqs, host_headers+phage_headers

def write_seq(filepath, seqs, headers):
    f = open(filepath, 'w')
    for i in range(len(seqs)):
        f.write('>promoter' + str(i) + ' ' + headers[i][headers[i].find('complement'):] + '\n')
        f.write(seqs[i] + '\n')
    f.close()

# filepath1 = 'after_type_result'
# filepath2 = ''
def show_allseqs(filepath1):
    seqs, headers = read_after_merge_result(filepath1)
    #print(headers)
    # write_seq(filepath2, seqs, headers)
    return seqs, headers


if __name__ == '__main__':
    path = '/home/wangc/Desktop/website/DPProm/'
    aftertypepath = path + 'after_type_result'
    seqs, headers = show_allseqs(aftertypepath)
    #print(seqs)

