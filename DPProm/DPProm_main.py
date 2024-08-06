import os
import warnings

import sys
# sys.path.append('/data4/wangchen/DPProm_website/DPProm/')
#sys.path.append('/data1/WWW/flask_website/DPProm/DPProm/')

from merge_seqs import merge_seqs
from prokka.cut_genome import cut_genome_seqs
from prokka.run_prokka import run_prokka_main
from cdhit import runCDHIT
from type import predict_type
from remove_file import remove_file
from read_and_write import show_allseqs
from predict_independ import predict_independ
warnings.filterwarnings("ignore", category=Warning)
os.environ["TF_CPP_MIN_LOG_LEVEL"] = '3'



# def get_filename(filepath):
    # filelist = os.listdir(filepath)
    # filename = filelist[0]
    # return filename


# if __name__ == '__main__':
def ensure_directory_exists(directory_path):
    """确保给定的目录存在。如果不存在，则创建它。"""
    if not os.path.exists(directory_path):
        os.makedirs(directory_path)
        
        
def genome_predict(path, gbk_file_path, genome_file_path):
    
 
    # path = '/data4/wangchen/website/DPProm/'
    #path = '/home/hanzequan/DPProm/DPProm/'
    independpath = path + 'independ_test'
    resultpath = path + 'result'
    aftermergepath = path + 'after_merge_result'
    aftertypepath = path + 'after_type_result'

    # 确保所有需要的文件夹都存在2023.10.16
    paths_to_check = [independpath, resultpath, aftermergepath, aftertypepath]
    for p in paths_to_check:
        ensure_directory_exists(p)



    posfile = path + 'data/after_catch_promoters.txt'
    negfile = path + 'data/non_promoters.txt'
   # modelfile1 = path + 'model/model1.h5'
    current_dir = os.path.dirname(os.path.abspath(__file__))

# 构建 modelfile1 的相对路径
    modelfile1 = os.path.join(current_dir, 'best_random_forest_model.joblib')
    #modelfile1 = '/home/hanzequan/test_bectiral/rf_model/best_random_forest_model.joblib' 
    resultfile = path + 'result/print'
    aftermergefile = path + 'after_merge_result/print'
    

    
    
    
    
    prokka_filepath = path + 'prokka'
    gbk_file = gbk_file_path
    genome_file = genome_file_path
#     gbk_file = path + 'prokka/prokka_file/prokka_results_L5.gbk'#可修改
#     genome_file = path + 'prokka/genome/L5.fasta'#可修改
    independ_test_seqs_file = path + 'independ_test/independ_test_data'
    # show_file = aftertypepath + '/show.txt'

    modelfile2 =  path + 'model/model2.h5'
    hostfile =  path + 'data/host.fasta'
    phagefile =  path + 'data/phage.fasta'

    hosttypefile = path + 'after_type_result/Host'
    phagetypefile = path + 'after_type_result/Phage'
    

    cdhit = path + 'cdhit'
    cdhit_seqs = path + 'cdhit_seqs'
    

    fileList = [independpath, resultpath, aftermergepath, hosttypefile, phagetypefile, cdhit, cdhit_seqs]
    #确保所有需要的文件夹都存在2023.10.16
    for p in fileList:
        ensure_directory_exists(p)
    remove_file(fileList)
    # os.unlink(show_file)
  
    run_prokka_main(prokka_filepath)
 
    cut_genome_seqs(gbk_file, genome_file, independ_test_seqs_file)

 
    predict_independ(independpath, modelfile1, 99, 7, resultfile)
   
    merge_seqs(resultpath, aftermergepath)
    runCDHIT(aftermergepath, cdhit, cdhit_seqs)
    
    predict_type(cdhit_seqs, aftertypepath, modelfile2, 99)
    seqs, headers = show_allseqs(aftertypepath)
    for i in range(len(headers)):
        headers[i] = 'promoter ' + str(i) + ' ' + headers[i][headers[i].find('complement'):].replace('>', '')
    #os.unlink(genome_file)
    # print(headers)
    return seqs, headers
