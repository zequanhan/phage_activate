import os
import shutil
import subprocess


def runCDHIT(aftermergepath, cdhitpath, cdhit_seqs_path):
    fileList = os.listdir(aftermergepath)
    fileList.sort(key=lambda x: int(x[x.find('print') + 5: x.find('.')]))
    for i in fileList:
        # cdhit_seqs, cdhit_headers = [], []
        filepath1 = aftermergepath + '/' + i
        # seqs, header = getData(filepath1, True)
        filepath2 = cdhitpath + '/' + i
        command = "cd /data1/tools/cdhit-master; cd-hit -i %s -o %s -c 0.7 -n 5" % (filepath1, filepath2)
        subprocess.run(command, shell=True)
        filepath3 = cdhit_seqs_path + '/' + i
        shutil.move(filepath2, filepath3)

    # shutil.rmtree('/home/wangc/Desktop/predict_promoter4/cdhit')




