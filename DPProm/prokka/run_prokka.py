import os
import shutil
import subprocess


def run_prokka(prokka_filepath):
    input = prokka_filepath + '/genome'
    output = prokka_filepath + '/prokka_file'
    prokka_file = prokka_filepath + '/prokka.py'
    command = 'python %s -i %s -o %s -t 4' % (prokka_file, input, output)
    subprocess.run(command, shell=True)
    move_file(output + '/prokka', output)
    
   
   
   
   
def move_file(before_move_path, after_move_path):
    fileList = os.listdir(before_move_path)
    for i in fileList:
        filename = 'prokka_results_' + i + '.gbk'
        target_file = before_move_path + '/' + i + '/' + filename
        move_way = after_move_path + '/' + filename
        # print(target_file, move_way)
        shutil.move(target_file, move_way)



def run_prokka_main(prokka_filepath):
    run_prokka(prokka_filepath)
    shutil.rmtree(prokka_filepath + '/prokka_file/prokka')
   
    
