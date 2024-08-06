
import os




def remove_file(fileList):
    # 清除文件
    # shutil.rmtree('independ_test/')
    # shutil.rmtree('result/')
    # shutil.rmtree('after_merge_result/Host/')
    # shutil.rmtree('after_merge_result/Phage/')
    print('clearing')
    for i in range(len(fileList)):
        independ_test = fileList[i]
        independ_test_fileList = os.listdir(independ_test)
        for name in independ_test_fileList:
            os.remove(fileList[i] + '/' + name)




# if __name__ == '__main__':
#     fileList = ['independ_test', 'result', 'after_merge_result/Host', 'after_merge_result/Phage']
#     remove_file(fileList)