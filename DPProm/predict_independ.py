import os
from joblib import load

from dataprocess import getData
from dataencoder import number_encoder
from feature import com_seq_feature
from models import base_feature
import keras
os.environ["TF_CPP_MIN_LOG_LEVEL"] = '3'


def get_site(header):
    site = []
    for j in header:
        site.append(j[j.index('complement') - 1:])
    return site



def write_predict(resultfile, seqs, header, y_pred, index):
    f = open(resultfile + str(index) + '.txt', 'w')
    for i in range(len(seqs)):
        f.write('>promoter' + str(i) + header[i] + ' score = ' + str(y_pred[i]) + '\n')
        f.write(seqs[i] + '\n')
    f.close()



def read_independ(filepath):
    independ_seqs, independ_header = [], []
    filelist = os.listdir(filepath)

   
    filelist.sort(key=lambda x: int(x[x.find('data') + 4: x.find('.')]))

    for i in range(len(filelist)):
        seqs, header = [], []
        path = filepath + '/' + filelist[i]
        #print(path)
        #print(os.path.exists(path))
        if os.path.exists(path):
            seqs, header = getData(path, True)
            #print('read the %d file with %d sequences' % (i+1, len(seqs)))
            site = get_site(header)
            independ_seqs.append(seqs)
            independ_header.append(site)
            #print('reading end')

    return independ_seqs, independ_header



from joblib import load
import numpy as np

# def predict_independ(independpath, modelfile1, max_len, feature_len, resultfile):
#     # 读取独立数据集
#     independ_seqs, independ_header = read_independ(independpath)
#     print('Start data processing and prediction')
#     loaded_model = load(modelfile1)  # 加载模型

#     # 初始化用于存储特征的列表
#     all_independ_features = []
#     all_independ_data = []
    
#     # 批量处理序列
#     for seq in independ_seqs:
#         # 提取每个序列的特征和编码
#         seq_feature = com_seq_feature(seq)  # 假设返回shape为(feature_len,)
#         seq_data = number_encoder(seq, max_len)  # 假设返回shape为(max_len,)

#         # 存储特征
#         all_independ_features.append(seq_feature)
#         all_independ_data.append(seq_data)

#     # 将列表转换为NumPy数组，并调整形状
#     all_independ_features = np.array(all_independ_features)  # shape为(n_samples, feature_len)
#     all_independ_data = np.array(all_independ_data)  # shape为(n_samples, max_len)

#     # 确保数据形状正确
#     all_independ_data = np.reshape(all_independ_data, (len(independ_seqs), -1))
#     all_independ_features = np.reshape(all_independ_features, (len(independ_seqs), -1))

#     # 合并特征用于预测
#     X_independ = np.hstack([all_independ_data, all_independ_features])

#     print('End of data processing')
#     y_pred = loaded_model.predict_proba(X_independ)[:, 1]  # 获取属于正类的概率
#     print(f'Prediction: {y_pred}')

 
#     # 写入预测结果到文件
#     if print_seqs:
#         write_predict(resultfile, print_seqs, print_header, print_y_pred)

def predict_independ(independpath, modelfile1, max_len, feature_len, resultfile):
    # 读取独立数据集
    independ_seqs, independ_header = read_independ(independpath)
    print('Start data processing and prediction')
    loaded_model = load(modelfile1)  # 加载模型

    for idx, seq in enumerate(independ_seqs):
        # 提取每个序列集的特征和编码
        seq_feature = com_seq_feature(seq)  # 假设返回形状为 (n_samples, feature_len)
        seq_data = number_encoder(seq, max_len)  # 假设返回形状为 (n_samples, max_len)
        
        # 合并特征和数据用于预测
        X_batch = np.hstack([seq_data, seq_feature])
        
        # 对这批数据进行预测
        y_pred = loaded_model.predict_proba(X_batch)[:, 1]  # 获取属于正类的概率
        
        # 根据预测结果处理每个序列
        print_seqs, print_header, print_y_pred = [], [], []
        for i, pred_prob in enumerate(y_pred):
            if pred_prob >= 0.5:
                print_seqs.append(seq[i])  # 假设seq存储的是原始序列字符串
                print_header.append(independ_header[idx][i])  # 适当调整以匹配您的header结构
                print_y_pred.append(pred_prob)

        # 根据预测结果将这批序列写入文件
        if print_seqs:
#             print(idx)
            write_predict(resultfile, print_seqs, print_header, print_y_pred,idx)

    print('End of data processing and prediction')
