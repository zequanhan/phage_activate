import os
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
        print(path)
        print(os.path.exists(path))
        if os.path.exists(path):
            seqs, header = getData(path, True)
            print('read the %d file with %d sequences' % (i+1, len(seqs)))
            site = get_site(header)
            independ_seqs.append(seqs)
            independ_header.append(site)
            print('reading end')

    return independ_seqs, independ_header



def predict_independ(independpath, modelfile1, max_len, feature_len, resultfile):

    independ_seqs, independ_header = read_independ(independpath)
    print(len(independ_seqs))
    independ_feature, independ_data = [], []
    for i in range(len(independ_seqs)):
        print('dealing with the %d data' % i)
        independ_feature = com_seq_feature(independ_seqs[i])
        independ_data = (number_encoder(independ_seqs[i], max_len))
        print('end of data processing')
        keras.backend.clear_session()
        model = base_feature(max_len, feature_len, 1)
        model.load_weights(filepath=modelfile1)
        print_seqs, print_header, print_y_pred = [], [], []
        y_pred = model.predict([independ_data, independ_feature])
        print(y_pred)
        for j in range(len(y_pred)):
            if y_pred[j] >= 0.5:
                print_seqs.append(independ_seqs[i][j])
                print_header.append(independ_header[i][j])
                print_y_pred.append(y_pred[j])

        if print_seqs != []:
            write_predict(resultfile, print_seqs, print_header, print_y_pred, i)


