from dataprocess import getData
from dataencoder import number_encoder
from models import base
import os
import keras

def write(filepath, seqs, headers):
    f = open(filepath, 'w')
    for seq, header in zip(seqs, headers):
        f.write(header + '\n')
        f.write(seq + '\n')


def predict_type(resultpath, aftertypepath, modelfile2, max_len):
    # filelist = os.listdir(after_merge_path)
    filelist = os.listdir(resultpath)
    filelist.sort(key=lambda x: int(x[x.find('print') + 5: x.find('.')]))
    host, phage = 0, 0
    for i in filelist:
        host_seqs, host_headers, phage_seqs, phage_headers = [], [], [], []
        filepath0 = resultpath + '/' + i
        filepath1 = aftertypepath + '/' + 'Host' + '/' + i
        filepath2 = aftertypepath + '/' + 'Phage' + '/' + i
        seqs, headers = getData(filepath0, True)
        # if seqs == []:
        #     continue
        datas = number_encoder(seqs, max_len)
        keras.backend.clear_session()
        model = base(max_len, 1)
        model.load_weights(modelfile2)

        y_preds = model.predict(datas)
        #print(i)
        #print(y_preds)

        for j in range(len(datas)):
            if y_preds[j] >= 0.5:
                headers[j] += ' type = host'
                host += 1

                host_seqs.append(seqs[j])
                host_headers.append(headers[j])
            else:
                headers[j] += ' type = phage'
                phage += 1

                phage_seqs.append(seqs[j])
                phage_headers.append(headers[j])
        if host_seqs != []:
            write(filepath1, host_seqs, host_headers)
        if phage_seqs != []:
            write(filepath2, phage_seqs, phage_headers)
    #print('host:%d, phage:%d' % (host, phage))
