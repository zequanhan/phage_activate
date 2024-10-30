import torch
from Bio import SeqIO
import random
import numpy as np
import matplotlib.pyplot as plt

# load the model
device = 'cuda' if torch.cuda.is_available() else 'cpu'

model_path = "/home/hanzequan/megaDNA/megaDNA_phage_145M.pt"  # model name
model = torch.load(model_path, map_location=torch.device(device))
model.eval()  # Set the model to evaluation mode


class Alphabet:
    def __init__(self):
        self.nt = ['**', 'A', 'T', 'C', 'G', '#']  # Alphabet
        self.ln = 4  # 核苷酸数量
        self.map_encode = OrderedDict((n, i) for i, n in enumerate(self.nt))
        self.map_decode = OrderedDict((v, k) for k, v in self.map_encode.items())
        self.start = [0]
        self.end = [5]
        self.codon = list(self.map_encode.values())[1:-1]

    def encode(self, sequence):
        # 将序列编码为数字
        x = [self.map_encode[s] if s in self.map_encode.keys() else 1 for s in sequence]
        return self.start + x + self.end

    def decode(self, x):
        # 将数字解码为序列
        return [self.map_decode[i] for i in x]


# 计算雅可比矩阵的函数
def get_categorical_jacobian(seq, model, device, alphabet):
    with torch.no_grad():
        x, ln = torch.tensor(alphabet.encode(seq)).unsqueeze(0).to(device), len(seq)
        f = lambda k: model(k, return_value='logits')[..., 1:-1, 1:5].cpu().numpy()

        # 计算原始的 f(x)
        fx = f(x)[0]
        fx_h = np.zeros((ln, alphabet.ln, ln, alphabet.ln))

        # 克隆多个 x 并逐一修改
        xs = torch.tile(x, [alphabet.ln, 1])
        for n in range(ln):
            x_h = torch.clone(xs)
            x_h[:, n + 1] = torch.tensor(alphabet.codon)
            fx_h[n] = f(x_h)

        return fx_h - fx  # 返回 f(x+h) - f(x)


# 去均值或奇异值分解修正函数
def do_apc(x, rm=1):
    x = np.copy(x)
    if rm == 0:
        return x
    elif rm == 1:
        # 使用均值修正
        a1 = x.sum(0, keepdims=True)
        a2 = x.sum(1, keepdims=True)
        y = x - (a1 * a2) / x.sum()
    else:
        # 使用奇异值分解修正
        u, s, v = np.linalg.svd(x)
        y = s[rm:] * u[:, rm:] @ v[rm:, :]

    np.fill_diagonal(y, 0)  # 对角线置零
    return y


# 将雅可比矩阵转换为接触矩阵
def get_contacts(jacobian, symm=True, center=True, rm=1):
    j = jacobian.copy()

    # 中心化
    if center:
        for i in range(4):
            j -= j.mean(i, keepdims=True)

    # 计算 Frobenius 范数并修正
    j_fn = np.sqrt(np.square(j).sum((1, 3)))
    np.fill_diagonal(j_fn, 0)
    j_fn_corrected = do_apc(j_fn, rm=rm)

    # 对称化
    if symm:
        j_fn_corrected = (j_fn_corrected + j_fn_corrected.T) / 2

    return j_fn_corrected


# 将矩阵缩小为 N x N
def shrink_matrix(matrix, N):
    assert matrix.shape[0] == 3 * N and matrix.shape[1] == 3 * N, "Matrix must be of size 3N x 3N"
    shrunk_matrix = np.zeros((N, N))

    # 计算每个 3x3 块的均值
    for i in range(N):
        for j in range(N):
            shrunk_matrix[i, j] = np.mean(matrix[3 * i:3 * i + 3, 3 * j:3 * j + 3])

    return shrunk_matrix
