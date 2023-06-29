import os
os.environ["CUDA_DEVICE_ORDER"]="PCI_BUS_ID"
os.environ["CUDA_VISIBLE_DEVICES"]="1"
# os.environ["CUDA_VISIBLE_DEVICES"]="-1" # debug on cpu
import tensorflow as tf
from Fingerprint.mlp import MLP
from scipy.sparse import csr_matrix
from scipy.sparse import load_npz
import numpy as np
import pickle

def train():
    spec = load_npz('DeepEI/data/paper/nist_train_spec.npz')
    fps = load_npz('DeepEI/data/paper/nist_train_fingerprint.npz')
    spec_test = load_npz('DeepEI/data/paper/nist_test_spec.npz')
    fps_test = load_npz('DeepEI/data/paper/nist_test_fingerprint.npz')

    fp_list = map(int, open('DeepEI/data/paper/fpkeep_01_09.txt', 'r').read().strip().split())
    spec = spec.todense()
    spec_test = spec_test.todense()
    fps = csr_matrix(fps)
    fps_test = csr_matrix(fps_test)
    # training the model one-by-one, it will be time-consuming
    # the fp positions to train are fixed, so read them from a list
    for i, fp_idx in enumerate(fp_list):
        y = fps[:,i].todense()
        y_test = fps_test[:,i].todense()
        y = np.squeeze(np.asarray(y))
        y_test = np.squeeze(np.asarray(y_test))
        
        # check: 0.1 < bias < 0.9
        # the fp positions to train are known, so we don't filter here
        # fr = np.sum(y) / len(y)
        # if (fr < 0.1) or (fr > 0.9):
        #     continue
        Y = np.vstack((y, (1-y))).transpose()
        Y_test = np.vstack((y_test, (1-y_test))).transpose()
        
        # for write the evaluation results 
        mlp_result = open('saved/model/paper/nist_only/mlp_result.txt', 'a+')
        mlp = MLP(spec, Y, spec_test, Y_test)
        history = mlp.train()
        mlp_res = mlp.test()
        mlp_result.write("\t".join([str(fp_idx)] + [str(j) for j in mlp_res]))
        mlp_result.write("\n")
        mlp.save('saved/model/paper/nist_only/weight/{}.h5'.format(fp_idx), 'saved/model/paper/nist_only/model.json')
        with open('saved/model/paper/nist_only/history/{}.history'.format(fp_idx), 'wb') as file_pi:
            pickle.dump(history.history, file_pi)

def predict():
    pass

def config_gpu_mem(mem_size=1024):
    gpus = tf.config.list_physical_devices('GPU')
    if gpus:
        try:
            tf.config.set_logical_device_configuration(
                gpus[0],
                [tf.config.LogicalDeviceConfiguration(memory_limit=mem_size)])
            logical_gpus = tf.config.list_logical_devices('GPU')
            print(len(gpus), "Physical GPUs,", len(logical_gpus), "Logical GPUs")
        except RuntimeError as e:
            print(e)

if __name__ == "__main__":
    config_gpu_mem()
    train()