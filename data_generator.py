import tensorflow as tf
import h5py
import numpy as np
import math


chunk_len = 1000


def load_hdf_to_predict(file_name, keys):
    """
    add: un-cut train val test dataset and use the most simple list to load the data
    :param file_name: file_name
    :param keys: features labels
    :param type: tensor_dtypes = [np.float32, np.int64]
    :return:  a list of dataset for traini
    """
    with h5py.File(file_name, 'r') as f:
        features = []
        read_ids = []
        for ids in range(0, len(f.get('features'))):
            features.append(np.array(f[keys[0]][ids]))
            read_ids.append(np.array(f[keys[1]][ids]))
    return features, read_ids

def load_hdf_to_train(file_name, keys, type):
    """
    1.0: un-cut train val test dataset and use the most simple list to load the data
    :param file_name: file_name
    :param keys: features labels
    :param type: tensor_dtypes = [np.float32, np.int64]
    :return:  a list of dataset for traini
    """
    with h5py.File(file_name, 'r') as f:
        features = []
        labels = []
        for ids in range(0, len(f.get('features'))):
            features.append(np.array(f[keys[0]][ids], dtype=type[0]))
            labels.append(np.array(f[keys[1]][ids], dtype=type[1]))
    return features, labels


def load_hdf_with_percent(file_name, keys,  mode):
    """
    3.0: load the data per chunk, and flatten per list to ndarray
    :param file_name: the name of hdf file created by the long read
    :param keys: tensor_keys = ["features", "labels"]
    :param mode: train, val, test
    :return: a list include all the chunk and it is proved not suitable for training
    """
    with h5py.File(file_name, 'r') as f:
        features = []
        labels = []
        for ids in range(0, len(f.get('features'))): #
            if mode == "train":
                features.append(f[keys[0]][ids][0:int(chunk_len*0.8)])
                labels.append(f[keys[1]][ids][0:int(chunk_len*0.8)])

            if mode == "val":
                features.append(f[keys[0]][ids][int(chunk_len * 0.8): int(chunk_len * 0.9)])
                labels.append(f[keys[1]][ids][int(chunk_len * 0.8): int(chunk_len * 0.9)])

            if mode == "test":
                features.append(f[keys[0]][ids][int(chunk_len * 0.9):])
                labels.append(f[keys[1]][ids][int(chunk_len * 0.9):])

    return features, labels


class load_hdf_iter_percent(tf.keras.utils.Sequence):
    """
    4.0: load the data per batch, include batch_size chunks(32*800*11)
    :return a iter batch_size data for training, val and test
    """
    def __init__(self, hdf_files, keys, mode, batch_size=256):
        """
        initial the parameters
        :param hdf_files:
        :param keys:
        :param mode: train val test
        :param batch_size:
        """
        self.file = h5py.File(hdf_files,'r')
        self.keys = keys
        self.mode = mode
        self.batch_size = batch_size

    def __len__(self):
        """
        即数据第一维的数目/batch_size向上取整 (109/32)=4  至少得是1，为0会无法进去getitem函数
        且这里range的范围是__getitem__(self, idx) idx的循环范围
        :return:
        """
        return math.ceil(self.file['features'].shape[0] / self.batch_size) # 向上取整输入的是batch的数目

    def __getitem__(self, idx):
        """生成一个batch的数据"""

        while 1:
            f = self.file
            features = []
            labels = []
            if self.mode == "train":

                features = f[self.keys[0]][idx*self.batch_size:(idx+1)*self.batch_size, 0:int(chunk_len * 0.8)]
                labels = f[self.keys[1]][idx*self.batch_size:(idx+1)*self.batch_size, 0:int(chunk_len * 0.8)]

            if self.mode == "val":
                features = f[self.keys[0]][idx*self.batch_size:(idx+1)*self.batch_size, int(chunk_len * 0.8): int(chunk_len * 0.9)]
                labels = f[self.keys[1]][idx*self.batch_size:(idx+1)*self.batch_size, int(chunk_len * 0.8): int(chunk_len * 0.9)]

            if self.mode == "test":
                features = f[self.keys[0]][idx*self.batch_size:(idx+1)*self.batch_size, int(chunk_len * 0.9):chunk_len]
                labels = f[self.keys[1]][idx*self.batch_size:(idx+1)*self.batch_size, int(chunk_len * 0.9):chunk_len]

            return tf.convert_to_tensor(features), tf.convert_to_tensor(labels)  # 偶尔np.array()送进去会报错

class load_hdf_iter_predict(tf.keras.utils.Sequence):
    """
    5.0: load the data per batch, include batch_size chunks(32*800*11)
    :return a iter batch_size data for training, val and test
    """
    def __init__(self, hdf_files, keys, batch_size=256):
        """
        initial the parameters
        :param hdf_files:
        :param keys:
        :param mode: train val test
        :param batch_size:
        """
        self.file = h5py.File(hdf_files, 'r')
        self.keys = keys
        self.batch_size = batch_size

    def __len__(self):
        """
        即数据第一维的数目/batch_size向上取整 (109/32)=4  至少得是1，为0会无法进去getitem函数
        且这里range的范围是__getitem__(self, idx) idx的循环范围
        :return:
        """
        return math.ceil(self.file['features'].shape[0] / self.batch_size) # 向上取整输入的是batch的数目

    def __getitem__(self, idx):
        """生成一个batch的数据"""

        while 1:
            f = self.file
            features = []
            # positions = []
            # read_ids = []

            features = f[self.keys][idx*self.batch_size:(idx+1)*self.batch_size, :]
            # positions = f[self.keys[1]][idx*self.batch_size:(idx + 1)*self.batch_size, :]
            # read_ids = f[self.keys[2]][idx*self.batch_size:(idx+1)*self.batch_size, :]

            return tf.convert_to_tensor(features),
            # tf.convert_to_tensor(positions), tf.convert_to_tensor(read_ids)  #








