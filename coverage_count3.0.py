"""
count the coverage of long reads version 2.0
count the coverage states of pileup, and try to make it in multi-processes

"""


import os
import pandas as pd
import numpy as np
import h5py

from functools import partial
import multiprocessing as mp

def get_file_list(file_path):
    fileList = []
    fileName = []
    for root, dirs, files in os.walk(file_path):

        for fileObj in files:
            # 空列表写入遍历的文件名称，目录路径拼接文件名称
            fileName.append(fileObj)
            path = os.path.join(root, fileObj).replace('\\', '/')
            fileList.append(path)
        return fileList

def count_coverage(mpileup_file):
    pileup = pd.read_csv(mpileup_file, delimiter="\t", header=None, quoting=3).values # read it first
    sr_coverage = pileup[:,3].astype("int")

    count_0 = 0
    count_1 = 0

    for i in range(len(pileup)):
        if sr_coverage[i] != 0:
            count_1 += 1
        if sr_coverage[i] == 0:
            count_0 += 1


    return count_0, count_1, len(pileup)



def collect_coverage():

    count_coverages = partial(count_coverage) # 固定一部分参数

    # output data
    counts_0 = 0
    counts_1 = 0
    len_pileups = 0

    label_idx = 0
    pool = mp.Pool()
    for out in pool.imap(count_coverages, file_list):

        (count_0, count_1, len_pileup) = out
        counts_0 += count_0
        counts_1 += count_1

        len_pileups += len_pileup
        #print("Generated {} pileup files".format(label_idx))
        label_idx += 1
    pool.close()

    print("The rate of onlyS: ", (counts_1 / len_pileups))
    print("The rate of noS: ", (counts_0 / len_pileups))


if __name__ == '__main__':


    folders_to_encode = "./split_mpileup"  # the path that mpileup file is saved
    file_list = get_file_list(folders_to_encode)
    collect_coverage()
