"""
Version 2.0
make the SR as the only label in multi process
"""
import numpy as np
import os
import pandas as pd
import h5py
import argparse

from functools import partial
import multiprocessing as mp

label_symbols = ["*", "A", "C", "G", "T"]
feature_symbols = ["a", "c", "g", "t", "A", "C", "G", "T", "#", "*"]


def write_wrs(id,read,file):
    with open(file, mode='a') as f:
        seq_name = ">" + str(id) + "\n"
        f.write(seq_name)
        f.write(read)
        f.write("\n")
    f.close()



def get_file_list(file_path):
    fileList = []
    fileName = []
    for root, dirs, files in os.walk(file_path):

        for fileObj in files:
            #
            fileName.append(fileObj)
            path = os.path.join(root, fileObj).replace('\\', '/')
            fileList.append(path)
        return fileList


def find_insertions(base_pileup):
    """Finds all of the insertions in a given base's pileup string.

    Args:
        base_pileup: Single base's pileup string output from samtools mpileup
    Returns:
        insertions: list of all insertions in pileup string
        next_to_del: whether insertion is next to deletion symbol (should be ignored)
    """
    insertions = []
    idx = 0
    next_to_del = []
    while idx < len(base_pileup):
        if base_pileup[idx] == "+" and base_pileup[idx+1].isdigit():
            end_of_number = False
            insertion_bases_start_idx = idx+1
            while not end_of_number:
                if base_pileup[insertion_bases_start_idx].isdigit():
                    insertion_bases_start_idx += 1
                else:
                    end_of_number = True
            insertion_length = int(base_pileup[idx:insertion_bases_start_idx])
            inserted_bases = base_pileup[insertion_bases_start_idx:insertion_bases_start_idx+insertion_length]
            insertions.append(inserted_bases)
            next_to_del.append(True if base_pileup[idx - 1] in '*#' else False)
            idx = insertion_bases_start_idx + insertion_length + 1  # skip the consecutive base after insertion
        else:
            idx += 1
    return insertions, next_to_del


def calculate_positions(start_pos, end_pos, subreads, sr_coverage):
    """Calculates positions array from read pileup columns.

    Args:
        start_pos: Starting index of pileup columnsyjf
        end_pos: Ending index of pileup columns
        subreads: Array of subread strings from pileup file
        exclude_no_coverage_positions: Boolean specifying whether to include 0 coverage
        positions during position calculation
    Returns:
        positions: Array of tuples containing major/minor positions of pileup
    """
    positions = []
    # Calculate ref and insert positions

    for i in range(start_pos, end_pos):

        if sr_coverage[i] == 0:
            lr_pos = []
            lr_pos.append([i, 0])
            positions += lr_pos  #
            continue
        else:
            if type(subreads[i]) == float:  #
                subreads[i] = str('*')
            base_pileup = subreads[i].strip("^]").strip("$")

            # Get all insertions in pileup
            insertions, next_to_del = find_insertions(base_pileup)

            # Find length of maximum insertion
            longest_insertion = len(max(insertions, key=len)) if insertions else 0

            # Keep track of ref and insert positions in the pileup and the insertions
            # in the pileup.
            ref_insert_pos = []  # ref position for ref base pos in pileup, insert for additional inserted bases
            ref_insert_pos.append([i, 0])
            for j in range(longest_insertion):
                ref_insert_pos.append([i, j + 1])
            positions += ref_insert_pos

    return positions

def reencode_base_pileup(ref_base, pileup_str):
    """Re-encodes mpileup output of list of special characters to list of nucleotides.

    Args:
        ref_base : Reference nucleotide
        pileup_str : mpileup encoding of pileup bases relative to reference
    Returns:
        A pileup string of special characters replaced with nucleotides.
    """
    pileup = []
    for c in pileup_str:
        if c == "." or c == "*":  # wrs code or c == "*" 只要没有的一律拿长读本身代替
            pileup.append(ref_base)
        elif c == ",":
            pileup.append(ref_base.lower())
        else:
            pileup.append(c)
    return "".join(pileup)


def feature_encode(pileup):
    """
        feature generate for each seq in LR set with the information in SR reads which aligned to it.
        :param pileup:
        :return:
    """
    start_pos = 0
    end_pos = len(pileup)
    subreads = pileup[:, 4]
    sr_coverage = pileup[:, 3].astype("int")
    positions = calculate_positions(start_pos, end_pos, subreads, sr_coverage)
    positions = np.array(positions)

    num_features = len(feature_symbols)

    pileup_counts = np.zeros(shape=(len(positions), num_features + 1), dtype=np.float32)  #
    labels = np.zeros((len(positions),))  # gap, A, C, G, T (sparse format)
    for i in range(len(positions)):
        ref_position = positions[i][0]
        insert_position = positions[i][1]
        if type(subreads[ref_position]) == float:  #
            print("float error attention", pileup[ref_position,:])
            subreads[ref_position] = str('*')
        base_pileup = subreads[ref_position].strip("^]").strip("$")
        base_pileup = reencode_base_pileup(pileup[ref_position, 2], base_pileup)  #
        insertions, next_to_del = find_insertions(base_pileup)
        insertions_to_keep = []
        if sr_coverage[ref_position] != 0.:
            pileup_counts[i, 10] = 1

        # Remove all insertions which are next to delete positions in pileup
        for k in range(len(insertions)):
            if next_to_del[k] is False:
                insertions_to_keep.append(insertions[k])

        # Replace all occurrences of insertions from the pileup string
        for insertion in insertions:
            base_pileup = base_pileup.replace("+" + str(len(insertion)) + insertion, "")

        if insert_position == 0:  # No insertions for this position
            for j in range(len(feature_symbols)):
                pileup_counts[i, j] = base_pileup.count(feature_symbols[j])
            # Add draft base and base quality to encoding

        elif insert_position > 0:
            # Remove all insertions which are smaller than minor position being considered
            # so we only count inserted bases at positions longer than the minor position
            insertions_minor = [x for x in insertions_to_keep if len(x) >= insert_position]
            for j in range(len(insertions_minor)):
                inserted_base = insertions_minor[j][insert_position - 1]
                if inserted_base not in feature_symbols:
                    pileup_counts[i, feature_symbols.index("*")] += 1
                else:
                    pileup_counts[i, feature_symbols.index(inserted_base)] += 1
                # pileup_counts[i, feature_symbols.index(inserted_base)] += 1 # original code wrs

        label_base = feature_symbols[np.argmax(pileup_counts[i, 0:9])].upper()
        if label_base == "#":
            label_base = "*"
        labels[i] = label_symbols.index(label_base)
        # index = np.argmax(pileup_counts)

    return pileup_counts, labels, positions



def pack_feature_encode(mpileup_file):
    """
    split the pileup for per seq and pack it together again，
    then encode per seq_region separately and add the seq number info to each region,
    finally pack the nd.array together and return
    :return:
    """
    pileup_total = pd.read_csv(mpileup_file, delimiter="\t", header=None, quoting=3).values
    seq_name = pileup_total[:, 0]  # the first column of pileup is the seq_name set
    sequence = pileup_total[:, 2]
    truth_cov = pileup_total[:, 3]  # #
    sr_cov = pileup_total[:, 3]

    pileup_counts = np.empty(shape=(0, feature_size + 1), dtype=np.float32)
    positions = np.empty(shape=(0, 3), dtype=np.int64)
    labels = np.empty(shape=(0,), dtype=np.int64)

    cur_start = 0
    seq_num = 1  #

    for i in range(0, len(pileup_total) - 1):
        if seq_name[i] == seq_name[i + 1]:
            continue

        truth = truth_cov[cur_start: i + 1]
        sr = sr_cov[cur_start: i + 1]
        if not (np.any(sr)):  # 
            cur_end = i + 1
            read = ''.join([e for e in sequence[cur_start: cur_end]])
            write_wrs(seq_name[i], read, parsed_args.uncovered_fasta)  #
            read = ''
        else:
            cur_end = i + 1

            # test begin
            read = ''.join([e for e in sequence[cur_start: cur_end]])
            # print("covered: ", seq_name[i])
            # print(read)  covered.fasta no need write out
            read = ''
            # test end

            pileup = pileup_total[cur_start:cur_end, :]  # split pileup into per seq

            feature, label, feature_position = feature_encode(pileup)  # region for per seq

            if len(feature) == 0 or len(feature_position) == 0:  # if the seq has no coverage from SR or truth,ignore it
                cur_start = i + 1
                seq_num += 1  # seq_num解码的时候说不定还可以拿去还原序列号呢，所以先自加1
                continue

            feature_position = np.insert(feature_position, [0], seq_num, axis=1)

            assert len(feature) == len(label)

            pileup_counts = np.append(pileup_counts, feature, axis=0)  # 拼接
            labels = np.append(labels, label, axis=0)
            positions = np.append(positions, feature_position, axis=0)

        cur_start = cur_end
        seq_num += 1
        #print(seq_num)

    assert len(positions) == len(pileup_counts) or len(positions) == len(labels)

    print("position array are equal")

    return pileup_counts, positions, labels


def sliding_window(array, window, step=1, axis=0):
    """Generate chunks for encoding and labels.

    Args:
        array: Numpy array with the pileup counts
        window: Length of output chunks
        step: window minus chunk overlap
        axis: defaults to 0
    Returns:
        Iterator with chunks
    """
    chunk = [slice(None)] * array.ndim
    end = 0
    chunks = []
    for start in range(0, array.shape[axis] - window + 1, step):
        end = start + window
        chunk[axis] = slice(start, end)
        chunks.append(array[tuple(chunk)])
    if array.shape[axis] > end:
        start = array.shape[axis] - window
        chunk[axis] = slice(start, array.shape[axis])
        chunks.append(array[tuple(chunk)])
    return chunks


def encode(chunk_len, chunk_ovlp, input_file):  # chunk_len=1000, chunk_ovlp=200
    """
    the function inculde three parts: encode the feature and label,
    count the position according the alignment file
    slice the feature_encode and take account the overlaps=200
    next is feeding into the network
    :return:
    """

    encoding, encoding_positions, label = pack_feature_encode(input_file)  # generate the temp pileup data and position from the aligned data
    read_ids = encoding_positions[:, 0]
    max_test = 0  # test the max length before training
    length = 0
    for i in range(len(read_ids) - 1):
        if read_ids[i] == read_ids[i + 1]:
            length = length + 1
        else:
            if length > max_test:
                max_test = length
            length = 0
    print("max is:", max_test)  # test the max length before training
    assert (len(encoding) == len(label)), print("Encoding and label dimensions not as expected:", )
    assert (len(encoding_positions) == len(encoding)), print("Encoding and positions not as expected:", )

    encoding_chunks = sliding_window(encoding, chunk_len, step=chunk_len - chunk_ovlp)
    position_chunks = sliding_window(encoding_positions, chunk_len, step=chunk_len - chunk_ovlp)
    label_chunks = sliding_window(label, chunk_len, step=chunk_len - chunk_ovlp)
    read_id_chunks = sliding_window(read_ids, chunk_len, step=chunk_len - chunk_ovlp)

    return encoding_chunks, position_chunks, label_chunks, read_id_chunks


def generate_hdf5(chunk_len, chunk_overlap, output_file):
    """Generate encodings in multiprocess loop and save tensors to HDF5."""
    encode_func = partial(encode, chunk_len, chunk_overlap)

    # output data
    features = []  # features in column
    labels = []  # correct labeling
    positions = []  # track match/insert for stitching
    read_ids = []  # track folder name and windows

    label_idx = 0
    pool = mp.Pool(processes=60)
    for out in pool.imap(encode_func, file_list):
        if (label_idx + 1) % 100 == 0:
            print('Generated {} pileups'.format(label_idx + 1))
        (encoding_chunks, position_chunks, label_chunks, read_id_chunks) = out
        if encoding_chunks and position_chunks and label_chunks:
            if encoding_chunks[0].shape[0] == chunk_len \
                    and label_chunks[0].shape[0] == chunk_len \
                    and position_chunks[0].shape[0] == chunk_len:
                features += (encoding_chunks)
                labels += (label_chunks)
                positions += (position_chunks)
                read_ids += (read_id_chunks)
                label_idx += 1
    pool.close()
    print('Generated {} pileup files'.format(label_idx))
    features = np.stack(features, axis=0)
    labels = np.stack(labels, axis=0)
    positions = np.stack(positions, axis=0)
    h5_file = h5py.File(output_file, 'w')
    h5_file.create_dataset('features', data=features)
    h5_file.create_dataset('positions', data=positions)
    h5_file.create_dataset('labels', data=labels)
    h5_file.create_dataset('read_ids', data=read_ids)
    h5_file.close()

def check_file(file):
    if os.path.exists(file):
        print(file, " is exists!")
        os.remove(file)


def build_parser():
    """Setup option parsing for sample."""
    parser = argparse.ArgumentParser(description='Multi proprecess mpileup files and store encoded data in one HDF5 format.')
    parser.add_argument('--mpileup-folder', type=str, help='The path of split_mpileup foder')
    parser.add_argument('--output-file', type=str, help='Path to output HDF5 file.')
    parser.add_argument('--uncovered-fasta', type=str, help='The path of uncovered-fasta')

    args = parser.parse_args()

    return args

if __name__ == '__main__':
    feature_size = 10
    parsed_args = build_parser()
    check_file(parsed_args.uncovered_fasta)  # clear the output file to avoid overlying
    file_list = get_file_list(parsed_args.mpileup_folder)
    generate_hdf5(chunk_len=1000, chunk_overlap=200, output_file=parsed_args.output_file)




## 调试参数：--mpileup-folder ./testmpileup --output-file  ./testcovered_fasta.hdf --uncovered-fasta ./testuncovered.fasta