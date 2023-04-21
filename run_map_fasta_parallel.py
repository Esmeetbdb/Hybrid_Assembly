import sqlite3
import time
from joblib import Parallel, delayed
import sys
import numpy as np
import functions_map_fasta_parallel as r


def parallel_init_cmap_dict(cmap_list, np_positions_cmap, np_name_id_cmap, cmap_pos_dict, n_threads):
    """
    Call the init_cmap_dict function in parallel
    :param cmap_list: List of all cmap names
    :param np_positions_cmap: numpy array with format [[pos_1 pos_2 pos_3 ... pos_n]]. The array has 1 row per k-mer
        per cmap
    :param np_name_id_cmap: numpy array with format [[cmap_name position_id]] with a row per cmap and position_id
    :param cmap_pos_dict: dictionary with format {cmap_name:[list of all position_ids in this cmap]}
    :param n_threads: the number of threads that can be used for computation
    :return: a list with format [[cmap_name, {cmap_location:False}, [sorted_locations]]]
    """
    cmap_temp = Parallel(n_jobs=int(n_threads))(
        delayed(r.init_cmap_dict)(cmap_name, np_positions_cmap, np_name_id_cmap, cmap_pos_dict[cmap_name]) for cmap_name
        in cmap_list)

    return cmap_temp


def get_mapping_pos_fasta(min_kmer_consecutive, fasta_list, db_name, n_threads, np_positions_fasta, np_name_id_fasta,
                          np_positions_cmap, np_name_id_cmap):
    """
    Run get_mapping_location_fasta in parallel
    :param min_kmer_consecutive: The minimum number of consecutive overlapping k-mers for an overlap to be considered
    :param fasta_list: list containing the names of all fasta contigs
    :param db_name: name of the database
    :param n_threads: the number of threads that can be used for computation
    :param np_positions_fasta: numpy array with format [[pos_1 pos_2 pos_3 ... pos_n]]. The array has 1 row per k-mer
        per fasta
    :param np_name_id_fasta: numpy array with format [[fasta_name position_id]] with a row per fasta and position_id
    :param np_positions_cmap: numpy array with format [[pos_1 pos_2 pos_3 ... pos_n]]. The array has 1 row per k-mer
        per cmap
    :param np_name_id_cmap: numpy array with format [[cmap_name position_id]] with a row per cmap and position_id
    :return: list with format [[{cmap_name:{"cmap_start:cmap_stop":["fasta_start:fasta_stop", fasta_name]}},
    remove_fasta, remove_cmap]]
    """
    mapped_fasta_temp = Parallel(n_jobs=int(n_threads))(
        delayed(r.get_mapping_location_fasta)(fasta, int(min_kmer_consecutive), db_name, np_positions_cmap,
                                              np_name_id_cmap, np_positions_fasta,
                                              np_name_id_fasta) for fasta in fasta_list)
    return mapped_fasta_temp


def get_all_mapping_info(db_name, cmap_list, min_kmer_consecutive, fasta_list, n_threads):
    """
    Function to obtain all mapping information, of which sequence regions overlap with which cmap regions
    :param db_name: name of the database
    :param cmap_list: list of all the cmaps used
    :param min_kmer_consecutive: The minimum number of consecutive overlapping k-mers for an overlap to be considered
    :param fasta_list: list of all the fasta contigs
    :param n_threads:
    :return: the number of threads that can be used for computation
    """
    connection = sqlite3.connect("{}.db".format(db_name))
    cursor = connection.cursor()

    np_positions_fasta, np_name_id_fasta = r.get_positions('fasta', cursor)
    np_positions_cmap, np_name_id_cmap = r.get_positions('cmap', cursor)

    cmap_pos_dict = r.get_cmap_pos_dict(cmap_list, cursor)

    init_cmap_dict = parallel_init_cmap_dict(cmap_list, np_positions_cmap, np_name_id_cmap, cmap_pos_dict, n_threads)
    del cmap_pos_dict
    empty_cmap_dict, sorted_pos_cmap_lict = r.get_empty_cmap_dict(init_cmap_dict)
    del init_cmap_dict

    mapped_fasta_temp = get_mapping_pos_fasta(min_kmer_consecutive, fasta_list, db_name, n_threads,
                                              np_positions_fasta, np_name_id_fasta, np_positions_cmap, np_name_id_cmap)
    del np_positions_fasta
    del np_name_id_fasta
    del np_positions_cmap
    del np_name_id_cmap

    filled_cmap_dict, remove_fasta, remove_cmap = r.enter_info_cmap_dict(mapped_fasta_temp, empty_cmap_dict)

    return filled_cmap_dict, sorted_pos_cmap_lict, remove_fasta, remove_cmap


def merge(db_name, enzyme_sequence, fasta_sequence_dict, min_kmer_consecutive, n_threads, all_seq={}):
    print("started merging")
    t = time.time()
    mapped_pos_number = 0

    fasta_list = list(fasta_sequence_dict.keys())
    cmap_list = r.get_cmap_list(db_name)

    cmap_mapping, sorted_pos, remove_fasta, remove_cmap = get_all_mapping_info(db_name, cmap_list,
                                                                                 min_kmer_consecutive,
                                                                                 fasta_list, n_threads)

    for cmap in cmap_list:
        if cmap not in all_seq:
            all_seq[cmap] = ''

    sequences_temp = Parallel(n_jobs=int(n_threads))(
        delayed(r.merge_fa_cmap)(cmap_mapping, cmap, sorted_pos, enzyme_sequence, fasta_sequence_dict, all_seq[cmap])
        for cmap in cmap_list)

    for item in sequences_temp:
        all_seq[item[0]] = item[1]
        mapped_pos_number += item[2]
    print("everything together took: {}".format(time.time() - t))
    return all_seq, remove_fasta, remove_cmap, mapped_pos_number
