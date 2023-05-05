import time
import sqlite3
import numpy as np
cimport numpy as np
from joblib import Parallel, delayed
import os
import psutil
cimport cython
@cython.boundscheck(False)
@cython.wraparound(False)

def get_positions(table_type, cursor):
    """
    Extract position ids and enzyme sites from sql table and store into memory as np
    arrays
    :param table_type: fasta or cmap
    :param cursor: connection to database
    :return: np_positions: 2D numpy array with format [[pos_1 pos_2 ... pos_n]]
    :return: np_name_id: 2D numpy array with format [[contig/cmap_name position_id]]
    """
    cmd = """
        SELECT * FROM {}_contigs;
        """.format(table_type)

    tmp_pos = cursor.execute(cmd).fetchall()

    contig_array = np.array([v[1] for v in tmp_pos])
    pos_id_array = np.array([v[2] for v in tmp_pos], dtype=np.int32)

    pos_dict = {}
    pos_np_list = []
    for i in range(3, len(tmp_pos[1])):
        pos_dict["position_{}".format(i - 2)] = np.array([v[i] for v in tmp_pos], dtype=np.int32)
        pos_np_list.append(pos_dict["position_{}".format(i - 2)])
    pos_np = tuple(pos_np_list)
    np_positions = np.column_stack(pos_np)
    np_name_id = np.column_stack((contig_array,pos_id_array))
    del pos_id_array
    del contig_array
    del pos_np
    del pos_dict
    del tmp_pos
    del pos_np_list

    return np_positions, np_name_id


def get_cmap_pos_dict(cmap_list, cursor):
    """
    Function to obtain a dictionary that contains a key for each cmap, with as value a list of all the position_ids in
    that cmap
    :param cmap_list: list of all cmap names
    :param cursor: Connection to the database
    :return: {cmap_name:[position_ids]}
    """
    cmap_pos_dict = {}
    for cmap in cmap_list:
        cmap_pos_dict[cmap] = get_all_positions_in_cmap(cmap, cursor)

    return cmap_pos_dict


def get_empty_cmap_dict(cmap_temp):
    """
    Changes the format of cmap_temp to two dictionaries. One is a dictionary with a cmap as each key and the value is
    a dictionary with format {"cmap_start:cmap_end":False}. The other dictionary also has a cmap as each key, the value
    is a list that contains the sorted intervals between enzyme sites
    :param cmap_temp: a list with format [[cmap_name, {cmap_location:False}, [sorted_locations]]]
    :return: {cmap:{"cmap_start:cmap_end":False}}, {cmap:[sorted_intervals]}
    """
    cmap_dict = {}
    sorted_cmap_pos = {}
    for i in range(len(cmap_temp)):
        cmap_dict[cmap_temp[i][0]] = cmap_temp[i][1]
        sorted_cmap_pos[cmap_temp[i][0]] = cmap_temp[i][2]

    return cmap_dict, sorted_cmap_pos


def get_overlapping_fasta(fasta_contig, db_name):
    """
    Function to extract the overlaps of 1 fasta location to all cmaps from the overlap sql table
    :param fasta_contig: name of the fasta contig of which overlapping sequences will be extracted
    :param db_name: name of the sql database
    :return: set containing the overlap information in the format((fasta_loc, cmap_loc, cmap_contig))
    """
    connection = sqlite3.connect("{}.db".format(db_name))
    cursor = connection.cursor()
    cmd = """SELECT fasta_loc, cmap_loc, cmap_contig FROM overlap
        WHERE fasta_contig=\"{}\"
        ORDER BY fasta_loc;""".format(fasta_contig)

    mapped_positions_temp = set(cursor.execute(cmd).fetchall())
    return mapped_positions_temp


def change_format_overlap_fasta(mapped_positions_temp):
    """
    This funcion changes the output of get_overlapping_fasta to a dictionary.
    :param mapped_positions_temp: output of get_overlapping_fasta
    :return: dictionary containing the overlap information of each fasta in the format:
     {fasta_loc:("<cmap_loc>,<cmap_contig>")}
    """
    mapped_positions = {}

    for mapping_pos in mapped_positions_temp:
        mapped_positions[mapping_pos[0]] = []

    for mapping_pos in mapped_positions_temp:
        mapped_positions[mapping_pos[0]].append("{},{}".format(mapping_pos[1], mapping_pos[2]))

    for fasta_loc in mapped_positions:
        mapped_positions[fasta_loc] = set(mapped_positions[fasta_loc])

    return mapped_positions


def sort_overlap_output_lists(fasta_locations_list, cmap_locations_list, names_length_list):
    """
    Function that puts all different output lists in one big lists and sorts from longest to shortest overlap
    while making sure all lists are sorted in the same way so data that belongs together has the same index
    :param fasta_locations_list:
    :param cmap_locations_list:
    :param names_length_list:
    :return: sorted list (sorted by overlap len reverse) with format [[names_length_list],[fasta_locations_list],[cmap_locations_list]]
    """
    sorted_data = []

    for i in range(len(fasta_locations_list)):
        sorted_data.append([names_length_list[i], fasta_locations_list[i], cmap_locations_list[i]])

    sorted_data.sort(reverse=True, key=lambda sorted_data: sorted_data[0][2])

    return sorted_data


def get_consecutive_overlapping_kmers(mapped_positions, min_kmer_consecutive, fasta_contig):
    """
    Function that extracts consecutive overlapping k-mers from the overlap dictionary created in
    change_format_overlap_fasta
    :param mapped_positions: Output of change_format_overlap_fasta
    :param min_kmer_consecutive: The minimum number of consecutive overlapping k-mers for an overlap to be considered
    :param fasta_contig: The fasta contig of which the overlaps are examined
    :return: final_fasta_loc: list that contains a list for each consecutive overlap, the list contains each position_id
    of the fasta that overlaps with a cmap
    :return: final_cmap_loc: list that contains a list for each consecutive overlap, the list contains each position_id
    of the cmap that overlaps with the fasta
    :return: final_names_overlap_len: list that contains a list for each consecutive overlap,
    the list contains the fasta contig name, the cmap name and the length of the overlap between these.
    """
    final_fasta_loc = []
    final_cmap_loc = []
    final_names_overlap_len = []
    skip_dict = {}

    for fa_pos in mapped_positions:
        for overlap_cmaps in mapped_positions[fa_pos]:
            content = overlap_cmaps.split(',')
            cmap_loc = int(content[0])
            cmap_name = content[1]
            k = "{},{},{}".format(fa_pos, cmap_loc, cmap_name)
            if skip_dict.get(k) != None:
                continue
            skip_dict["{},{},{}".format(fa_pos, cmap_loc, cmap_name)] = 1
            temp_fasta_loc = [int(fa_pos)]
            temp_cmap_loc = [cmap_loc]
            temp_names_overlap_len = [fasta_contig, cmap_name, 0]

            temp_fasta_loc.append(fa_pos)
            temp_cmap_loc.append(cmap_loc)
            nx = "{},{}".format(temp_cmap_loc[-1] + 1, cmap_name)
            while temp_fasta_loc[-1] + 1 in mapped_positions and nx in mapped_positions[temp_fasta_loc[-1] + 1]:
                temp_fasta_loc.append(temp_fasta_loc[-1] + 1)
                temp_cmap_loc.append(temp_cmap_loc[-1] + 1)
                nx = "{},{}".format(temp_cmap_loc[-1] + 1, cmap_name)
                skip_dict["{},{},{}".format(temp_fasta_loc[-1] + 1, temp_cmap_loc[-1] + 1, cmap_name)] = 1
                temp_names_overlap_len[2] += 1
            if temp_fasta_loc[-1] + 1 not in mapped_positions or nx not in mapped_positions[temp_fasta_loc[-1] + 1]:
                if temp_names_overlap_len[2] >= min_kmer_consecutive:
                    final_fasta_loc.append(temp_fasta_loc)
                    final_cmap_loc.append(temp_cmap_loc)
                    final_names_overlap_len.append(temp_names_overlap_len)
    return final_fasta_loc, final_cmap_loc, final_names_overlap_len


def get_idx_overlapping_fasta_locations(sorted_data):
    """
    Function that obtains the indexes of overlaps in the sorted data list that have overlapping fasta locations
    :param sorted_data:
    :return:
    """
    delete_list = set([])

    for i in range(len(sorted_data)):
        if i in delete_list:
            continue
        for j in range(len(sorted_data)):
            if (j in delete_list) or (i >= j):
                continue
            elif any(np.isin(sorted_data[i][1], sorted_data[j][1])):
                delete_list.add(j)
    return delete_list



def get_pos_from_loc(contig_name, loc, np_positions, np_name_id):
    """
    Takes the position id and contig name and returns the actual location on the cmap/sequence
    :param contig_name: name of the contig
    :param loc: position_id
    :param np_positions: numpy array with format [[pos_1 pos_2 pos_3 ... pos_n]]. The array has 1 row per k-mer
    :param np_name_id: numpy array with format [[name position_id]] with a row per cmap/fasta and position_id
    :return: list of positions with format ["start:end", "start:end"] for each pair of enzyme sites with the position id
    """
    conditions = (np_name_id[:,0] == contig_name) & (np_name_id[:,1] == str(loc))
    row_idx = np.where(conditions)[0]
    positions = []
    for i in range(np_positions.shape[1]-1):
        positions.append("{},{}".format(np_positions[row_idx,i][0], np_positions[row_idx, i + 1][0]))
    return positions


def change_output_format(sorted_data, np_positions_cmap, np_name_id_cmap, np_positions_fasta, np_name_id_fasta):
    """
    Changes the format of the consecutive overlap from 3 lists to a dictionary
    :param sorted_data: list of lists that contains the filtered consecutive overlap information
    :param np_positions_cmap: numpy array with format [[pos_1 pos_2 pos_3 ... pos_n]]. The array has 1 row per k-mer
        per cmap
    :param np_name_id_cmap: numpy array with format [[cmap_name position_id]] with a row per cmap and position_id
    :param np_positions_fasta: numpy array with format [[pos_1 pos_2 pos_3 ... pos_n]]. The array has 1 row per k-mer
        per fasta
    :param np_name_id_fasta: numpy array with format [[fasta_name position_id]] with a row per fasta and position_id
    :return: final_output: dictionary with format {cmap_name:{"cmap_star:cmap_end":["fasta_start:fasta_end", contig]}}
    :return: remove_fasta: List of fasta position ids that have been mapped and should be removed before further
    alignment rounds
    :return: remove_cmap: List of cmap position ids that have been mapped and should be removed before further
    alignment rounds
    """
    final_output = {}
    remove_fasta = []
    remove_cmap = []
    for i in range(len(sorted_data)):
        cmap_name = sorted_data[i][0][1]
        fasta_name = sorted_data[i][0][0]
        alignment_length = sorted_data[i][0][2]
        cmap_pos = []
        fasta_pos = []
        t = time.time()
        for j in range(len(sorted_data[i][1])):
            remove_fasta.append([fasta_name,sorted_data[i][1][j]])
            remove_cmap.append([cmap_name,sorted_data[i][2][j]])
            if j == 0:
                cmap_pos += get_pos_from_loc(cmap_name, sorted_data[i][2][j], np_positions_cmap, np_name_id_cmap)
                fasta_pos += get_pos_from_loc(fasta_name, sorted_data[i][1][j], np_positions_fasta, np_name_id_fasta)

            else:
                cmap_pos.append(get_pos_from_loc(cmap_name, sorted_data[i][2][j], np_positions_cmap, np_name_id_cmap)[-1])
                fasta_pos.append(get_pos_from_loc(fasta_name, sorted_data[i][1][j], np_positions_fasta, np_name_id_fasta)[-1])
        t = time.time()
        if final_output.get(cmap_name) == None:
            final_output[cmap_name] = {}
        for j in range(len(cmap_pos)):
            if final_output[cmap_name].get(cmap_pos[j]) == None:
                final_output[cmap_name][cmap_pos[j]] = [fasta_pos[j], fasta_name, alignment_length]
            else:
                if alignment_length > final_output[cmap_name][cmap_pos[j]][2]:
                    final_output[cmap_name][cmap_pos[j]] = [fasta_pos[j], fasta_name, alignment_length]
    return final_output, remove_fasta, remove_cmap


def get_mapping_location_fasta(fasta_contig, min_kmer_consecutive, db_name, np_positions_cmap, np_name_id_cmap,
                               np_positions_fasta, np_name_id_fasta):
    """
    Function to obtain all the locations where a fasta contig overlaps with any cmap. The function then filters these
    and returns the filtered output.
    :param fasta_contig: Name of the fasta contig of interest
    :param min_kmer_consecutive: The minimum number of consecutive overlapping k-mers for an overlap to be considered
    :param db_name: name of the database
    :param np_positions_cmap: numpy array with format [[pos_1 pos_2 pos_3 ... pos_n]]. The array has 1 row per k-mer
        per cmap
    :param np_name_id_cmap: numpy array with format [[cmap_name position_id]] with a row per cmap and position_id
    :param np_positions_fasta: numpy array with format [[pos_1 pos_2 pos_3 ... pos_n]]. The array has 1 row per k-mer
        per fasta
    :param np_name_id_fasta: numpy array with format [[fasta_name position_id]] with a row per fasta and position_id
    :return: list with format [{cmap_name:{"cmap_start:cmap_stop":["fasta_start:fasta_stop", fasta_name]}},
    remove_fasta, remove_cmap]
    """
    t = time.time()
    mapped_pos_temp = get_overlapping_fasta(fasta_contig, db_name)
    mapped_positions = change_format_overlap_fasta(mapped_pos_temp)
    del mapped_pos_temp
    
    final_fasta_loc, final_cmap_loc, final_names_overlap_len = get_consecutive_overlapping_kmers(mapped_positions,
                                                                                                 min_kmer_consecutive,
                                                                                                 fasta_contig)

    sorted_data = sort_overlap_output_lists(final_fasta_loc, final_cmap_loc, final_names_overlap_len)

    del final_fasta_loc
    del final_cmap_loc
    del final_names_overlap_len

    overlapping_idx = get_idx_overlapping_fasta_locations(sorted_data)

    for i in sorted(overlapping_idx, reverse=True):
        del sorted_data[i]

    final_output, remove_fasta, remove_cmap = change_output_format(sorted_data,
                                                                   np_positions_cmap,
                                                                   np_name_id_cmap,
                                                                   np_positions_fasta,
                                                                   np_name_id_fasta)
    return [final_output, remove_fasta, remove_cmap]


def init_cmap_dict(cmap_name, np_positions_cmap, np_name_id_cmap, positions):
    """
    Function that initializes a dictionary, each interval between two enzyme sites is a key. The format is:
    {"cmap_start:cmap_end":False}. The function also creates a list with the ordered positions.
    :param cmap_name: Name of the cmap
    :param np_positions_cmap: numpy array with format [[pos_1 pos_2 pos_3 ... pos_n]]. The array has 1 row per k-mer
        per cmap
    :param np_name_id_cmap: numpy array with format [[cmap_name position_id]] with a row per cmap and position_id
    :param positions: list of all position_ids in the cmap
    :return: a list with format [cmap_name, {cmap_location:False}, [sorted_locations]]
    """
    mapped_cmap = {}
    sorted_pos = []
    mapped_cmap[cmap_name] = {}
    for i in range(len(positions)):
        if i == 0:
            sorted_pos += get_pos_from_loc(cmap_name, positions[i], np_positions_cmap, np_name_id_cmap)
            cmap_pos = get_pos_from_loc(cmap_name, positions[i], np_positions_cmap, np_name_id_cmap)
            for j in range(len(cmap_pos)):
                mapped_cmap[cmap_pos[j]] = False
        else:
            cmap_pos = get_pos_from_loc(cmap_name, positions[i], np_positions_cmap, np_name_id_cmap)[-1]
            sorted_pos.append(cmap_pos)
            mapped_cmap[cmap_pos] = False
    return [cmap_name, mapped_cmap, sorted_pos]


def get_all_positions_in_cmap(cmap, cursor):
    cmap_positions = []
    cmd = """SELECT pos_contig FROM cmap
            WHERE contig=\"{}\"
            ORDER BY pos_contig ASC;""".format(cmap)

    cmap_positions_temp = cursor.execute(cmd).fetchall()
    for pos in cmap_positions_temp:
        cmap_positions.append(pos[0])

    return cmap_positions


def get_cmap_list(db_name):
    connection = sqlite3.connect("{}.db".format(db_name))
    print("{}.db".format(db_name))
    cursor = connection.cursor()
    cmap_list = []
    get_cmaps = """
        SELECT contig FROM cmap_contigs;"""
    cmap_list_cmd = cursor.execute(get_cmaps)
    for item in cmap_list_cmd.fetchall():
        if item[0] not in cmap_list:
            cmap_list.append(item[0])
    return cmap_list


def enter_info_cmap_dict(mapped_fasta_temp, cmap_dict):
    """
    Fills up the empty initialized cmap dictionary using the fasta mapping information
    :param mapped_fasta_temp: list with format [[{cmap_name:{"cmap_start:cmap_stop":["fasta_start:fasta_stop", fasta_name]}},
    remove_fasta, remove_cmap]]
    :param cmap_dict: {cmap:{"cmap_start:cmap_end":False}}
    :return: {cmap:{"cmap_start:cmap_end":["Fasta_strat:Fasta_end", fasta_contig]}},
            [["fasta_contig",fasta_pos_id]], [["cmap_contig",cmap_pos_id]]
    """
    remove_fasta = []
    remove_cmap = []
    for i in range(len(mapped_fasta_temp)):
        remove_fasta += mapped_fasta_temp[i][1]
        remove_cmap += mapped_fasta_temp[i][2]
        tmp_dict = mapped_fasta_temp[i][0]

        for cmap_name in tmp_dict:
            for pos in tmp_dict[cmap_name]:
                if not cmap_dict[cmap_name].get(pos):
                    cmap_dict[cmap_name][pos] = tmp_dict[cmap_name][pos]
                else:
                    if tmp_dict[cmap_name][pos][2] > cmap_dict[cmap_name][pos][2]:
                        cmap_dict[cmap_name][pos] = tmp_dict[cmap_name][pos]
    return cmap_dict, remove_fasta, remove_cmap


def merge_fa_cmap(mapped_cmaps, cmap, sorted_pos, enzyme_sequence, fasta_sequence_dict, sequence):
    mapped_pos_number = 0
    if sequence == '':
        for position in sorted_pos:
            if not mapped_cmaps[cmap][position]:
                content = position.split(",")
                n_count = int(content[1]) - int(content[0])-len(enzyme_sequence)
                sequence += enzyme_sequence
                for i in range(n_count):
                    sequence += 'N'
            else:
                content = mapped_cmaps[cmap][position][0].split(",")
                seq = fasta_sequence_dict[mapped_cmaps[cmap][position][1]][int(content[0]):int(content[1])]
                sequence += seq
                mapped_pos_number += 1
    else:
        #change so that the seq is a list that is changed
        for position in sorted_pos[cmap]:
            if mapped_cmaps[cmap][position] != False:
                content = mapped_cmaps[cmap][position][0].split(",")
                seq = fasta_sequence_dict[mapped_cmaps[cmap][position][1]][int(content[0]):int(content[1])]
                pos = position.split(',')
                sequence = sequence[0:int(pos[0])] + seq + sequence[int(pos[1])::]
                mapped_pos_number += 1
    return [cmap, sequence, mapped_pos_number]
