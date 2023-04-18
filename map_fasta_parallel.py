import sqlite3
import time
from joblib import Parallel, delayed
import sys

def get_mapping_location_fasta(fasta_contig, min_kmer_consecutive, cursor):
    connection = sqlite3.connect("{}.db".format(cursor))
    cursor = connection.cursor()
    t = time.time()
    mapped_positions = {}
    cmd = """SELECT fasta_loc, cmap_loc, cmap_contig FROM overlap
    WHERE fasta_contig=\"{}\"
    ORDER BY fasta_loc;""".format(fasta_contig)

    mapped_positions_temp = set(cursor.execute(cmd).fetchall())
    for item in mapped_positions_temp:
        if item[0] not in mapped_positions:
             mapped_positions[item[0]] = []
        mapped_positions[item[0]].append("{},{}".format(item[1], item[2]))
    for fa in mapped_positions:
        mapped_positions[fa] = set(mapped_positions[fa])
    final_dict = {}
    final_list = []
    j = 0
    for fa_pos in mapped_positions:
        for cmap_pos in mapped_positions[fa_pos]:
            content = cmap_pos.split(',')
            temp = [fasta_contig, [], content[1], [], 0]

            temp[1].append(fa_pos)
            temp[3].append(int(content[0]))
            nx = "{},{}".format(temp[3][-1]+1, content[1])
            i = 0
            while temp[1][-1]+1 in mapped_positions and nx in mapped_positions[temp[1][-1]+1]:
                temp[1].append(temp[1][-1]+1)
                temp[3].append(temp[3][-1]+1)
                nx = "{},{}".format(temp[3][-1]+1, content[1])
                i += 1
            if temp[1][-1]+1 not in mapped_positions or  nx not in mapped_positions[temp[1][-1]+1]:
                if len(temp[1]) > min_kmer_consecutive:
                    temp[4] = len(temp[1])
                    final_dict[j] = temp
                    final_list.append(temp)
                    j += 1
    return [fasta_contig, final_list]


def get_kmer_n_fasta(fasta_contig, cursor):
    cmd = """SELECT pos_contig from fasta
        WHERE contig=\"{}\";""".format(fasta_contig)

    fasta_pos = cursor.execute(cmd).fetchall()
    return len(fasta_pos)


def get_mapping_pos_fasta(min_kmer_consecutive, fasta_list, cursor, db_name, n_threads):
    mapped_fasta = {}
    remove_fasta = []
    remove_cmap = []
    print("started getting fasta mapping positions")
    t = time.time()
    if min_kmer_consecutive == 'max':
        mapped_fasta_temp = Parallel(n_jobs=int(n_threads))(delayed(get_mapping_location_fasta)(fasta, get_kmer_n_fasta(fasta,cursor), db_name) for fasta in fasta_list)
    else:
        mapped_fasta_temp = Parallel(n_jobs=int(n_threads))(delayed(get_mapping_location_fasta)(fasta, int(min_kmer_consecutive), db_name) for fasta in fasta_list)

    for i in range(len(mapped_fasta_temp)):
        mapped_fasta[mapped_fasta_temp[i][0]] = mapped_fasta_temp[i][1]
    print("Parallel method took: {}".format(time.time()-t))
    for fasta in mapped_fasta:
        for item in mapped_fasta[fasta]:
            for i in range(len(item[1])):
                remove_fasta.append([fasta,item[1][i]])
                remove_cmap.append([item[2],item[3][i]])
    return mapped_fasta, remove_fasta, remove_cmap


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


def get_cmap_positions_mapped(mapped_fasta, cmap_list, cursor):
    print("started getting cmap positions mapped")
    t = time.time()
    mapped_cmap = {}
    for cmap in cmap_list:
        mapped_cmap[cmap] = {}
        positions = get_all_positions_in_cmap(cmap, cursor)
        for pos in positions:
            mapped_cmap[cmap][pos] = False

    for fasta in mapped_fasta:
        for overlaps in mapped_fasta[fasta]:
            #Overlaps: ['fasta_name', [list of fasta positions], 'cmap name', [list of cmap positions], 'overlap len']
            cmap = overlaps[2]
            for i in range(overlaps[4]):
                if mapped_cmap[cmap][overlaps[3][i]] == False:
                    mapped_cmap[cmap][overlaps[3][i]] = (fasta, overlaps[1][i], overlaps[4])
                else:
                    if mapped_cmap[cmap][overlaps[3][i]][2] < overlaps[4]:
                        mapped_cmap[cmap][overlaps[3][i]] = (fasta, overlaps[1][i], overlaps[4])
    print('getting cmap positions mapped took: {}'.format(time.time()-t))
    return mapped_cmap


def get_fasta_loc(fasta, pos_fasta, cursor):
    cmd = """
        SELECT * FROM fasta_contigs
        WHERE contig=\"{}\"
        AND pos_contig={}""". format(fasta, pos_fasta)

    tmp_pos = cursor.execute(cmd).fetchall()[0]
    fasta_positions = []
    for i in range(3,len(tmp_pos)-1):
        fasta_positions.append("{},{}".format(tmp_pos[i],tmp_pos[i+1]))
    return fasta_positions

def get_cmap_loc(cmap, pos_cmap, cursor):
    cmd = """
        SELECT * FROM cmap_contigs
        WHERE contig=\"{}\"
        AND pos_contig={}""". format(cmap, pos_cmap)

    tmp_pos = cursor.execute(cmd).fetchall()[0]
    cmap_positions = []
    for i in range(3,len(tmp_pos)-1):
        cmap_positions.append("{},{}".format(tmp_pos[i],tmp_pos[i+1]))
    return cmap_positions

def change_cmap_mapped_format(mapped_cmap, cursor):
    print("started changing format")
    t = time.time()
    final_dict = {}
    sorted_pos = {}
    for cmap in mapped_cmap:
        final_dict[cmap] = {}
        sorted_pos[cmap] = []

    for cmap in mapped_cmap:
        cmap_positions = get_all_positions_in_cmap(cmap, cursor)
        for position in cmap_positions:
            cmap_locations = get_cmap_loc(cmap, position, cursor)
            for loc in cmap_locations:
                if loc not in sorted_pos[cmap]:
                    sorted_pos[cmap].append(loc)
            if position not in mapped_cmap[cmap] or mapped_cmap[cmap][position] == False:
                for loc in cmap_locations:
                    if loc not in final_dict[cmap]:
                        final_dict[cmap][loc] = False
            else:
                fasta_locations = get_fasta_loc(mapped_cmap[cmap][position][0], mapped_cmap[cmap][position][1], cursor)
                for i in range(len(cmap_locations)):
                    if cmap_locations[i] not in final_dict[cmap]:
                        final_dict[cmap][cmap_locations[i]] = (fasta_locations[i], mapped_cmap[cmap][position][0])
                    else:
                        if final_dict[cmap][cmap_locations[i]] == False:
                            final_dict[cmap][cmap_locations[i]] = (fasta_locations[i], mapped_cmap[cmap][position][0])
    print("Changing mapping format took: {}".format(time.time()-t))
    return final_dict, sorted_pos


def merge_fa_cmap(mapped_cmaps, cmap, enzyme_sequence, fasta_sequence_dict, sorted_pos, sequence):
    if sequence == '':
        for position in sorted_pos[cmap]:
            if mapped_cmaps[cmap][position] == False:
                content= position.split(",")
                n_count = int(content[1]) - int(content[0])
                sequence += enzyme_sequence
                for i in range(n_count):
                    sequence += 'N'
            else:
                content = mapped_cmaps[cmap][position][0].split(",")
                seq = fasta_sequence_dict[mapped_cmaps[cmap][position][1]][int(content[0]):int(content[1])]
                sequence += seq
    else:
        for position in sorted_pos[cmap]:
            if mapped_cmaps[cmap][position] != False:
                content = mapped_cmaps[cmap][position][0].split(",")
                seq = fasta_sequence_dict[mapped_cmaps[cmap][position][1]][int(content[0]):int(content[1])]
                pos = position.split(',')
                sequence = sequence[0:int(pos[0])]+seq+sequence[int(pos[1])::]
    return [cmap, sequence]

def merge(db_name, enzyme_sequence, fasta_sequence_dict, min_kmer_consecutive, n_threads, all_seq = {}):
    print("started merging")
    t = time.time()
    connection = sqlite3.connect("{}.db".format(db_name))
    cursor = connection.cursor()
    fasta_list = list(fasta_sequence_dict.keys())
    cmap_list = get_cmap_list(db_name)

    fasta_mapping, remove_fasta, remove_cmap = get_mapping_pos_fasta(min_kmer_consecutive, fasta_list, cursor, db_name, n_threads)
    mapped_cmap = get_cmap_positions_mapped(fasta_mapping, cmap_list, cursor)
    cmap_mapping, sorted_cmap_pos = change_cmap_mapped_format(mapped_cmap, cursor)
    for cmap in cmap_list:
        if cmap not in all_seq:
            all_seq[cmap] = ''

    sequences_temp = Parallel(n_jobs=int(n_threads))(delayed(merge_fa_cmap)(cmap_mapping, cmap, enzyme_sequence, fasta_sequence_dict, sorted_cmap_pos, all_seq[cmap]) for cmap in cmap_list)

    for item in sequences_temp:
        all_seq[item[0]] = item[1]

    print("everything together took: {}".format(time.time()-t))
    return all_seq, remove_fasta, remove_cmap

