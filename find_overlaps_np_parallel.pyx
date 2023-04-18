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

def compare_1(fasta_list, cmap_array, deviation, id, contig, position_id):
    t = time.time()
    ranges = []
    output_list = []
    conditions = []
    for j in range(3, len(fasta_list)):
        min_val = round((fasta_list[j] - deviation), 0)
        max_val = round((fasta_list[j] + deviation), 0)
        conditions.append(((cmap_array[:, j-3] >= min_val) & (cmap_array[:, j-3] <= max_val)))

    conditions_array = np.column_stack(conditions)
    row_idx = np.where(np.all(conditions_array, axis = 1))[0]

    for row in row_idx:
        output_list.append([None, fasta_list[1], fasta_list[0], fasta_list[2], contig[row], id[row], position_id[row]])


#    print('One overlap: {}'.format(time.time()-t))
    return output_list


def find_overlaps(db_name, deviation, n_threads):
    print(db_name)
    connection = sqlite3.connect("{}.db".format(db_name))
    cursor = connection.cursor()
    cursor.execute("""DROP TABLE IF EXISTS overlap""")
    sqlite3.register_adapter(np.int64, int)
    sqlite3.register_adapter(np.int32, int)

    create_overlap_table = """
    CREATE TABLE overlap (
    id INT PRIMARY KEY,
    fasta_contig VARCHAR(255) NOT NULL,
    fasta_id INT NOT NULL,
    fasta_loc INT NOT NULL,
    cmap_contig VARCHAR(255) NOT NULL,
    cmap_id INT NOT NULL,
    cmap_loc INT NOT NULL
    );
    """
    t = time.time()
    cursor.execute(create_overlap_table)

    GET_CMAP = """
    SELECT * FROM cmap;
    """

    cmaps_temp = cursor.execute(GET_CMAP).fetchall()
    id_list = np.array([v[0] for v in cmaps_temp], dtype=np.int64)
    contig = np.array([v[1] for v in cmaps_temp])
    position_id = np.array([v[2] for v in cmaps_temp], dtype=np.int32)

    distances_dict = {}
    cmap_np_list = []
    for i in range(3, len(cmaps_temp[1])):
        distances_dict["distance_{}".format(i-2)] = np.array([v[i] for v in cmaps_temp], dtype=np.int32)
        cmap_np_list.append(distances_dict["distance_{}".format(i-2)])
    cmap_np = tuple(cmap_np_list)
    del(distances_dict)
    del(cmap_np_list)
    del(cmaps_temp)
    
    cmap = np.column_stack(cmap_np)
    del(cmap_np)
    print("hello {}".format(time.time() -t))
    FIND_OVERLAP = """
    SELECT * FROM fasta;
    """
    print("selecting all pos from FA tooks: {}".format(time.time()-t))
    t = time.time()
    id = 1
    x = cursor.execute(FIND_OVERLAP)
    y = x.fetchall()
    del(x)
    print(len(y))
    pop_overlap_temp = Parallel(n_jobs=int(n_threads))(delayed(compare_1)(fasta_list, cmap, deviation, id_list, contig, position_id) for fasta_list in y)
    print("Finding overlaps took: {}".format(time.time()-t))
    t2 = time.time()
    id = 1
    pop_overlap = []
    for item in pop_overlap_temp:
        for item2 in item:
            item2[0] = id
            id += 1
            pop_overlap.append(item2)
    pid = os.getpid()
    python_process = psutil.Process(pid)
    memoryUse = python_process.memory_info()[0]/2.**30  # memory use in GB...I think
    print(psutil.virtual_memory())
    print('memory use: {}'.format(memoryUse))
    cursor.executemany('INSERT INTO overlap VALUES (?,?,?,?,?,?,?)', pop_overlap)
    connection.commit()
    print("Inserting overlaps took: {}".format(time.time()-t2))	
    print("finding all overlaps took: {}.\nFound {} overlaps".format(time.time()-t,len(pop_overlap)))

