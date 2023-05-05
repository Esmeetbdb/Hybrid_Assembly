import sqlite3
import time

def get_positions(cmap):
    """
    Function that extracts all the positions that had enzyme binding as indicated by the cmap file
    :param cmap: String that contains the path to the cmap file to be used during alignment
    :return: positions: Dictionary with format: {'contig/cmap name' : [ordered list of enzyme site locations]}
    :return: contig_length: Dictionary with format: {'contig/cmap name' : length of the cmap}
    """
    t = time.time()
    print("getting positions")
    positions = {}
    contig_length = {}
    with open(cmap, 'r') as f:
        all_inf = f.read()
    print("read file")
    split_lines = all_inf.split("\n")
    del all_inf
    del split_lines[-1]
    for line in split_lines:
        if line.startswith('#'):
            continue
        content = line.strip().split('\t')
        if content[0] not in positions:
            positions[content[0]] = []
        pos = int(round(float(content[5]), 0))
        contig_length[content[0]] = int(round(float(content[1]), 0))
        positions[content[0]].append(pos)
    print("positions obtained, took: {}".format(time.time() -t))
    return positions, contig_length


def make_length_table(contig_length, db_name):
    """
    Function that creates a sql table that contains the length of each cmap
    :param contig_length: Dictionary with format: {'contig/cmap name' : length of the cmap}
    :param db_name: Name of the sql database
    :return: None
    """
    connection = sqlite3.connect("{}.db".format(db_name))
    cursor = connection.cursor()

    create_len_table = """
    CREATE TABLE cmap_len (
    id VARCHAR PRIMARY KEY,
    contig VARCHAR(255) NOT NULL,
    length INT NOT NULL
    );"""

    cursor.execute(create_len_table)
    i = 1
    for contig in contig_length:
        pop_len_table = """
        INSERT INTO cmap_len VALUES(
        {},{},{});""".format(i, contig, contig_length[contig])

        cursor.execute(pop_len_table)
        connection.commit()
        i += 1
