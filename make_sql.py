import sqlite3
import time


def get_distance(sites):
    """
    This function loops over all contigs/cmaps and calculates the distances between each consecutive enzyme site.
    :param sites: dictionary with format: {'contig/cmap name' : [ordered list of enzyme site locations]}
    :return: distances: dictionary with format {'congtig/cmap name' : [list of distances (ordered)]}
    :return: distance_list: contains all distances in one list
    """
    distances = {}
    distance_list = []
    for seq in sites:
        distances[seq] = []
        for i in range(len(sites[seq]) - 1):
            distance = sites[seq][i + 1] - sites[seq][i]
            distances[seq].append(distance)
            distance_list.append(distance)
    print("Finished get_distances")
    return distances, distance_list


def kmer_table(cursor, distance_columns, table_type, index_column_count, col_number, table_content):
    """
    This function creates the table that contains all distances between enzyme sites on the cmap/fasta contig.
    The table has 3 standard columns that contain the ID (primary key of the entry), contig (name of the cmap/fasta),
    and the pos_contig (which positions all the distances in this row are linked to), there is also one column for
    each of the distances per row, the number of distances depends on the chosen k-mer length,
    where number distances = k-mer length
    :param cursor: connection to the sql database
    :param distance_columns: list with the names of the distance columns
    :param table_type: fasta or cmap
    :param index_column_count: list with the names of the distance columns to index the table
    :param col_number: string with a ? for each column in the table
    :param table_content: The actual information to fill up the table after table creation
    :return: None
    """
    create_kmer_table = """
        CREATE TABLE {} (
        id VARCHAR PRIMARY KEY,
        contig VARCHAR(255) NOT NULL,
        pos_contig INT NOT NULL,
        {}
        );
        """.format(table_type, distance_columns)
    cursor.execute(create_kmer_table)
    cursor.executemany('INSERT INTO {} VALUES ({})'.format(table_type, col_number), table_content)
   
    create_index = """
        CREATE INDEX idx_distances_{}
        ON {} ({});""".format(table_type, table_type, index_column_count)
    cursor.execute(create_index)


def contig_table(cursor, position_columns, table_type, col_number, table_content):
    """
    This function creates the table that contains all positions of enzyme sites on the cmap/fasta contig.
    The table has 3 standard columns that contain the ID (primary key of the entry), contig (name of the cmap/fasta),
    and the pos_contig (which indicate that this row is the nth pos for this contig), there is also one column for
    each of the positions per row, the number of distances depends on the chosen k-mer length,
    where number distances = k-mer length + 1
    :param cursor: connection to the sql database
    :param position_columns: list with the names of the position columns
    :param table_type: fasta or cmap
    :param col_number: string with a ? for each column in the table
    :param table_content: The actual information to fill up the table after table creation
    :return: None
    """
    create_contigs_table = """
        CREATE TABLE {}_contigs (
        id INT PRIMARY KEY,
        contig VARCHAR(255) NOT NULL,
        pos_contig INT NOT NULL,
        {}
        );
        """.format(table_type, position_columns)
    cursor.execute(create_contigs_table)

    cursor.executemany('INSERT INTO {}_contigs VALUES (?,{})'.format(table_type, col_number), table_content)


def get_column_count(n_distances):
    """
    This function prepares all the column count related strings and lists for sql table creation
    :param n_distances: gives the number of distances per row
    :return: distance_columns: list with the names of the distance columns
    :return: index_column_count: list with the names of the distance columns to index the table
    :return: position_columns: list with the names of the position columns
    :return: col_number: string with a ? for each column in the table
    """
    distance_columns = ''
    index_column_count = ''
    position_columns = ''
    col_number = '?,?,?'

    for i in range(n_distances):
        distance_columns += 'distance_{} INT(255) NOT NULL,'.format(i+1)
        index_column_count += 'distance_{},'.format(i+1)
    index_column_count = index_column_count[:-1]
    distance_columns = distance_columns[:-1]

    for j in range(n_distances + 1):
        position_columns += 'position_{} INT NOT NULL,'.format(j+1)
    position_columns = position_columns[:-1]

    for i in range(n_distances):
        col_number += ',?'

    return distance_columns, index_column_count, position_columns, col_number


def get_table_content(distance_list, n_distances, id_counter, sites, contig, distances):
    """
    This function generates the input for the distance and position sql tables.
    :param distance_list: A list of ordered distances between consecutive enzyme sites for 1 contig
    :param n_distances: number of distances per row
    :param id_counter: counter to have a unique ID per table entry
    :param sites: dictionary with format: {'contig/cmap name' : [ordered list of enzyme site locations]}
    :param contig: name of the contig
    :param distances: dictionary with format {'congtig/cmap name' : [list of distances (ordered)]}
    :return: values: list
    """
    k = 1
    values = []
    positions = []
    for i in range(len(distance_list) - n_distances):
        values_temp = [str(id_counter), contig, str(k)]
        positions_temp = [str(id_counter), contig, str(k)]
        k += 1
        id_counter += 1

        for j in range(n_distances+1):
            values_temp.append(str(distances[contig][i + j]))
            positions_temp.append(str(sites[contig][i + j]))
        del (values_temp[-1])
        values.append(values_temp)
        positions.append(positions_temp)
    return values, positions, id_counter


def make_db(distances, db_name, n_distances, sites, table_type):
    """
    This function calls all functions required to initiate the sql database and connects to the database to enter
    all table information.
    :param distances: dictionary with format {'congtig/cmap name' : [list of distances (ordered)]}
    :param db_name: name of the sql database
    :param n_distances: number of distances per row
    :param sites: dictionary with format: {'contig/cmap name' : [ordered list of enzyme site locations]}
    :param table_type: fastq ir cmap
    :return: None
    """
    t = time.time()
    print("Start making db, {}".format(table_type))
    connection = sqlite3.connect("{}.db".format(db_name))
    cursor = connection.cursor()
    id_counter = 1

    distance_columns, index_column_count, position_columns, col_number = get_column_count(n_distances)

    pop_table = []
    pop_contigs = []
    for contig in distances:
        distance_list = distances[contig]
        values, positions, id_counter = get_table_content(distance_list, n_distances, id_counter,
                                                          sites, contig, distances)
        pop_table += values
        pop_contigs += positions
    kmer_table(cursor, distance_columns, table_type, index_column_count, col_number, pop_table)
    contig_table(cursor, position_columns, table_type, col_number, pop_contigs)
    connection.commit()
    print("Finished populating database, took; {}".format(time.time()-t))

