import sqlite3
def remove_fasta(remove_fa_list, db_name):
    connection = sqlite3.connect("{}.db".format(db_name))
    cursor = connection.cursor()

    remove_cmd = """
        DELETE FROM fasta
        WHERE contig=? AND pos_contig=?"""
    cursor.executemany(remove_cmd,remove_fa_list)
    connection.commit()

def remove_cmap(remove_cmap_list, db_name):
    connection = sqlite3.connect("{}.db".format(db_name))
    cursor = connection.cursor()

    remove_cmd = """
            DELETE FROM cmap
            WHERE contig=? AND pos_contig=?"""
    cursor.executemany(remove_cmd, remove_cmap_list)
    connection.commit()
