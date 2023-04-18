import sqlite3
from multiprocessing import Process
import argparse
import time
import sys

def cmap2sql(cmap, db_name, k):
    import process_cmap as pc
    import make_sql
    positions_cmap, contig_length = pc.get_positions(cmap)
    print(len(positions_cmap))
    pc.make_length_table(contig_length, db_name)
    distances_cmap,distance_list_cmap = make_sql.get_distance(positions_cmap)
    make_sql.make_db(distances_cmap, db_name, k, positions_cmap, 'cmap')


def fa2sql(fasta, enzyme_site, db_name, k):
    import process_fasta as pf
    import make_sql

    contigs_fasta = pf.get_contigs(fasta)
    positions_fasta = pf.get_recognition_site(contigs_fasta, enzyme_site)
    distances_fasta, distance_list_fasta = make_sql.get_distance(positions_fasta)
    make_sql.make_db(distances_fasta, db_name, k, positions_fasta, "fasta")

    return contigs_fasta

def assemble(args):
    cmap2sql(args.cmap, args.prefix, args.k_mer)
    fa2sql(args.Fasta, args.Enzyme_site, args.prefix, args.k_mer)

    import map_fasta_parallel as mf
    import remove_pos as rp
    import find_overlaps_np_parallel as fo
    import process_fasta as pf
    import make_output as mo

    rounds_list = args.overlap_len.split(',')
    print(rounds_list)
    deviation_list = args.deviationlist.split(',')
    fasta_seq_dict = pf.get_contigs(args.Fasta)
    seq_dict = {}

    for i in range(len(rounds_list)):
        print("Round: {}".format(i))
        fo.find_overlaps(args.prefix, int(deviation_list[i]), args.n_threads)
        seq_dict, remove_fa, remove_cmap_list = mf.merge(args.prefix, args.Enzyme_site, fasta_seq_dict, rounds_list[i], args.n_threads, seq_dict)
        print("Number of fasta locations mapped: {}".format(len(remove_fa)))
        rp.remove_fasta(remove_fa, args.prefix)
        rp.remove_cmap(remove_cmap_list, args.prefix)
    mo.make_fasta(args.prefix, seq_dict)



def main():
    parser = argparse.ArgumentParser(description = "Make a de novo assembly by combining Optical mapping data with long read sequencing data.")


    parser.add_argument("Fasta", type = str, help = "Path to fasta file that contains long read sequencing data.")
    parser.add_argument("cmap", type = str, help = "Path to cmap file with optical mapping data.")
    parser.add_argument("Enzyme_site", type = str, help = "Recognition site of the enzyme using in optical mapping.")
    parser.add_argument("prefix", type=str, help="Prefix of the final output file. Will be in fasta format.")
    parser.add_argument("--n_threads", "-c", type=int, default = 16, help="Number of threads to use during mapping.")
    parser.add_argument("--database_name","-d", type = str, default = "Hybrid_assembly", help = "Name of the sql database that will be used during assembly step. Must be unique.")
    parser.add_argument("--k_mer", "-k", type = int, default = 5, help = "K-mer size used for overlapping fasta and optical mapping sites. Larger genomes will require larger k-mer sizes. Larger k-mer sizes lower sensitivity.")
    parser.add_argument("--deviationlist", "-D", type=str, default="500,1000,2000", help = "The max deviation applied to find overlaps for each round of alignment")
    parser.add_argument("--overlap_len", "-o", type=str, default="6,8,10", help = "The number of positions that have to overlap at least for an alignment to be considered correct")
    parser.set_defaults(func=assemble)
    args = parser.parse_args(sys.argv[1:])
    args.func(args)

if __name__ == '__main__':
    main()
