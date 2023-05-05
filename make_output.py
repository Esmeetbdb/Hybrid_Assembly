def make_fasta(prefix, seq_dict):
    with open('{}.fa'.format(prefix), 'w') as file:
        for cmap in seq_dict:
            file.writelines('>{}\n'.format(cmap))
            file.writelines('{}\n'.format(seq_dict[cmap]))


def counter_file(prefix, mapped_positions_dictionary):
    with open('{}.mapped_pos.csv'.format(prefix), 'w') as file:
        for round in mapped_positions_dictionary:
            file.writelines("{}, {}\n".format(round, mapped_positions_dictionary[round]))

