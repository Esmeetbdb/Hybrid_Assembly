def make_fasta(prefix, seq_dict):
    with open('{}.fa'.format(prefix), 'w') as file:
        for cmap in seq_dict:
            file.writelines('>{}\n'.format(cmap))
            file.writelines('{}\n'.format(seq_dict[cmap]))
