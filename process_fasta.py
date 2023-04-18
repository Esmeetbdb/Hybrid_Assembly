import time
def get_contigs(fasta):
    """
    This function extracts all contigs from a fasta file and puts them in a dictionary. For each fasta contig a
    reverse is also made
    :param fasta: string that contains the path to the fasta file to be used during alignment
    :return: contigs: dictionary with format {'contig name' : 'sequence information'}
    """
    t = time.time()
    contigs = {}
    with open(fasta, 'r') as f:
        sequence = f.read()

    split_contigs = sequence.split('>')
    del sequence
    del split_contigs[0]

    for contig in split_contigs:
        content = contig.split("\n", 1)
        contigs[content[0].strip()] = content[1].replace("\n", "").upper()
    temp_reverse = {}
    for contig in contigs:
        temp_reverse['{}_reverse'.format(contig)] = contigs[contig][::-1]
    contigs.update(temp_reverse)
    print('Finished obtaining contigs, took: {}'.format(time.time()-t))
    return contigs


def get_reverse_complement(string):
    """
    Function that creates the complementary sequence of a DNA sequence
    :param DNA_seq: String that contains the DNA sequence of which the complement should be made.
    :return: String that contains the complement of input sequence
    """
    reverse_complement = ""
    for char in string:
        if char == "G":
            reverse_complement += "C"
        elif char == "C":
            reverse_complement += "G"
        elif char == "A":
            reverse_complement += "T"
        elif char == "T":
            reverse_complement += "A"
        else:
            print("Enzyme site incorrect. Please try again.")
    return reverse_complement

def get_recognition_site(contigs, enzyme_site):
    """
    Function that finds the positions of all enzyme recognition sites per fasta contig
    :param contigs: Dictionary with format {'contig name' : 'sequence information'}
    :param enzyme_site: String that contains the DNA sequence recognised by the Enzyme used in cmap creation.
    :return: Dictionary with format: {'contig/cmap name' : [ordered list of enzyme site locations]}
    """
    t = time.time()
    sites = {}
    len_enzyme = len(enzyme_site)
    reverse_enzyme = get_reverse_complement(enzyme_site)
    for seq in contigs:
        sites[seq] = []
        for i in range(len(contigs[seq])):
            if contigs[seq][i:i + len_enzyme] == enzyme_site or contigs[seq][i:i + len_enzyme] == reverse_enzyme[::-1]:
                sites[seq].append(i + 1)
    print("Finished get_recognition_sites, took: {}".format(time.time()-t))
    return sites


