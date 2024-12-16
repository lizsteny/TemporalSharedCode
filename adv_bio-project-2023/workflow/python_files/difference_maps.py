#!usr/bin/env python3
""" compares TAIR10 and Araport11 gene mapping results from tx2gene TSV files
Author: Mathijs Balk
Usage: python3 difference_maps.py [TAIR10 file] [Araport11 file]
With both files being a tx2gene TSV file
"""
from sys import argv
from pathlib import Path

def parse_geneids(path_tx2gene_file):
    """ opens a tx2gene TSV file and parses out the gene IDs

    :param path_tx2gene_file: path to the file to be scanned
    :return: set containing the found gene IDs
    """
    geneid_list = []
    # open the tx2gene file
    with open(path_tx2gene_file, 'r') as file:
        for line in file:
            # parse the gene ID out of the second column
            geneid = line.strip().split()[1]
            # put it into a list
            geneid_list.append(geneid)
    # return the found gene IDs as a set
    return set(geneid_list)


def main():
    """ takes in two tx2gene files, compares gene IDs, prints the unique ones

    :return: None
    """
    if len(argv) != 2:
        raise ValueError('use python3 '
              'difference_maps.py [tx2gene file directory]')

    # make paths to the tx2gene TSV files
    path_tair10 = Path(argv[1] + 'tx2geneTAIR10.tsv')
    path_araport11 = Path(argv[1] + 'tx2geneAraport11.tsv')

    # parse the gene IDs out of the tx2gene files
    tair10_set = parse_geneids(path_tair10)
    araport11_set = parse_geneids(path_araport11)

    # see which gene IDs are unique to one of the sets
    difference_set = tair10_set ^ araport11_set

    # print every unique gene ID to the screen
    for gene_id in difference_set:
        print(gene_id)

    # print the amount of unique gene IDs
    print(f'amount of unique gene IDs: {len(difference_set)}')

# make sure the script is run as a whole
if __name__ == "__main__":
    main()
