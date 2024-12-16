#!/usr/bin/env python3
"""
Name: Liz Stephany Villabona Arenas
Exercise: Predict
Usage: python3 TMHMM.py filename
    where: filename: str, path to file in .3line format
"""
from sys import argv

def store_data(filename):
    """
    Parse the .3line format and returns a list of data

    :param filename: str, path to file in .3line format
    :return: list of data
    """
    with open(filename, 'r') as f:
        lines = f.readlines()
    list_data=[]
    for line in range(0, len(lines),3):
        list_line = lines[line:line + 3]
        label = list_line[0]
        protein_seq = list_line[1]
        symbol_class = list_line[2]
        list_data.append([label, protein_seq, symbol_class])
    return list_data

def number_records(list_data):
    """
    Calculate number of records in the list of data

    :param list_data: list of tuples
    :return: tuple, number_seq and transmembrane records
    """
    number_seq = len(list_data)
    trans = 0
    for labels in list_data[:]:
        pos = labels[0].count('M')
        if pos != 0:
            trans += 1
    return number_seq, trans

def count_transmembrane_segments(symbol_class):
    """
    Calculate proteins classified as transmembrane for a segment

    :param symbol_class: str, protein symbol class
    :return: int, number of proteins classified as transmembrane
    """
    symbol_class = symbol_class.replace('O', 'I')
    segments = symbol_class.split('I')
    transmembrane_count = 0
    for seg in segments:
        if 'M' in seg:
            transmembrane_count += 1
    return transmembrane_count

def proteins_transmembrane(list_data):
    """
    Calculate proteins classified as transmembrane and return a TSV file

    :param list_data: list of tuples
    :return: str, path to file in TSV format
    """
    with open('TMHMM.tsv', 'w') as output_file:
        for data in list_data:
            protein_id = data[0].split('|')[0].lstrip('>').strip()
            transmembrane = count_transmembrane_segments(data[2])
            line = (f"{protein_id}\t{transmembrane}\n")
            output_file.write(line)

def main():
    """ Main function of the script """
    filename = argv[1]

    list_data = store_data(filename)
    number_seq, trans = number_records(list_data)
    print(number_seq)
    print(trans)
    proteins_transmembrane(list_data)

if __name__ == "__main__":
    main()