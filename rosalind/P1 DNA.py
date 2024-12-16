#!/usr/bin/env python3
"""
Name: Liz Stephany Villabona Arenas
Exersice reading data and basic operations on DNA strings
"""
from sys import argv
def GC_content(DNA_string):
    """
    Return GC content as a percentage
    DNA_string: str, DNA sequence, only ACGT characters allowed
    """
    GC_count = DNA_string.count('G') + DNA_string.count('C')
    GC_percentage = (GC_count / len(DNA_string)) * 100
    return GC_percentage
def trans_DNA_RNA(DNA_string):
    """
    Return the transcription of the string
    DNA_string: str, DNA sequence, only ACGT characters allowed
    """
    RNA_string = ""
    for nucleotide in DNA_string:
        if nucleotide == 'C':
            RNA_string += 'C'
        elif nucleotide == 'G':
            RNA_string += 'G'
        elif nucleotide == 'A':
            RNA_string += 'A'
        elif nucleotide == 'T':
            RNA_string += 'U'
        else:
            print('Please use a DNA sequence, only ACGT characters are allowed')
    return RNA_string
def reverse_com_DNA(DNA_string):
    """
    Return the reverse complement of the string
    DNA_string: str, DNA sequence, only ACGT characters allowed
    """
    reverse_string = ""
    for nucleotide in reversed(DNA_string):
        if nucleotide == 'A':
            reverse_string += 'T'
        elif nucleotide == 'T':
            reverse_string += 'A'
        elif nucleotide == 'C':
            reverse_string += 'G'
        elif nucleotide == 'G':
            reverse_string += 'C'
        else:
            print('Please use a DNA sequence, only ACGT characters are allowed')
    return reverse_string
def reverse_trans_DNA_RNA(DNA_string):
    """
    Return the reverse complement transcription of the string
    DNA_string: str, DNA sequence, only ACGT characters allowed
    """
    reverse_trans_string = ""
    for nucleotide in reversed(DNA_string):
        if nucleotide == 'A':
            reverse_trans_string += 'U'
        elif nucleotide == 'T':
            reverse_trans_string += 'A'
        elif nucleotide == 'C':
            reverse_trans_string += 'G'
        elif nucleotide == 'G':
            reverse_trans_string += 'C'
        else:
            print('Please use a DNA sequence, only ACGT characters are allowed')
    return reverse_trans_string
def main():
    """
    Main code for the script
    """
    DNA_string = argv[1]
    print(GC_content(DNA_string))
    print(trans_DNA_RNA(DNA_string))
    print(reverse_com_DNA(DNA_string))
    print(reverse_trans_DNA_RNA(DNA_string))
if __name__ == "__main__":
    main()