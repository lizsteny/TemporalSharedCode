#!/usr/bin/env python3
"""
Name: Liz Stephany Villabona Arenas
Exercise: Translate RNA to protein
"""
from sys import argv
def translate_RNA_to_protein(RNA_string):
    """
    Codon table
    RNA: str, RNA sequence, only ACGU characters allowed
    """
    codon_table = {
        "UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L",
        "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
        "AUU": "I", "AUC": "I", "AUA": "I", "AUG": "M",
        "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
        "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S",
        "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
        "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        "UAU": "Y", "UAC": "Y", "UAA": "*", "UAG": "*",
        "CAU": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
        "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K",
        "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",
        "UGU": "C", "UGC": "C", "UGA": "*", "UGG": "W",
        "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R",
        "AGU": "S", "AGC": "S", "AGA": "R", "AGG": "R",
        "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    }
    protein = ""
    stop_codon = 0
    if RNA_string[:3] != 'AUG':
        # First codon is not the start codon
        return None
    if RNA_string[-3:] not in ['UGA','UAG','UAA']:
        # Last codon is not a stop codon
        return None
    if len(RNA_string) % 3 != 0:
        # String not a multiple of 3
        return None
    for nucleotide in range(0, len(RNA_string), 3):
        codon = RNA_string[nucleotide:nucleotide + 3]
        if codon in codon_table:
            amino_acid = codon_table[codon]
            if amino_acid == '*':
                stop_codon += 1
                if stop_codon > 1:
                    return None
            else:
                protein += amino_acid
        else:
            # Codon not found in codon table
            return None
    return protein
def main():
    """
     Main code for the script
     """
    RNA_strings = argv[1:]
    for RNA_string in RNA_strings:
        protein_seq = translate_RNA_to_protein(RNA_string)
        if protein_seq is not None:
            print(protein_seq)
if __name__ == "__main__":
    main()