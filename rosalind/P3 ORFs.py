#!/usr/bin/env python3
"""
Name: Liz Stephany Villabona Arenas
Exersice to find ORFs in a DNA string and then, translate them.
"""
from sys import argv
def ORFs_finding(DNA_string):
    """
    Function to calculate the ORFs in a DNA sequence
    DNA_string: str, RNA sequence, only ACGT characters allowed
    """
    orfs = []
    for codon in range(3):
        for nucleotide_start in range(codon, len(DNA_string), 3):
            start_codon = DNA_string[nucleotide_start:nucleotide_start + 3]
            if start_codon == "ATG":
                for nucleotide_end in range(nucleotide_start, len(DNA_string), 3):
                    stop_codon = DNA_string[nucleotide_end:nucleotide_end + 3]
                    if stop_codon in ("TAA", "TGA", "TAG"):
                        orf = DNA_string[nucleotide_start:nucleotide_end + 3]
                        orfs.append(orf)
                        break
    return orfs
def translate_DNA_to_protein(DNA_string):
    """
    Codon table
    DNA: str, DNA sequence, only ACGT characters allowed
    """
    codon_table = {
        'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
        'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
        'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
        'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
        'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
        'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
        'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
    }
    protein = []
    for nucleotide in range(0, len(DNA_string), 3):
        codon = DNA_string[nucleotide:nucleotide + 3]
        amino_acid = codon_table.get(codon, '')
        if amino_acid == '*':
            break
        protein.append(amino_acid)

    return ''.join(protein)

def print_translations(DNA_string):
    for orf in ORFs_finding(DNA_string):
        protein_sequence = translate_DNA_to_protein(orf)
        if protein_sequence:
            print(protein_sequence)

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
    return reverse_string

def main():
    """
    Main code for the script
    """
    DNA_string = argv[1]

    forward_strand_orfs = ORFs_finding(DNA_string)
    reverse_string = reverse_com_DNA(DNA_string)
    reverse_strand_orfs = ORFs_finding(reverse_string)

    print("Total forward strand ORFs:", len(forward_strand_orfs))
    print("Total reverse strand ORFs:", len(reverse_strand_orfs))

    print("Forward strand ORFs:")
    print_translations(DNA_string)

    print("\nReverse strand ORFs:")
    print_translations(reverse_string)

if __name__ == "__main__":
    main()