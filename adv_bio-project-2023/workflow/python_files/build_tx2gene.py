# !/usr/bin/env python3
"""
Author: Jeroen Persoon, Mathijs Balk adapted for Snakemake by Rik Aikes
Description: builds a tx2gene.tsv file for the specified genome. Used by
tximport in R.
Usage: Snakemake using script rule
"""
from pathlib import Path


def select_TSVfiles():
    """

    :return: paths to TSV files
    """
    group_dir = '/prog/BIF30806/project/groen_team/'

    # find samples in TAIR10 dir, put in list
    TAIR10_path = Path(f'{group_dir}results_TAIR10/')
    TAIR10_sampledirs = list(TAIR10_path.iterdir())

    # add abundance.tsv to every path
    for sample_dir in TAIR10_sampledirs:
        tsv_path = Path(f'{sample_dir}/abundance.tsv')
        yield tsv_path

def parser(tsv_path):
    """makes a new file with 2 columns: first column contains the
    transcript_id and second column (tab delimited) contains the gen_id.
    input_file: the kallisto output file
    """
    transcript_ids = []
    with open(tsv_path, 'r') as fo:
        for line in fo:

            # skip empty lines
            if not line.strip():
                continue

            # takes all input file line starting with transcript id
            elif line.startswith('AT'):
                transcript_id = line.split('\t')[0]
                transcript_ids.append(transcript_id)

        # return the transcript IDs as a set
        return set(transcript_ids)



def write_file(complete_set, outf):
    transcript_id_list = list(complete_set)

    # make list with gene IDs
    gene_id_list = []
    for transcript_id in transcript_id_list:
        # throw away everything after the dot
        gene_id = transcript_id.split('.')[0]
        gene_id_list.append(gene_id)

    # make path to output file
    output_file = Path(outf)

    # write in a tab-delimited file
    with open(output_file, 'w') as pf: # open the output file

        # write header based on tximport manual
        pf.write('TXNAME\tGENEID\n')
        pf.write('<chr>\t<chr>\n')

        # write every transcriptID and its geneID
        for pos in range(len(transcript_id_list)):
            sentence = '{}\t{}\n'\
                .format(transcript_id_list[pos], gene_id_list[pos])
            pf.write(sentence)


def main():
    # get output parameter from snakemake
    try:
        outf=snakemake.output[0]
        files= snakemake.input
        complete_set = set()
        for tsv_path in files:
            set_transcript_ids = parser(tsv_path)
            complete_set = complete_set | set_transcript_ids
    except:
        outf= 'build_tx2gene_text.txt'
        complete_set = set()
        for tsv_path in select_TSVfiles():
            set_transcript_ids = parser(tsv_path)
            complete_set = complete_set | set_transcript_ids
    write_file(complete_set, outf)

if __name__ == '__main__':
    main()
