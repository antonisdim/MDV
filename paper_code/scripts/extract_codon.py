#!/usr/bin/env python
# -*- coding: utf-8 -*

__author__ = "Evangelos A. Dimopoulos"
__copyright__ = "Copyright 2021, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

from Bio import SeqIO


def extract_codon(input_fasta, codon_position, output_fasta):
    """Get the correct codon position from a fasta alignment"""

    # convert to 0 indexed position for biopython
    codon = int(codon_position) - 1

    # store the extracted codons here
    codon_list = []

    # This extract every third character, starting at codon_position - 1
    for record in SeqIO.parse(input_fasta, "fasta"):
        codon_rec = record
        codon_rec.seq = record.seq[codon::3]
        codon_list.append(codon_rec)

    # write output
    with open(output_fasta, "w") as outfile:
        SeqIO.write(codon_list, outfile, "fasta")


if __name__ == "__main__":
    # noinspection PyUnresolvedReferences
    extract_codon(
        input_fasta=snakemake.input[0],
        codon_position=snakemake.wildcards.codon,
        output_fasta=snakemake.output[0],
    )
