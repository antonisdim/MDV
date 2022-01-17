#!/usr/bin/env python
# -*- coding: utf-8 -*

__author__ = "Evangelos A. Dimopoulos"
__copyright__ = "Copyright 2021, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

from Bio import SeqIO


def clean_from_BACs(input_file, output_file):
    """Clean residual BAC sequences from modern genomes"""

    fasta_file = open(input_file, "rU")

    record_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))

    if "KT833851.1" in record_dict:
        record_dict["KT833851.1"].seq = (
            record_dict["KT833851.1"].seq[:157057]
            + record_dict["KT833851.1"].seq[163447:]
        )
    elif "KT833852.1" in record_dict:
        record_dict["KT833852.1"].seq = (
            record_dict["KT833852.1"].seq[:157057]
            + record_dict["KT833852.1"].seq[163447:]
        )
    elif "FJ436097.1" in record_dict:
        record_dict["FJ436097.1"].seq = (
            record_dict["FJ436097.1"].seq[:156920]
            + record_dict["FJ436097.1"].seq[164269:]
        )

    SeqIO.write(record_dict.values(), output_file, "fasta")


if __name__ == "__main__":
    # noinspection PyUnresolvedReferences
    clean_from_BACs(
        input_file=snakemake.input[0],
        output_file=snakemake.output[0],
    )
