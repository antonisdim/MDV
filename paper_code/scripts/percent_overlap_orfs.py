#!/usr/bin/env python
# -*- coding: utf-8 -*

__author__ = "Evangelos A. Dimopoulos"
__copyright__ = "Copyright 2021, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

from Bio import SeqIO


def percent_overlap_orfs(gene_alns, output_file):
    """Count the percent overlap in each gene alignment"""

    files_to_keep = []

    for file in gene_alns:

        record_nucleotide_dict = {}

        for record in SeqIO.parse(file, "fasta"):
            A = record.seq.count("A")
            C = record.seq.count("C")
            T = record.seq.count("T")
            G = record.seq.count("G")
            bases_covered = A + C + T + G
            record_nucleotide_dict[record.id] = bases_covered / float(len(record.seq))
            print(record.id, len(record.seq), sep=" ")
            print("A:", A, "T:", T, "C:", C, "G:", G)

        print(record_nucleotide_dict)

        if all(value >= 0.10 for value in record_nucleotide_dict.values()):
            files_to_keep.append(file)

    with open(output_file, "w") as outfile:
        print("\n".join(files_to_keep), file=outfile)


if __name__ == "__main__":
    # noinspection PyUnresolvedReferences
    percent_overlap_orfs(
        gene_alns=snakemake.input,
        output_file=snakemake.output[0],
    )
