#!/usr/bin/env python
# -*- coding: utf-8 -*

__author__ = "Evangelos A. Dimopoulos"
__copyright__ = "Copyright 2021, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

import pandas as pd


def fix_bed_coord(input_file, output_file):
    """Convert 0 indexed bed coords into 1 indexed"""

    bed = pd.read_csv(input_file, sep="\t", names=["Chrom", "Start", "End"])

    # we add 1 so the starting coord is 1 indexed

    bed["Start"] = bed["Start"] + 1

    # we subtract 1 because the intervals from bedtools are half closed

    bed["End"] = bed["End"] - 1

    bed["new_Chrom"] = [
        "_".join(["Intergenic", str(num + 1)]) for num in range(len(bed))
    ]

    bed[["new_Chrom", "Start", "End"]].to_csv(
        output_file, sep="\t", index=False, header=False
    )


if __name__ == "__main__":
    # noinspection PyUnresolvedReferences
    fix_bed_coord(
        input_file=snakemake.input[0],
        output_file=snakemake.output[0],
    )
