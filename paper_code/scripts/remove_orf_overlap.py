#!/usr/bin/env python
# -*- coding: utf-8 -*

__author__ = "Evangelos A. Dimopoulos"
__copyright__ = "Copyright 2021, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

import pandas as pd
import numpy as np


def remove_orf_overlap(input_file, output_file):
    """Remove overlapping regions among ORFs"""

    pd.set_option("precision", 0)

    # read bed file

    bed_raw = pd.read_csv(
        filepath_or_buffer=input_file,
        sep="\t",
        names=["chrom", "start", "end", "gene", "something", "strand"],
    )

    # define new cols
    bed_file = bed_raw.copy()

    bed_file["new_start"] = np.nan
    bed_file["overlap"] = np.nan
    bed_file["new_start"] = np.nan
    bed_file["new_end"] = np.nan
    bed_file.sort_values(by=["start"], inplace=True)

    # give new coordinates
    for idx, row in bed_file.iterrows():
        if idx == 0:
            print("first line, not useful")
            continue
        else:
            if bed_file.iloc[idx, 1] < bed_file.iloc[idx - 1, 2]:
                print("overlap found", bed_file["gene"].iloc[idx], sep=" ")
                bed_file.loc[idx, "overlap"] = (
                    bed_file.iloc[idx - 1, 2] - bed_file.iloc[idx, 1]
                )
                bed_file.loc[idx, "new_start"] = bed_file.iloc[idx - 1, 2] + 1
                bed_file.loc[idx - 1, "new_end"] = (
                    bed_file.iloc[idx - 1, 2] - bed_file["overlap"].iloc[idx] - 1
                )

    # print(bed_file)

    # define new df to store the new coordinates
    new_bed_file = pd.DataFrame()

    new_bed_file["chrom"] = bed_file.iloc[:, 0]
    new_bed_file["start"] = bed_file["new_start"].fillna(bed_file["start"])
    new_bed_file["end"] = bed_file["new_end"].fillna(bed_file["end"])
    new_bed_file["gene"] = bed_file.iloc[:, 3]
    new_bed_file["something"] = bed_file.iloc[:, 4]
    new_bed_file["strand"] = bed_file.iloc[:, 5]

    for idx, row in new_bed_file.iterrows():
        if new_bed_file.loc[idx, "start"] >= new_bed_file.loc[idx, "end"]:
            new_bed_file.drop([idx], axis=0, inplace=True)

    # make sure the coordinates still give complete codons
    for idx, row in new_bed_file.iterrows():
        length = len(
            range(int(new_bed_file["start"][idx]), int(new_bed_file["end"][idx]) + 1)
        )
        if length % 3 == 0:
            print(
                new_bed_file["gene"][idx],
                "at",
                length,
                "corresponds to complete codons",
                sep=" ",
            )
        elif (length - 1) % 3 == 0:
            new_bed_file.loc[idx, "end"] = new_bed_file["end"][idx] - 1
            print(
                new_bed_file["gene"][idx],
                "at",
                length,
                "bp does not correspond to complete codons, we subtract one bp for complete codons:",
                length - 1,
                sep=" ",
            )
        elif (length - 2) % 3 == 0:
            new_bed_file.loc[idx, "end"] = new_bed_file["end"][idx] - 2
            print(
                new_bed_file["gene"][idx],
                "at",
                length,
                "bp does not correspond to complete codons, we subtract two bp for complete codons:",
                length - 2,
                sep=" ",
            )

    # print(new_bed_file)

    # print the updated bed file
    new_bed_file["start"] = new_bed_file["start"].astype(int)
    new_bed_file["end"] = new_bed_file["end"].astype(int)

    new_bed_file.to_csv(output_file, sep="\t", index=False, header=False)


if __name__ == "__main__":
    # noinspection PyUnresolvedReferences
    remove_orf_overlap(
        input_file=snakemake.input[0],
        output_file=snakemake.output[0],
    )
