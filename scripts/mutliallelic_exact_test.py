import pandas as pd 
from sys import argv 
from scipy.stats import fisher_exact

script, snp_table = argv

# from https://link.springer.com/article/10.1007%2Fs11262-007-0157-1
# (UL) region is 113,610 bp in length and extends from positions 14,696 to 128,305
# (US) region is 11,668 bp in length and extends from positions 154,349 to 166,016

# the combined length of the IRL, a-like sequence and IRS regions is 26,133 bp
# The lengths of the terminal long (TRL) and short (TRS) repeats (13,903 and 12,184 bp, respectively) 


# total length
total_length = 178246 

# coordinates of the unique regions
unique_long = range(14696, 128305)
unique_short = range(154349, 166016)

# total length of the unique and duplicated regions
unique_regions_bp = len(unique_long) + 1 + len(unique_short) + 1
dup_regions_bp = total_length - unique_regions_bp

# counters
multiallele_unique = 0
multiallele_dup = 0

# read table with mutliallelic sites
mutliallelic_sites = pd.read_csv(snp_table, sep='\t')


for site in mutliallelic_sites['POS'].unique():
    # if in unique regions counts add it to the counter
    if site in unique_long or site in unique_short:
        multiallele_unique += 1
    # if not in there it means it is in the duplicated regions, added there
    else:
        multiallele_dup += 1 

oddsratio, pvalue = fisher_exact([[multiallele_unique, multiallele_dup], [unique_regions_bp, dup_regions_bp]])


print(f'Mutliallelic sites found in unique regions of the MDV genome: {multiallele_unique}')
print(f'Mutliallelic sites found in duplicated regions of the MDV genome: {multiallele_dup}')
print(f'Oddsratio of Fisher exact test: {oddsratio}')
print(f'Pvalue of Fisher exact test: {pvalue}')
