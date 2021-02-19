import pandas as pd 
from sys import argv 

script, input_file, output_file = argv 

df = pd.read_csv(input_file, sep='\t', names=['Chrom', 'Start', 'End'])

# we add 1 so the starting coord is 1 indexed

df['Start'] = df['Start'] + 1

# we subtract 1 because the intervals from bedtools are half closed 

df['End'] = df['End'] - 1

df.to_csv(output_file, sep='\t', index=False, header=False)

