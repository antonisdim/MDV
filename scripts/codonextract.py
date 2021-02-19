#!/usr/bin/env python

### <-- codonextract.py made by Audrey T. Lin (audrey.lin@gtc.ox.ac.uk) --> ###
### <-- Department of Zoology, University of Oxford --> ###

from Bio import SeqIO
from sys import argv

# python scripts/codonextract.py full_CDS.fasta > codon_pos_1_2_or_3.fasta 
script, aln_file, codon_pos = argv 

# convert to 0 indexed position for biopython 
position = int(codon_pos) - 1

# 0 for first codon, 1 for second codon, 2 for third codon
which_codon = position 

for record in SeqIO.parse(aln_file, "fasta"):
    print('>'+record.id)
    # 3 is the codon length. This extract every third character, starting at which_codon
    print(record.seq[which_codon::3]) 

