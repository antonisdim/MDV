from Bio import SeqIO
from Bio.Seq import MutableSeq

from sys import argv

# script, reference, sample, output = argv

script, alignment_file, output = argv

def mutation_type(alleles):
    """
    Is this mutation a transition (A <-> G and C <-> T) or a transversion (everything else).
    """
    alleles = set([a.upper() for a in alleles])
    assert len(alleles) == 2 and alleles.issubset({"A", "C", "G", "T", "N"})
    return "ts" if alleles in [{"A", "G"}, {"C", "T"}] else "tv"


# read the files 
alignment = list(SeqIO.parse(alignment_file, 'fasta'))
reference = [seq for seq in alignment if seq.id == 'EF523390.1'][0]

ts_postitions = []

for seq in alignment:
    
    # check they have they same amount of bases
    if len(reference.seq) != len(seq.seq):
        raise RuntimeError(
            "The reference and the pileup sample fasta have different seq lens. Fix that and revisit."
        )
    
    # start a counter for position indexing
    position = 0
    
    for ref_base, test_base in zip(reference.seq, seq.seq):
        # if extension skip 
        if ref_base.upper() == '-' or test_base.upper() == '-':
            position += 1
            continue
            
        # if the same skip 
        if ref_base.upper() == test_base.upper():
            position += 1
            continue
        
        # if mutation count 
        if mutation_type([ref_base, test_base]) == "ts":
            print(seq.id, position, ref_base, test_base)
            # seq.seq = MutableSeq(str(seq.seq))
            # seq.seq[position] = "N"
            ts_postitions.append(position)
        position += 1

ts_postitions = list(set(ts_postitions))

# set all the columns with Ts as N
for seq in alignment:

    for position in ts_postitions:
        print(f"Substituting position {position} for seq {seq.id}")
        seq.seq = MutableSeq(str(seq.seq))
        seq.seq[position] = "N"
    

# output the new fasta file
with open(output, "w") as out_fasta:
    SeqIO.write(alignment, out_fasta, "fasta")
    
    
    
    