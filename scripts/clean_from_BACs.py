from Bio import SeqIO
from sys import argv 

script, input_file = argv 

fasta_file = open(argv[1], "rU")

record_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))

record_dict["KT833851.1"].seq = record_dict["KT833851.1"].seq[:157057] + record_dict["KT833851.1"].seq[163447:]

record_dict["KT833852.1"].seq = record_dict["KT833852.1"].seq[:157057] + record_dict["KT833852.1"].seq[163447:]

record_dict["FJ436097.1"].seq = record_dict["FJ436097.1"].seq[:156920] + record_dict["FJ436097.1"].seq[164269:]

SeqIO.write(record_dict.values(), "all_modern_plus_HVT_clean.fasta", "fasta")
