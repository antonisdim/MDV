from Bio import SeqIO
from sys import argv

script = argv[0]
regions = argv[1]
input_files = argv[2:]

files_to_keep = []

for file in input_files:

	record_nucleotide_dict = {}

	for record in SeqIO.parse(file, "fasta"):
		A = record.seq.count('A')
		C = record.seq.count('C')
		T = record.seq.count('T') 
		G = record.seq.count('G') 
		bases_covered = A + C + T + G
		record_nucleotide_dict[record.id] = bases_covered / float(len(record.seq))
		print(record.id, len(record.seq), sep=' ')
		print('A:', A, 'T:', T, 'C:', C, 'G:', G)

	print(record_nucleotide_dict)

	if all(value >= 0.10 for value in record_nucleotide_dict.values()):
		files_to_keep.append(file)

if regions == 'coding':
    with open('selected_orfs_coding.txt', 'w') as fout:
	    print("\n".join(files_to_keep), file=fout)
elif regions == 'non_coding':
    with open('selected_orfs_non_coding.txt', 'w') as fout:
            print("\n".join(files_to_keep), file=fout)

print(len(files_to_keep))
