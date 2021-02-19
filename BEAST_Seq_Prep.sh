### BEAST SEQUENCE PREP STEPS ###

# 1) Clean the sequences from the BACs like we did in step 4) previously 

python scripts/clean_from_BACs.py mod_gen_files/all_modern_plus_HVT.fa

mv all_modern_plus_HVT_clean.fasta mod_gen_files/all_modern_plus_HVT_clean.fa

# 2) Use pyfasta to split the modern genomes into 500 bp kmers with an overlap of 5 bp and align them to the masked genome separately. 


# split it with Biopython SeqIO.to_dict and SeqIO.parse and SeqIO.write

parallel --jobs 4 "pyfasta split -n 1 -k 500 -o 5 {1}" ::: mod_gen_files/AF147806.2.fasta mod_gen_files/AF243438.1.fasta mod_gen_files/AY510475.1.fasta mod_gen_files/DQ530348.1.fasta mod_gen_files/EF523390.1.fasta mod_gen_files/EU499381.1.fasta mod_gen_files/FJ436096.1.fasta mod_gen_files/FJ436097.1.fasta mod_gen_files/JF742597.1.fasta mod_gen_files/JQ314003.1.fasta mod_gen_files/JQ806361.1.fasta mod_gen_files/JQ806362.1.fasta mod_gen_files/JQ809691.1.fasta mod_gen_files/JQ809692.1.fasta mod_gen_files/JQ820250.1.fasta mod_gen_files/JQ836662.1.fasta mod_gen_files/JX844666.1.fasta mod_gen_files/KT833851.1.fasta mod_gen_files/KT833852.1.fasta mod_gen_files/KU173115.1.fasta mod_gen_files/KU173116.1.fasta mod_gen_files/KU173117.1.fasta mod_gen_files/KU173118.1.fasta mod_gen_files/KU173119.1.fasta mod_gen_files/KU744555.1.fasta mod_gen_files/KU744556.1.fasta mod_gen_files/KU744557.1.fasta mod_gen_files/KU744558.1.fasta mod_gen_files/KU744559.1.fasta mod_gen_files/KU744560.1.fasta mod_gen_files/KU744561.1.fasta mod_gen_files/KX290013.1.fasta mod_gen_files/KX290014.1.fasta mod_gen_files/KX290015.1.fasta mod_gen_files/KX290016.1.fasta mod_gen_files/MF431493.1.fasta mod_gen_files/MF431494.1.fasta mod_gen_files/MF431495.1.fasta mod_gen_files/MF431496.1.fasta mod_gen_files/MG432697.1.fasta mod_gen_files/MG518371.1.fasta mod_gen_files/NC_002229.3.fasta 


bowtie2-build ref_genome/EF523390.1_mask.fasta ref_genome/EF523390.1_mask
             
             
parallel --jobs 5 "bowtie2 -f --very-fast-local -p 5 -x ref_genome/EF523390.1_mask -U mod_gen_files/{1}.500mer.5overlap.fasta | /home/antony/bin/bin/samtools view -Shu > mod_gen_files/{1}.bam" ::: AF147806.2.split AF243438.1.split AY510475.1.split DQ530348.1.split EF523390.1.split EU499381.1.split FJ436096.1.split FJ436097.1.split JF742597.1.split JQ314003.1.split JQ806361.1.split JQ806362.1.split JQ809691.1.split JQ809692.1.split JQ820250.1.split JQ836662.1.split JX844666.1.split KT833851.1.split KT833852.1.split KU173115.1.split KU173116.1.split KU173117.1.split KU173118.1.split KU173119.1.split KU744555.1.split KU744556.1.split KU744557.1.split KU744558.1.split KU744559.1.split KU744560.1.split KU744561.1.split KX290013.1.split KX290014.1.split KX290015.1.split KX290016.1.split MF431493.1.split MF431494.1.split MF431495.1.split MF431496.1.split MG432697.1.split MG518371.1.split NC_002229.3.split 


### HVT: mod_gen_files/NC_002641.1.fasta ###

# 3) Produce fasta pileups with HTSBox like before 

parallel --jobs 20 "/home/antony/bin/bin/samtools sort -o mod_gen_files/{1}_sorted.bam mod_gen_files/{1}.bam" ::: AF147806.2.split AF243438.1.split AY510475.1.split DQ530348.1.split EF523390.1.split EU499381.1.split FJ436096.1.split FJ436097.1.split JF742597.1.split JQ314003.1.split JQ806361.1.split JQ806362.1.split JQ809691.1.split JQ809692.1.split JQ820250.1.split JQ836662.1.split JX844666.1.split KT833851.1.split KT833852.1.split KU173115.1.split KU173116.1.split KU173117.1.split KU173118.1.split KU173119.1.split KU744555.1.split KU744556.1.split KU744557.1.split KU744558.1.split KU744559.1.split KU744560.1.split KU744561.1.split KX290013.1.split KX290014.1.split KX290015.1.split KX290016.1.split MF431493.1.split MF431494.1.split MF431495.1.split MF431496.1.split MG432697.1.split MG518371.1.split NC_002229.3.split 

for i in mod_gen_files/*sorted.bam; do echo $i; /home/alex/bin/htsbox pileup -f ref_genome/EF523390.1_mask.fasta -l15 -M $i > masked_mod_files/`basename $i`.fasta; done 

for i in masked_mod_files/*bam.fasta; do sed -i "s/>EF523390.1/>`basename -s .split_sorted.bam.fasta $i`/g" $i; done
 
for i in masked_mod_files/*bam.fasta; do echo "mv $i masked_mod_files/`basename -s .split_sorted.bam.fasta $i`_masked.fasta"; done


# 4) cat the fastas to get the *clean_masked* files 

# get all modern 
cat masked_mod_files/*_masked.fasta > masked_mod_files/all_modern_minus_HVT_clean_masked.fa
cat masked_mod_files/all_modern_minus_HVT_clean_masked.fa HVT/NC_002641.1.fasta > masked_mod_files/all_modern_plus_HVT_clean_masked.fa

# get all ancient 
cat pileup_fastas/*fasta > pileup_fastas/all_anc_clean_masked.fa

# get total
cat pileup_fastas/all_anc_clean_masked.fa masked_mod_files/all_modern_minus_HVT_clean_masked.fa > all_seqs/all_seqs_minus_HVT_clean_masked.fa
cat all_seqs/all_seqs_minus_HVT_clean_masked.fa HVT/NC_002641.1.fasta > all_seqs/all_seqs_plus_HVT_clean_masked.fa


# 5) Align all sequences with MAFFT and build a RAxML tree


# find how many sequences have Ns
seqtk comp all_seqs/all_seqs_minus_HVT_clean_masked.fa | column -t | awk '{print $1"\t"$9/$2*100}'


OL1099	58.6094
OL1385	17.0871
OL1389	40.2825
OL1986	46.1615
OL2272	35.088
AF147806.2	16.3695
AF243438.1	16.0873
AY510475.1	16.7218
DQ530348.1	15.3995
EF523390.1	15.5089
EU499381.1	15.965
FJ436096.1	15.9656
FJ436097.1	15.9684
JF742597.1	16.8307
JQ314003.1	15.983
JQ806361.1	15.9106
JQ806362.1	16.0368
JQ809691.1	16.0464
JQ809692.1	16.0548
JQ820250.1	16.5025
JQ836662.1	16.9423
JX844666.1	16.0811
KT833851.1	16.1311
KT833852.1	16.4009
KU173115.1	16.4133
KU173116.1	17.0641
KU173117.1	16.2248
KU173118.1	16.31
KU173119.1	16.2652
KU744555.1	16.0542
KU744556.1	16.0969
KU744557.1	15.9891
KU744558.1	15.988
KU744559.1	15.9847
KU744560.1	15.9908
KU744561.1	16.0963
KX290013.1	15.8567
KX290014.1	15.9835
KX290015.1	16.1473
KX290016.1	16.1524
MF431493.1	16.1058
MF431494.1	16.1165
MF431495.1	16.1041
MF431496.1	15.868
MG432697.1	15.9179
MG518371.1	15.7069
NC_002229.3	15.859


# The ancient Sequences Have massive amounts of Ns so they can't align (mafft running for more than a week). So we are doing it in steps. Aligning all the modern ones in one go and then sequentially adding every one of them with the add option 


mafft --maxiterate 1000 --thread 5 --nwildcard masked_mod_files/all_modern_plus_HVT_clean_masked.fa 1> aln_files/all_modern_plus_HVT_clean_masked_aln.fasta 2> aln_files/log-all_modern_plus_HVT_clean_masked_aln.txt &

mafft --maxiterate 1000 --thread 5 --nwildcard masked_mod_files/all_modern_minus_HVT_clean_masked.fa 1> aln_files/all_modern_minus_HVT_clean_masked_aln.fasta 2> aln_files/log-all_modern_minus_HVT_clean_masked_aln.txt &



# 6) Build alignments for MLE trees with HVT as an outgroup


# Start adding the ancient sequences 

# FOR BEAST MLE TREE we will be sequentially adding as many ancient seqs as possible and then we will prune the HVT outgroup

# plus HVT and we have 5 ancient sequences

mafft --maxiterate 1000 --thread 5 --nwildcard --add pileup_fastas/OL1099.fasta aln_files/all_modern_plus_HVT_clean_masked_aln.fasta 1> aln_files/all_modern_plus_HVT_clean_masked_aln_1anc.fasta 2> aln_files/log-all_modern_plus_HVT_clean_masked_aln_1anc.txt 

mafft --maxiterate 1000 --thread 5 --nwildcard --add pileup_fastas/OL1385.fasta aln_files/all_modern_plus_HVT_clean_masked_aln_1anc.fasta 1> aln_files/all_modern_plus_HVT_clean_masked_aln_2anc.fasta 2> aln_files/log-all_modern_plus_HVT_clean_masked_aln_2anc.txt 

mafft --maxiterate 1000 --thread 5 --nwildcard --add pileup_fastas/OL1389.fasta aln_files/all_modern_plus_HVT_clean_masked_aln_2anc.fasta 1> aln_files/all_modern_plus_HVT_clean_masked_aln_3anc.fasta 2> aln_files/log-all_modern_plus_HVT_clean_masked_aln_3anc.txt 

mafft --maxiterate 1000 --thread 5 --nwildcard --add pileup_fastas/OL1986.fasta aln_files/all_modern_plus_HVT_clean_masked_aln_3anc.fasta 1> aln_files/all_modern_plus_HVT_clean_masked_aln_4anc.fasta 2> aln_files/log-all_modern_plus_HVT_clean_masked_aln_4anc.txt 

mafft --maxiterate 1000 --thread 5 --nwildcard --add pileup_fastas/OL2272.fasta aln_files/all_modern_plus_HVT_clean_masked_aln_4anc.fasta 1> aln_files/all_modern_plus_HVT_clean_masked_aln_5anc.fasta 2> aln_files/log-all_modern_plus_HVT_clean_masked_aln_5anc.txt 

cp aln_files/all_modern_plus_HVT_clean_masked_aln_5anc.fasta  aln_files/all_modern_plus_HVT_clean_masked_aln_anc_BT.fasta 


# FOR A MLE TREE WITH AS MANY SEQS AS POSSIBLE 

# I'll continue adding samples so that I can produce one final alignment with a final outgroup for a general RAxML tree 

# we have 11 ancient sequences

mafft --maxiterate 1000 --thread 5 --nwildcard --add pileup_fastas/OL1599.fasta aln_files/all_modern_plus_HVT_clean_masked_aln_5anc.fasta 1> aln_files/all_modern_plus_HVT_clean_masked_aln_6anc.fasta 2> aln_files/log-all_modern_plus_HVT_clean_masked_aln_6anc.txt 

mafft --maxiterate 1000 --thread 5 --nwildcard --add pileup_fastas/OL1934.fasta aln_files/all_modern_plus_HVT_clean_masked_aln_6anc.fasta 1> aln_files/all_modern_plus_HVT_clean_masked_aln_7anc.fasta 2> aln_files/log-all_modern_plus_HVT_clean_masked_aln_7anc.txt 

mafft --maxiterate 1000 --thread 5 --nwildcard --add pileup_fastas/OL1936.fasta aln_files/all_modern_plus_HVT_clean_masked_aln_7anc.fasta 1> aln_files/all_modern_plus_HVT_clean_masked_aln_8anc.fasta 2> aln_files/log-all_modern_plus_HVT_clean_masked_aln_8anc.txt 

mafft --maxiterate 1000 --thread 5 --nwildcard --add pileup_fastas/OL1984.fasta aln_files/all_modern_plus_HVT_clean_masked_aln_8anc.fasta 1> aln_files/all_modern_plus_HVT_clean_masked_aln_9anc.fasta 2> aln_files/log-all_modern_plus_HVT_clean_masked_aln_9anc.txt 

mafft --maxiterate 1000 --thread 5 --nwildcard --add pileup_fastas/OL1987.fasta aln_files/all_modern_plus_HVT_clean_masked_aln_9anc.fasta 1> aln_files/all_modern_plus_HVT_clean_masked_aln_10anc.fasta 2> aln_files/log-all_modern_plus_HVT_clean_masked_aln_10anc.txt 

mafft --maxiterate 1000 --thread 5 --nwildcard --add pileup_fastas/OL2178.fasta aln_files/all_modern_plus_HVT_clean_masked_aln_10anc.fasta 1> aln_files/all_modern_plus_HVT_clean_masked_aln_11anc.fasta 2> aln_files/log-all_modern_plus_HVT_clean_masked_aln_11anc.txt 

cp aln_files/all_modern_plus_HVT_clean_masked_aln_11anc.fasta  aln_files/all_modern_plus_HVT_clean_masked_aln_anc_MLE.fasta 

# cleanup the intermediate files and the logs
rm aln_files/all_modern_plus_HVT_clean_masked_aln_1anc.fasta aln_files/all_modern_plus_HVT_clean_masked_aln_2anc.fasta aln_files/all_modern_plus_HVT_clean_masked_aln_3anc.fasta aln_files/all_modern_plus_HVT_clean_masked_aln_4anc.fasta aln_files/all_modern_plus_HVT_clean_masked_aln_5anc.fasta
rm aln_files/all_modern_plus_HVT_clean_masked_aln_6anc.fasta aln_files/all_modern_plus_HVT_clean_masked_aln_7anc.fasta aln_files/all_modern_plus_HVT_clean_masked_aln_8anc.fasta aln_files/all_modern_plus_HVT_clean_masked_aln_9anc.fasta aln_files/all_modern_plus_HVT_clean_masked_aln_10anc.fasta aln_files/all_modern_plus_HVT_clean_masked_aln_11anc.fasta
rm aln_files/log*txt


# This is the final alignment for the tree all_modern_plus_HVT_clean_masked_aln_rMLE.fasta


# 7) Build alignments for BEAST

# Laurent said that mafft introduces noise, so if you bear in mind seqs have already been aligned to the same genomes you can cat everything you need together. Plus HVT dataset for the rooted trees stays the same as the HVT outgroup isn't of the same size as the rest of the genome

cat masked_mod_files/all_modern_minus_HVT_clean_masked.fa \
	pileup_fastas/OL1099.fasta \
	pileup_fastas/OL1385.fasta \
	pileup_fastas/OL1389.fasta \
	pileup_fastas/OL1986.fasta \
	pileup_fastas/OL2272.fasta > aln_files/all_modern_minus_HVT_clean_masked_aln_anc_BT.fasta




# 8) Build trees - Here is where I am building the trees: /home/antony/Marek-Capture/masked_data/BEAST_seq_prep_dir/genome_seq_archive/modern_genomes_chunks/BEAST_alns_and_trees

# BEAST TREE 

raxmlHPC-PTHREADS -f a -T 10 -x 12345 -k -# 100 -p 12345 -m GTRGAMMA -o NC_002641.1 -s aln_files/all_modern_plus_HVT_clean_masked_aln_anc_BT.fasta -n all_modern_plus_HVT_clean_masked_aln_anc_BT -w /media/jbod/home/antony/Marek-Capture/masked_data/NEW_ANAL_2021/MLE_trees &> MLE_trees/log-make_BEAST_TREE.txt &


# GENERAL RAxML tree with 8 ancient sequences

raxmlHPC-PTHREADS -f a -T 10 -x 12345 -k -# 100 -p 12345 -m GTRGAMMA -o NC_002641.1 -s aln_files/all_modern_plus_HVT_clean_masked_aln_anc_MLE.fasta -n all_modern_plus_HVT_clean_masked_aln_anc_MLE -w /media/jbod/home/antony/Marek-Capture/masked_data/NEW_ANAL_2021/MLE_trees &> MLE_trees/log-make_RAxML_TREE.txt &



# 9) Find the overlapping regions between the open reading frames. Use the custom script made to parse the bed file created by Steven

python scripts/remove_overlaps_from_orfs.py bed_files/MDV_gene_loci.bed bed_files/MDV_gene_loci_no_overlap.bed 


# 10) Extract the respective coding sequences from the BEAST alignment 

# a. Use msaview and the bed file with the non overlapping ORFs. Dir of BED File:Users/edimopoulos/Marek-Steven/Capture_Data/BEAST_seq_prep_dir/bed-files/New_MDV_gene_loci.bed , from this directory: /Users/edimopoulos/Marek-Steven/Capture_Data/BEAST_seq_prep_dir/genome_seq_archive/modern_genomes

mkdir indiv_genes

parallel --jobs 4 "seqkit subseq -r {1}:{2} aln_files/all_modern_minus_HVT_clean_masked_aln_anc_BT.fasta > indiv_genes/{3}_aln.fasta" ::: 108 570 1206 1574 1934 2195 2548 2807 3236 3624 4736 5361 5561 6152 7052 8575 8932 9538 9922 10123 10684 11342 13435 13712 14004 14706 15202 15368 18098 18410 18587 19396 19839 20308 21258 21552 21982 22476 23336 25980 28149 28936 33883 35211 35466 37016 38018 38549 39028 40111 41220 43553 44841 45926 50239 50479 51376 53163 55787 56846 57809 59601 61664 62024 64395 66968 70812 74475 75317 77243 77738 78655 79111 88832 89399 92540 92910 94544 97065 98144 100126 101395 102877 103529 104694 105080 105707 106014 107862 110527 111918 112813 113101 114395 114632 115229 118454 119643 121194 122171 123027 123612 124099 124339 125850 127107 127407 127629 128502 128942 129375 129881 131346 132076 132223 132516 133210 133534 135040 135782 136518 136874 137394 138395 138949 139132 139730 139838 140203 140573 141032 141405 141943 142585 143163 144318 145001 151967 152214 152502 152705 152933 153143 153193 153394 154029 154254 154775 155481 156310 157062 158344 159157 159268 159854 160586 161213 162516 163700 163991 165995 166840 167110 167207 167426 167484 167767 168073 168220 168559 175721 177059 177532 :::+ 398 935 1307 1768 2137 2236 2610 3208 3385 3947 4804 5468 5800 6559 7300 8805 9129 9585 10116 10263 10836 13120 13626 13771 14390 14966 15297 15487 18355 18493 18928 19671 19844 21249 21272 21959 22158 23282 25912 27983 28901 31245 35157 35438 37004 38002 38263 38995 40071 41193 43409 44722 45800 50107 50409 51114 53016 55604 56797 57712 59560 60554 61720 64318 66776 70543 74396 75299 77239 77644 78571 79050 81087 89299 92248 92611 94322 97012 98096 99469 101235 102657 102969 104383 104921 105244 105880 107720 110288 111810 112667 113079 114390 114475 115225 118429 119494 121064 122003 122671 123338 123617 124242 124836 126185 127154 127535 128234 128609 129184 129566 130999 131435 132216 132417 132536 133407 133764 135288 136189 136541 136981 137606 138718 139104 139533 139792 139966 140406 140767 141133 141770 142233 142833 143306 144644 151804 152143 152291 152597 152881 152938 153157 153255 153525 154070 154583 155314 156020 156381 158117 158643 159189 159345 160081 161029 161932 163583 163720 164671 166111 166971 167172 167221 167431 167660 167862 168150 168396 175362 176047 177202 177960 :::+ MDV000.5 MDV001 MDV002 MDV002.6 MDV003.e.1 MDV003.2 MDV003.e.3 MDV003.4 MDV003.6 MDV003.8 MDV005.1 MDV005.2 MDV005.3 MDV005.4 MDV005.5 MDV005.6 MDV005.7 MDV005.8 MDV006.2 MDV006.3 MDV006.4 MDV006.5 MDV006.6 MDV007 MDV009 MDV008 MDV010.e.8 MDV010.e.9 MDV011 MDV011.5 MDV012 MDV012.8 MDV013 MDV014 MDV014.5 MDV015 MDV015.5 MDV016 MDV017 MDV018 MDV019 MDV020 MDV022 MDV023 MDV024 MDV025.2 MDV025.1 MDV026 MDV027.e.10 MDV028 MDV029 MDV027.e.11 MDV030 MDV031 MDV031.5 MDV032 MDV033 MDV034 MDV036 MDV035 MDV037 MDV038 MDV039.5 MDV040 MDV041 MDV042 MDV043 MDV044 MDV046 MDV045 MDV047 MDV048 MDV049 MDV049.5 MDV050 MDV050.5 MDV051 MDV052 MDV053 MDV054 MDV055 MDV056 MDV057 MDV057.1 MDV057.4 MDV057.8 MDV058 MDV059 MDV060 MDV061 MDV062 MDV064 MDV063 MDV063.5 MDV065 MDV066 MDV067 MDV068 MDV069 MDV070 MDV071 MDV071.4 MDV071.8 MDV072 MDV072.4 MDV072.6 MDV072.8 MDV073 MDV073.4 MDV074 MDV075.1 MDV075.2 MDV075.4 MDV075.5 MDV075.6 MDV075.7 MDV075.8 MDV075.9 MDV075.91 MDV075.92 MDV076 MDV076.8 MDV077.5 MDV078.1 MDV078.2 MDV078.3 MDV078.e.16 MDV078.4 MDV078.e.18 MDV078.5 MDV079 MDV080 MDV080.5 MDV081 MDV081.5 MDV082 MDV084 MDV083 MDV084.5 MDV085 MDV085.3 MDV085.6 MDV085.9 MDV086 MDV086.1 MDV086.4 MDV086.6 MDV087 MDV088 MDV089 MDV090 MDV091 MDV091.5 MDV092 MDV092.8 MDV093 MDV094 MDV095 MDV095.5 MDV096 MDV097.3 MDV097.9 MDV098 MDV098.3 MDV098.6 MDV098.9 MDV099 MDV099.5 MDV101 MDV100 MDV102 MDV102.5 MDV103


# b. Correct the directionality of the ORFs in the - strand (thus find the reverse complement of those sequences). 

# b1. Use grep to initially find the ones that are reversed and then use seqtk to reverse complement them. Run the grep command in this directory: /Users/edimopoulos/Marek-Steven/Capture_Data/BEAST_seq_prep_dir/bed-files

grep '-' bed_files/MDV_gene_loci_no_overlap.bed  | cut -f 4 | paste -s -d ' ' 

# b2. Reverse complement them with seqtk. Run this command from this directory: /Users/edimopoulos/Marek-Steven/Capture_Data/BEAST_seq_prep_dir/genome_seq_archive/modern_genomes/indiv-genes

parallel --jobs 4 "seqtk seq -r -l70 indiv_genes/{1}_aln.fasta > indiv_genes/{1}_aln_rev_comp.fasta" ::: MDV000.5 MDV002 MDV003.e.1 MDV003.e.3 MDV003.4 MDV003.8 MDV005.1 MDV005.7 MDV006.2 MDV006.3 MDV006.5 MDV006.6 MDV007 MDV009 MDV011.5 MDV016 MDV017 MDV020 MDV023 MDV024 MDV025.2 MDV025.1 MDV026 MDV028 MDV029 MDV030 MDV031 MDV031.5 MDV032 MDV034 MDV036 MDV040 MDV041 MDV042 MDV044 MDV046 MDV049 MDV050 MDV054 MDV057.8 MDV059 MDV060 MDV061 MDV062 MDV064 MDV065 MDV069 MDV071 MDV072 MDV072.8 MDV073 MDV073.4 MDV075.4 MDV075.7 MDV075.9 MDV075.91 MDV075.92 MDV076.8 MDV078.2 MDV078.4 MDV078.5 MDV080 MDV082 MDV084 MDV084.5 MDV085 MDV085.6 MDV085.9 MDV086.1 MDV086.4 MDV086.6 MDV090 MDV091 MDV092.8 MDV098 MDV098.9 MDV101 MDV102.5 MDV103

# Remove the files with the wrong directionality. Run this command from this directory: /Users/edimopoulos/Marek-Steven/Capture_Data/BEAST_seq_prep_dir/genome_seq_archive/modern_genomes/indiv-genes

parallel --jobs 4 "rm indiv_genes/{1}_aln.fasta" ::: MDV000.5 MDV002 MDV003.e.1 MDV003.e.3 MDV003.4 MDV003.8 MDV005.1 MDV005.7 MDV006.2 MDV006.3 MDV006.5 MDV006.6 MDV007 MDV009 MDV011.5 MDV016 MDV017 MDV020 MDV023 MDV024 MDV025.2 MDV025.1 MDV026 MDV028 MDV029 MDV030 MDV031 MDV031.5 MDV032 MDV034 MDV036 MDV040 MDV041 MDV042 MDV044 MDV046 MDV049 MDV050 MDV054 MDV057.8 MDV059 MDV060 MDV061 MDV062 MDV064 MDV065 MDV069 MDV071 MDV072 MDV072.8 MDV073 MDV073.4 MDV075.4 MDV075.7 MDV075.9 MDV075.91 MDV075.92 MDV076.8 MDV078.2 MDV078.4 MDV078.5 MDV080 MDV082 MDV084 MDV084.5 MDV085 MDV085.6 MDV085.9 MDV086.1 MDV086.4 MDV086.6 MDV090 MDV091 MDV092.8 MDV098 MDV098.9 MDV101 MDV102.5 MDV103


# c. Find out which orfs have at least a 10% overlap among all the sequences in the alignment. Use the custom made script percent_overlap_orfs.py , from this directory: /Users/edimopoulos/Marek-Steven/Capture_Data/BEAST_seq_prep_dir

python scripts/percent_overlap_orfs.py coding `ls indiv_genes/*_aln.fasta`


# d. Use seqtk to concatenate all the sequences together. Run this command from this folder: /Users/edimopoulos/Marek-Steven/Capture_Data/BEAST_seq_prep_dir

mkdir BEAST_partitions

seqkit concat `cat selected_orfs_coding.txt` > BEAST_partitions/MDV_CDS.fasta



# 11) Find the intragenic/non-coding regions of the genome. 


# a. merge all the expanded ORFs into overlapping regions 

# We use the initial bed file to find the total coding regions

bedtools merge -i bed_files/MDV_gene_loci.bed > bed_files/MDV_gene_merged_loci_tmp.bed

# a1. find the complements of the merged coding regions 

cat ref_genome/EF523390.1_mask.fasta.fai | cut -f 1-2 > bed_files/EF523390.1_mask_gen_file.txt

bedtools complement -i bed_files/MDV_gene_merged_loci_tmp.bed -g bed_files/EF523390.1_mask_gen_file.txt > bed_files/MDV_gene_non_coding.bed

python scripts/fix_bedtools_coord.py bed_files/MDV_gene_non_coding.bed bed_files/MDV_gene_non_coding_1idx.bed

# b. run the following command from this dir: /Users/edimopoulos/Marek-Steven/Capture_Data/BEAST_seq_prep_dir/genome_seq_archive/modern_genomes

mkdir non_coding
parallel --jobs 4 "seqkit subseq -r {1}:{2} aln_files/all_modern_minus_HVT_clean_masked_aln_anc_BT.fasta > non_coding/{3}_aln.fasta" ::: 1 399 1769 2139 2504 2612 3387 3948 5826 6560 7301 9131 9826 10266 10996 13121 13627 15006 15298 17543 18356 19672 21250 21961 22159 23283 25913 28903 31246 33785 35158 38005 38998 40072 41194 43410 44723 45801 50108 51116 53017 55605 57715 59561 61593 64321 66777 70544 75301 77646 78572 79051 89300 92613 94323 97013 98097 99470 101236 102658 102970 104384 104922 105882 107721 110289 111811 112668 114391 115228 119497 121065 122004 122672 123619 126956 127538 128611 129290 129567 131000 131658 132420 132804 133766 135289 136190 137608 138719 139535 139794 140147 140408 141772 142234 142834 143307 144645 152146 152599 152882 153159 153256 153526 154586 155315 156021 156952 158118 159191 160477 161030 162410 163584 165216 166337 166972 167173 167433 167661 168152 175365 176048 177203 177961 :::+ 107 569 1933 2194 2547 2806 3623 4735 6151 7051 8574 9537 9921 10683 11341 13434 13711 15201 15367 18097 18409 19838 21257 21981 22475 23335 25979 28935 31258 33882 35210 38017 39027 40110 41219 43552 44840 45925 50238 51375 53162 55786 57808 59600 61663 64394 66967 70811 75316 77737 78654 79110 89398 92909 94543 97064 98143 100125 101394 102876 103528 104693 105079 106013 107861 110526 111917 112812 114394 115228 119642 121193 122170 123026 124098 127106 127628 128941 129374 129880 131345 132075 132515 133209 135039 135781 136517 138394 138948 139729 139837 140202 140572 141942 142584 143162 144317 145000 152213 152704 152932 153192 153393 154028 154774 155480 156309 157061 158343 159267 160585 161212 162515 163699 165994 166839 167109 167206 167483 167766 168219 175720 177058 177531 178245 :::+ Intragenic_1 Intragenic_2 Intragenic_3 Intragenic_4 Intragenic_5 Intragenic_6 Intragenic_7 Intragenic_8 Intragenic_9 Intragenic_10 Intragenic_11 Intragenic_12 Intragenic_13 Intragenic_14 Intragenic_15 Intragenic_16 Intragenic_17 Intragenic_18 Intragenic_19 Intragenic_20 Intragenic_21 Intragenic_22 Intragenic_23 Intragenic_24 Intragenic_25 Intragenic_26 Intragenic_27 Intragenic_28 Intragenic_29 Intragenic_30 Intragenic_31 Intragenic_32 Intragenic_33 Intragenic_34 Intragenic_35 Intragenic_36 Intragenic_37 Intragenic_38 Intragenic_39 Intragenic_40 Intragenic_41 Intragenic_42 Intragenic_43 Intragenic_44 Intragenic_45 Intragenic_46 Intragenic_47 Intragenic_48 Intragenic_49 Intragenic_50 Intragenic_51 Intragenic_52 Intragenic_53 Intragenic_54 Intragenic_55 Intragenic_56 Intragenic_57 Intragenic_58 Intragenic_59 Intragenic_60 Intragenic_61 Intragenic_62 Intragenic_63 Intragenic_64 Intragenic_65 Intragenic_66 Intragenic_67 Intragenic_68 Intragenic_69 Intragenic_70 Intragenic_71 Intragenic_72 Intragenic_73 Intragenic_74 Intragenic_75 Intragenic_76 Intragenic_77 Intragenic_78 Intragenic_79 Intragenic_80 Intragenic_81 Intragenic_82 Intragenic_83 Intragenic_84 Intragenic_85 Intragenic_86 Intragenic_87 Intragenic_88 Intragenic_89 Intragenic_90 Intragenic_91 Intragenic_92 Intragenic_93 Intragenic_94 Intragenic_95 Intragenic_96 Intragenic_97 Intragenic_98 Intragenic_99 Intragenic_100 Intragenic_101 Intragenic_102 Intragenic_103 Intragenic_104 Intragenic_105 Intragenic_106 Intragenic_107 Intragenic_108 Intragenic_109 Intragenic_110 Intragenic_111 Intragenic_112 Intragenic_113 Intragenic_114 Intragenic_115 Intragenic_116 Intragenic_117 Intragenic_118 Intragenic_119 Intragenic_120 Intragenic_121 Intragenic_122 Intragenic_123 Intragenic_124 Intragenic_125

# c. Find out which intragenic regions have at least a 10% overlap among all the sequences in the alignment

python scripts/percent_overlap_orfs.py non_coding `ls non_coding/*_aln.fasta`

# d. Use seqtk to concatenate all the sequences together. Run this command from this folder: /Users/edimopoulos/Marek-Steven/Capture_Data/BEAST_seq_prep_dir

seqkit concat `cat selected_orfs_non_coding.txt` > BEAST_partitions/MDV_non_coding.fasta




# 12) Partition MDV alignments for BEAST 

python ../../../codonextract.py ./MDV-CDS.fasta > codon1pos.fasta 
python ../../../codonextract.py ./MDV-CDS.fasta > codon2pos.fasta 
python ../../../codonextract.py ./MDV-CDS.fasta > codon3pos.fasta























