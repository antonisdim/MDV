#Require:
#msa_view (PHAST)
#seqtk
#union (EMBOSS)
#prank
#codeml


#bedfile of REF genome processed to separate exons which can be concatenated to form CDS.
bedfile=/media/jbod/home/steven/MDV/final_PAML_new_cutoffs/EF523390.1.bed
#genomic aln
aln=/media/jbod/home/steven/MDV/final_PAML_new_cutoffs/mdv_mod_anc_no_HVT_aln_BEAST.fasta
#output file for the alns
alnout=/media/jbod/home/steven/MDV/final_PAML_new_cutoffs/out_alns
#tree file without HVT, annotated for PAML
tree=/media/jbod/home/steven/MDV/final_PAML_new_cutoffs/RAxML_bipartitions.mdv_no_HVT_MLE




## SLICE ALN ##
#get position of the REF genome in the alignment
posn=`cat $aln | grep \> | awk '1;/EF523390.1/{exit}' | wc -l`
#use msa_view to slice the aln
while read i; do locus=`echo $i | awk '{print $4}'`; start=`echo $i | awk '{print $2}'`; end=`echo $i | awk '{print $3}'`; msa_view $aln --start $start --end $end --refidx $posn | sed 's/> />/' > $alnout/"$locus"_aln.fa; done < EF523390.1.bed



## CONCAT MULTI-EXONS ##
#get a list of the multi-exonic genes
cat $bedfile | grep \.e\. | awk '{print $4}' | awk -F ".e." '{print $1}' | sort | uniq > to_concat.txt
#get multi-exon loci into a single file for each locus
for i in `cat to_concat.txt`; do for j in `grep "$i".e. $bedfile | awk '{print $4}'`; do cat $alnout/"$j"_aln.fa >> $alnout/"$i"_aln_to_concat.fa; done; done
#for each locus, for each ID, find the sequences, concatenate, and then output into the concatenated file.
for i in `cat to_concat.txt`; do for j in `cat $aln | grep ">"`; do awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' $alnout/"$i"_aln_to_concat.fa | grep $j | tr "\t" "\n" | union -filter >> $alnout/"$i"_aln.fa; done; done





## RC ORFs ON REVERSE STRAND ##
#get list of ORFs on reverse strand, take into account that individual exon files are not used any more
cat $bedfile | awk '{print $6"\t"$4}' | awk -F ".e." '{print $1}' | sort -k 2 | uniq | grep -v + | awk '{print $2}' > to_RC.txt
#use seqtk to RC these
for i in `cat to_RC.txt`; do seqtk seq -r $alnout/"$i"_aln.fa > $alnout/"$i"_aln_rc.fa; done


## PAML ##
bedfile=/drives/ssd0/steven/MDV/final_final_new_cutoffs/EF523390.1.bed
tree=/drives/ssd0/steven/MDV/final_final_new_cutoffs/RAxML_bipartitions.mdv_no_HVT_MLE
#final list of loci to include
cat $bedfile | awk '{print $6"\t"$4}' | awk -F ".e." '{print $1}' | sort | uniq > all_loci.txt
#directory structure for results
for i in `cat all_loci.txt | awk '{print $2}'`; do mkdir -p results/"$i"/{alt,null}; done
#modify codeml.ctl template file for each locus
while read i; do locus=`echo $i | awk '{print $2}'`; strand=`echo $i | awk '{print $1}'`; if [[ $strand == "-" ]]; then aln=`readlink -f out_alns/"$locus"_aln_rc.fa`; else aln=`readlink -f out_alns/"$locus"_aln.fa`; fi; for j in alt null; do cp codeml.ctl results/"$locus"/"$j"; if [[ $j == "alt" ]]; then sed -i -e "s;SEQFILE;"$aln";g" -e "s/LOCUSNAME/"$locus"/g" -e  "s/OMEGATOCHANGE/0/g" results/"$locus"/"$j"/codeml.ctl; else sed -i -e "s;SEQFILE;"$aln";g" -e "s/LOCUSNAME/"$locus"/g" -e  "s/OMEGATOCHANGE/1/g" results/"$locus"/"$j"/codeml.ctl; fi; done; done < all_loci.txt
#generate list of jobs for PAML to add into next command
#for i in `cat all_loci.txt | awk '{print $2}'`; do for j in alt null; do printf $i"/"$j" "; done; done
#run PAML
parallel --jobs 60 "cd results/{1} && codeml" ::: MDV000.5/alt MDV000.5/null MDV001/alt MDV001/null MDV002/alt MDV002/null MDV002.3/alt MDV002.3/null MDV002.6/alt MDV002.6/null MDV003/alt MDV003/null MDV003.2/alt MDV003.2/null MDV003.4/alt MDV003.4/null MDV003.6/alt MDV003.6/null MDV003.8/alt MDV003.8/null MDV004/alt MDV004/null MDV005/alt MDV005/null MDV005.1/alt MDV005.1/null MDV005.2/alt MDV005.2/null MDV005.3/alt MDV005.3/null MDV005.4/alt MDV005.4/null MDV005.5/alt MDV005.5/null MDV005.6/alt MDV005.6/null MDV005.7/alt MDV005.7/null MDV005.8/alt MDV005.8/null MDV006/alt MDV006/null MDV006.1/alt MDV006.1/null MDV006.2/alt MDV006.2/null MDV006.3/alt MDV006.3/null MDV006.4/alt MDV006.4/null MDV006.5/alt MDV006.5/null MDV006.6/alt MDV006.6/null MDV007/alt MDV007/null MDV008/alt MDV008/null MDV008.4/alt MDV008.4/null MDV008.8/alt MDV008.8/null MDV009/alt MDV009/null MDV009.5/alt MDV009.5/null MDV010/alt MDV010/null MDV011/alt MDV011/null MDV011.5/alt MDV011.5/null MDV012/alt MDV012/null MDV012.4/alt MDV012.4/null MDV012.8/alt MDV012.8/null MDV013/alt MDV013/null MDV013.5/alt MDV013.5/null MDV014/alt MDV014/null MDV014.5/alt MDV014.5/null MDV015/alt MDV015/null MDV015.5/alt MDV015.5/null MDV016/alt MDV016/null MDV017/alt MDV017/null MDV018/alt MDV018/null MDV019/alt MDV019/null MDV020/alt MDV020/null MDV020.5/alt MDV020.5/null MDV021/alt MDV021/null MDV022/alt MDV022/null MDV023/alt MDV023/null MDV024/alt MDV024/null MDV025.1/alt MDV025.1/null MDV025.2/alt MDV025.2/null MDV026/alt MDV026/null MDV027/alt MDV027/null MDV028/alt MDV028/null MDV029/alt MDV029/null MDV030/alt MDV030/null MDV031/alt MDV031/null MDV031.5/alt MDV031.5/null MDV032/alt MDV032/null MDV033/alt MDV033/null MDV034/alt MDV034/null MDV035/alt MDV035/null MDV036/alt MDV036/null MDV037/alt MDV037/null MDV038/alt MDV038/null MDV038.5/alt MDV038.5/null MDV039/alt MDV039/null MDV039.5/alt MDV039.5/null MDV040/alt MDV040/null MDV041/alt MDV041/null MDV042/alt MDV042/null MDV043/alt MDV043/null MDV044/alt MDV044/null MDV045/alt MDV045/null MDV046/alt MDV046/null MDV047/alt MDV047/null MDV048/alt MDV048/null MDV049/alt MDV049/null MDV049.1/alt MDV049.1/null MDV049.5/alt MDV049.5/null MDV050/alt MDV050/null MDV050.5/alt MDV050.5/null MDV051/alt MDV051/null MDV052/alt MDV052/null MDV053/alt MDV053/null MDV054/alt MDV054/null MDV055/alt MDV055/null MDV056/alt MDV056/null MDV057/alt MDV057/null MDV057.1/alt MDV057.1/null MDV057.4/alt MDV057.4/null MDV057.8/alt MDV057.8/null MDV058/alt MDV058/null MDV059/alt MDV059/null MDV060/alt MDV060/null MDV061/alt MDV061/null MDV062/alt MDV062/null MDV063/alt MDV063/null MDV063.5/alt MDV063.5/null MDV064/alt MDV064/null MDV065/alt MDV065/null MDV066/alt MDV066/null MDV067/alt MDV067/null MDV068/alt MDV068/null MDV069/alt MDV069/null MDV070/alt MDV070/null MDV071/alt MDV071/null MDV071.4/alt MDV071.4/null MDV071.8/alt MDV071.8/null MDV072/alt MDV072/null MDV072.2/alt MDV072.2/null MDV072.4/alt MDV072.4/null MDV072.6/alt MDV072.6/null MDV072.8/alt MDV072.8/null MDV073/alt MDV073/null MDV073.4/alt MDV073.4/null MDV074/alt MDV074/null MDV075/alt MDV075/null MDV075.1/alt MDV075.1/null MDV075.2/alt MDV075.2/null MDV075.3/alt MDV075.3/null MDV075.4/alt MDV075.4/null MDV075.5/alt MDV075.5/null MDV075.6/alt MDV075.6/null MDV075.7/alt MDV075.7/null MDV075.8/alt MDV075.8/null MDV075.9/alt MDV075.9/null MDV075.91/alt MDV075.91/null MDV075.92/alt MDV075.92/null MDV076/alt MDV076/null MDV076.4/alt MDV076.4/null MDV076.8/alt MDV076.8/null MDV077/alt MDV077/null MDV077.5/alt MDV077.5/null MDV078/alt MDV078/null MDV078.1/alt MDV078.1/null MDV078.2/alt MDV078.2/null MDV078.3/alt MDV078.3/null MDV078.4/alt MDV078.4/null MDV078.5/alt MDV078.5/null MDV078.6/alt MDV078.6/null MDV079/alt MDV079/null MDV080/alt MDV080/null MDV080.5/alt MDV080.5/null MDV081/alt MDV081/null MDV081.5/alt MDV081.5/null MDV082/alt MDV082/null MDV083/alt MDV083/null MDV084/alt MDV084/null MDV084.5/alt MDV084.5/null MDV085/alt MDV085/null MDV085.3/alt MDV085.3/null MDV085.6/alt MDV085.6/null MDV085.9/alt MDV085.9/null MDV086/alt MDV086/null MDV086.1/alt MDV086.1/null MDV086.2/alt MDV086.2/null MDV086.4/alt MDV086.4/null MDV086.6/alt MDV086.6/null MDV087/alt MDV087/null MDV088/alt MDV088/null MDV089/alt MDV089/null MDV089.5/alt MDV089.5/null MDV090/alt MDV090/null MDV091/alt MDV091/null MDV091.5/alt MDV091.5/null MDV092/alt MDV092/null MDV092.4/alt MDV092.4/null MDV092.8/alt MDV092.8/null MDV093/alt MDV093/null MDV094/alt MDV094/null MDV094.5/alt MDV094.5/null MDV095/alt MDV095/null MDV095.5/alt MDV095.5/null MDV096/alt MDV096/null MDV096.5/alt MDV096.5/null MDV097.3/alt MDV097.3/null MDV097.6/alt MDV097.6/null MDV097.9/alt MDV097.9/null MDV098/alt MDV098/null MDV098.3/alt MDV098.3/null MDV098.6/alt MDV098.6/null MDV098.9/alt MDV098.9/null MDV099/alt MDV099/null MDV099.5/alt MDV099.5/null MDV100/alt MDV100/null MDV101/alt MDV101/null MDV102/alt MDV102/null MDV102.5/alt MDV102.5/null MDV103/alt MDV103/null

## PROCESS RESULTS ##
#outputs file with 8 columns: locus ID, null lnL, alt lnL, 2* lnL difference, p-value, no. PS sites (any significance level), no. PS sites (0.9 < x < 0.99), no. PS sites (>0.99)
while read i;  do  locus=`echo $i | awk '{print $2}'`;  if  [ ! -s results/$locus/alt/$locus ];  then  printf $locus"\t"NaN"\t"NaN"\t"NaN"\t"NaN"\t"NaN"\t"NaN"\t"NaN; echo "" ; fi;  if grep -q lnL results/$locus/null/$locus;  then  null=`cat results/$locus/null/$locus | grep lnL | awk -F ":" '{print $4}' | awk -F "+" '{print $1}' | xargs`; else  null=NaN;  fi;  if grep -q lnL results/$locus/alt/$locus;  then  alt=`cat results/$locus/alt/$locus | grep lnL | awk -F ":" '{print $4}' | awk -F "+" '{print $1}' | xargs` ; else  alt=NaN;  fi;  if [ $null == "NaN" ];  then  printf $locus"\t"NaN"\t"NaN"\t"NaN"\t"NaN"\t"NaN"\t"NaN"\t"NaN; echo "" ; else  diff=`echo "$alt - $null" | bc`;  lnl=`echo "$diff * 2" | bc`;  sites=`awk '/Nielsen/,/The/' results/$locus/alt/$locus | tail -n +3 | head -n -3 | wc -l`; sigsites=`awk '/Nielsen/,/The/' results/$locus/alt/$locus | tail -n +3 | head -n -3 | grep "\*" | grep -v "\*\*" | wc -l`; hisigsites=`awk '/Nielsen/,/The/' results/$locus/alt/$locus | tail -n +3 | head -n -3 | grep "\*\*" | wc -l`; if (( $(echo "$lnl > 0" | bc -l) ));  then p=`chi2 1 $lnl | awk -F "=" '{print $4}' | xargs`;  else p=1; fi; printf $locus"\t"$null"\t"$alt"\t"$lnl"\t"$p"\t"$sites"\t"$sigsites"\t"$hisigsites; echo "";  fi;  done < all_loci.txt > all_results.tab

#or this which outputs additional information about the significance levels of sites.
#while read i; do locus=`echo $i | awk '{print $1}'`; line=`grep -P $locus'\t' EF523390.1.bed`; start=`echo $line | awk '{print $2}'`; end=`echo $line | awk '{print $3}'`; awk '/Nielsen/,/The/' results/$locus/alt/$locus | tail -n +3 | head -n -3 | grep "\*" > tmp; while read j; do protein_site=`echo $j | awk -F " " '{print $1}' | xargs`; ntsite=`echo "$protein_site * 3" | bc`; genome_site=`echo "$ntsite + $start" | bc`; star=`echo $j | awk -F. '{print $2}' | sed 's/[[:digit:]]//g'`; if [[ $star == "*" ]]; then colour="color=yellow"; elif [[ $star == "**" ]]; then colour="color=orange"; fi; printf $locus"\t"$start"\t"$end"\t"$protein_site"\t"$ntsite"\t"$genome_site"\t"$colour; echo ""; done < tmp; done < all_results_with_PSS.tab
