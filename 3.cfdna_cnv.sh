#!/bin/sh

for i in sample
do
        cd path/1.CpG
        cd path/2.m_a
	cat /dev/null >$i.m_a_count.o
	for chr in "1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "X" "Y"
	do
		tmpchr="chr"$chr
        	awk -v OFS="\t" -v chr=$tmpchr '{if($1==chr){print $1,$2,$3}}' ../ref_tissue/cnv.1000k.chrY.bed >$i.cnv.tmp.bed
      		gunzip -c ../1.CpG/CpG*$i.*  |awk -v OFS="\t" -v chr=$tmpchr '{if($3==chr){print $3,$4,$4,$1,$2,$5}}'|awk -v OFS="\t" '{if(NF==6){print $0}}'|sort -k2,2n  >$i.cpg.tmp.sorted.txt
		#echo "tmp"
        	python ../scripts/m_a.py -i $i.cpg.tmp.sorted.txt --region $i.cnv.tmp.bed --reference ../ref_tissue/"Ref"$tmpchr"_len1000k.txt" -o $i.m_a_count.tmp.o
		cat $i.m_a_count.o $i.m_a_count.tmp.o >>$i.m_a_count.tmp2.o
		mv $i.m_a_count.tmp2.o $i.m_a_count.o
	done
	
        rm $i.cnv.tmp.bed $i.cpg.tmp.sorted.txt $i.m_a_count.tmp.o
done

