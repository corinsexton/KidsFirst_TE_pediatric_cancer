#!/bin/bash

# find /home/cos689/mei-pediatric-cancer/xtea -maxdepth 5 -name *filtered_*.vcf > KF_vcf_files_xtea.txt
#./get_EUR_sample.py > EUR_KF_vcfs.txt

while read -r f; do

	dataset=${f#/home/cos689/mei-pediatric-cancer/xtea/}
	dataset=${dataset%%/*}

	sample=${f%/*/*vcf}
	sample=${sample##*/}

	te=${f%/*vcf}
	te=${te##*/}

	#	1KG_ALU_gt0.9_EUR.bed 
	if [[ "$te" == "Alu" ]]; then
		count=$(bedtools intersect -b alu_gt.95_gnom_kg.bed \
			-a $f -wa -u  | wc -l)
		echo $dataset $sample $te $count
		
	#elif [[ "$te" == "L1" ]]; then
	#	count=$(bedtools intersect -b /home/cos689/data1/corinne/ref/1kgenomes/1KG_L1_gt0.5.bed \
	#		-a $f -wa -u  | wc -l)
	#	echo $dataset $sample $te $count

	fi


#done < EUR_KF_vcfs.txt
done < ALL_KF_vcfs.txt
