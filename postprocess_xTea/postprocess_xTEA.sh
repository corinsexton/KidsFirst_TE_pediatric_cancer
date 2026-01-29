#!/bin/bash

results_dir=$1
ids=$2
dataset_prefix=$3

gencode_gff=<path>/hg38_gencode_ann.gff3

kg_sva=<path>/1KG_SVA.bed.gz
kg_alu=<path>/1KG_L1.bed.gz
kg_l1=<path>/1KG_ALU.bed.gz

#### coded inside to do SVA, L1, and ALU
if [ ! -f ${results_dir}/../SVA_${dataset_prefix}.vcf ]; then
	./make_combined_bcf.sh ${results_dir} ${dataset_prefix} r
fi

./merge_within_35bp_vcf.py ${results_dir}/../SVA_${dataset_prefix}.vcf > ${results_dir}/../SVA_${dataset_prefix}_merged.vcf
./merge_within_35bp_vcf.py ${results_dir}/../ALU_${dataset_prefix}.vcf > ${results_dir}/../ALU_${dataset_prefix}_merged.vcf
./merge_within_35bp_vcf.py ${results_dir}/../L1_${dataset_prefix}.vcf > ${results_dir}/../L1_${dataset_prefix}_merged.vcf


bcftools annotate -c CHROM,FROM,TO,ID \
	-a ~/data1/corinne/ref/1kgenomes/1KG_ALU.bed.gz -o ${results_dir}/../ALU_${dataset_prefix}_1kg.vcf ${results_dir}/../ALU_${dataset_prefix}_merged.vcf
bcftools annotate -c CHROM,FROM,TO,ID \
	-a ~/data1/corinne/ref/1kgenomes/1KG_L1.bed.gz -o ${results_dir}/../L1_${dataset_prefix}_1kg.vcf ${results_dir}/../L1_${dataset_prefix}_merged.vcf
bcftools annotate -c CHROM,FROM,TO,ID \
	-a ~/data1/corinne/ref/1kgenomes/1KG_SVA.bed.gz -o ${results_dir}/../SVA_${dataset_prefix}_1kg.vcf ${results_dir}/../SVA_${dataset_prefix}_merged.vcf

./annotate_xtea.py ${results_dir}/../ALU_${dataset_prefix}_1kg.vcf ${results_dir}/../annotated_${dataset_prefix}_ALU.tsv
./annotate_xtea.py ${results_dir}/../L1_${dataset_prefix}_1kg.vcf ${results_dir}/../annotated_${dataset_prefix}_L1.tsv
./annotate_xtea.py ${results_dir}/../SVA_${dataset_prefix}_1kg.vcf ${results_dir}/../annotated_${dataset_prefix}_SVA.tsv

bgzip -f ${results_dir}/../ALU_${dataset_prefix}_1kg.vcf
bgzip -f ${results_dir}/../L1_${dataset_prefix}_1kg.vcf
bgzip -f ${results_dir}/../SVA_${dataset_prefix}_1kg.vcf

bcftools index -f ${results_dir}/../ALU_${dataset_prefix}_1kg.vcf.gz
bcftools index -f ${results_dir}/../L1_${dataset_prefix}_1kg.vcf.gz
bcftools index -f ${results_dir}/../SVA_${dataset_prefix}_1kg.vcf.gz


rm ${results_dir}/../SVA_${dataset_prefix}.vcf ${results_dir}/../ALU_${dataset_prefix}.vcf ${results_dir}/../L1_${dataset_prefix}.vcf
rm ${results_dir}/../ALU_${dataset_prefix}_merged.vcf ${results_dir}/../L1_${dataset_prefix}_merged.vcf ${results_dir}/../SVA_${dataset_prefix}_merged.vcf

