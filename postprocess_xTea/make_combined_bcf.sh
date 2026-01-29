#!/bin/bash

results_dir=$1
dataset_prefix=$2
bcf_dir=${results_dir}/../bcfs
remove_bcfs=$3

if [ ! -d $bcf_dir ]; then
	mkdir $bcf_dir
fi

for i in ${results_dir}/*/L1/*filtered_LINE1.vcf; do
    name=${i##*/}
    bcftools view -Ob -o ${bcf_dir}/${name}.bcf ${i}; bcftools index ${bcf_dir}/${name}.bcf;
done


bcftools merge -m none -0 ${bcf_dir}/*LINE1*bcf > ${bcf_dir}/${dataset_prefix}.vcf
bcftools +fill-tags ${bcf_dir}/${dataset_prefix}.vcf -- -t AN,AC,AF > ${results_dir}/../L1_${dataset_prefix}.vcf



for i in ${results_dir}/*/Alu/*filtered_ALU.vcf; do
    name=${i##*/}
    bcftools view -Ob -o ${bcf_dir}/${name}.bcf ${i}; bcftools index ${bcf_dir}/${name}.bcf;
done


bcftools merge -m none -0 ${bcf_dir}/*ALU*bcf > ${bcf_dir}/${dataset_prefix}.vcf
bcftools +fill-tags ${bcf_dir}/${dataset_prefix}.vcf -- -t AN,AC,AF > ${results_dir}/../ALU_${dataset_prefix}.vcf



for i in ${results_dir}/*/SVA/*filtered_SVA.vcf; do
    name=${i##*/}
    bcftools view -Ob -o ${bcf_dir}/${name}.bcf ${i}; bcftools index ${bcf_dir}/${name}.bcf;
done


bcftools merge -m none -0 ${bcf_dir}/*SVA*bcf > ${bcf_dir}/${dataset_prefix}.vcf
bcftools +fill-tags ${bcf_dir}/${dataset_prefix}.vcf -- -t AN,AC,AF > ${results_dir}/../SVA_${dataset_prefix}.vcf

if [[ "$remove_bcfs" == "remove" ]]; then
	rm -rf ${bcf_dir}
fi
