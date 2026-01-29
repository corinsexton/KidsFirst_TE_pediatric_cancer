#!/bin/bash

# example xtea run script


#SBATCH --job-name=<job_name>
#SBATCH -A <account>
#SBATCH --partition <partition>
#SBATCH --mem 5G
#SBATCH -c 1
#SBATCH -t 02:00:00
#SBATCH -o slurm_xtea_-%x.%j.out


<path to xtea executable>/xtea -i ids.txt -b id_cram.txt \
       -x null -p results/ \
       -o output_submit_script.sh -l /<path to xtea>/xTea/rep_lib_annotation/ \
       -r Homo_sapiens_assembly38.fa \
       -g gencode.v45.basic.annotation.gff3 \
       --xtea <path to xtea>/xTea/xtea/ -f 5907 -y 7  \
       --slurm -t 0-04:00 -n 8 -m 24
