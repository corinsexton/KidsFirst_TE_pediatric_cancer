#!/usr/bin/env python

import pyranges as pr, pandas as pd
import argparse
import gzip


parser = argparse.ArgumentParser()
#parser.add_argument("dataset",help = 'dataset (LFS or PCAWG)')
parser.add_argument("gff_file",help = 'gff3 annotation file')
parser.add_argument("vcf_file",help = 'vcf annotation file')
parser.add_argument("outfile",help = 'output_tsv file')
args = parser.parse_args()

ann_all = pr.read_gff3(args.gff_file)

if args.vcf_file.endswith('.gz'):
	vcf_file = gzip.open(args.vcf_file, 'rt') 
else:
	vcf_file = open(args.vcf_file, 'r')

line = vcf_file.readline()
while line[0:2] == '##':
	line = vcf_file.readline()
header = line.strip().split()


vcf = pd.read_csv(args.vcf_file, sep="\t", comment='#',names = header)
vcf = vcf.rename({"#CHROM": "Chromosome"}, axis=1)
vcf = vcf.rename({"POS": "Start"}, axis=1)
#vcf['End'] = vcf.INFO.str.split(';').str[1].replace('SVLEN=','',regex = True).astype(int) + vcf.Start
vcf['End'] = vcf.Start + 1

vcf = pr.PyRanges(vcf)


vcf_file.close()

# col names:
#Index(['Chromosome', 'Source', 'Feature', 'Start', 'End', 'Score', 'Strand',
#       'Frame', 'ID', 'Parent', 'gene_id', 'transcript_id', 'gene_type',
#       'gene_name', 'transcript_type', 'transcript_name', 'exon_number',
#       'exon_id', 'level', 'transcript_support_level', 'tag',
#       'havana_transcript', 'hgnc_id', 'ont', 'havana_gene', 'protein_id',
#       'ccdsid', 'Start_b', 'ID_b', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO',
#       'FORMAT', '2349', '2760A', '3425', '3432', '3506', '3687_2', '4028',
#       '4269', '4329', '4535', '4854', '4942', '5535', '5536', '5537', '5539',
#       '5540', '5541', '5542', '5543', '5544', '5567', 'End_b'],
#      dtype='object')

intersects = ann_all.join(vcf,how = 'right',slack = 35)
intersects = intersects.sort(by = "Start_b")

out = intersects.df

index_cols = [0,28,-1,30,31,29,2,14,11,20,21] + list(range(36,len(out.columns)-1))

out = out.iloc[:,index_cols]

out['Feature'] = out["Feature"].replace('-1','intergenic')
out = out.rename({"Start_b": "Start"}, axis=1)
out = out.rename({"End_b": "End"}, axis=1)
out = out.rename({"ID_b": "ID"}, axis=1)

out= out.drop_duplicates()


out.to_csv(args.outfile,sep = '\t',index = False)






