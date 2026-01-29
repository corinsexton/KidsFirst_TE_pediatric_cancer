#!/usr/bin/env python

import sys

if len(sys.argv) < 2:
	print('USAGE: ./merge_within_35bp_vcf.py input.vcf > output.vcf')
	sys.exit()

infile_name = sys.argv[1]

def write_out_line(line,genos):
	print('\t'.join(line) + '\t' + '\t'.join(genos))



prev_chr = 'chr1'
prev_start = 0
prev_line = ''

with open(infile_name) as infile:

	for line in infile:
		if line[0] == '#': 
			print(line.strip())
			sample_size = len(line.split()) - 9
			prev_genos = ['NA'] * sample_size
		else:
			for line in infile:
				ll = line.strip().split()
		
				curr_chr = ll[0]
				curr_start = int(ll[1])
				curr_line= ll[:9]
				curr_genos = ll[9:]
		
				if prev_chr == curr_chr and abs(prev_start - curr_start) <= 35:
						for i_gt in range(len(curr_genos)):
							if prev_genos[i_gt] == '1/1' or curr_genos[i_gt] == '1/1':
								curr_genos[i_gt] = '1/1'
							elif prev_genos[i_gt] == '0/1' or curr_genos[i_gt] == '0/1':
								curr_genos[i_gt] = '0/1'
							else:
								curr_genos[i_gt] = '0/0'
				else:
					if prev_start != 0:
						write_out_line(prev_line, prev_genos)
		
				prev_chr = curr_chr
				prev_start = curr_start
				prev_line = curr_line
				prev_genos = curr_genos
						
								
			write_out_line(prev_line, prev_genos)
								
		
							
						
						
				
		
