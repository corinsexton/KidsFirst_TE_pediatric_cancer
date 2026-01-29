#!/usr/bin/env python

import gzip
import sys
from collections import defaultdict

infile = gzip.open(sys.argv[1],'rb')

header = infile.readline().strip().split()
samples = header[3:]

small_bins_dict = defaultdict(int)
range_list = range(len(samples))

print("Iterating through bins...")

total = 0
for line in infile:
	nums = line.strip().split()[3:]
	total +=1

	for i in range_list:
		num = float(nums[i])
		if num < 0.85 or num > 1.15:
			small_bins_dict[samples[i]] += 1


total = float(total)

outfile = open(sys.argv[2],'w')
outfile.write("sample\todd_bins_perc\n")

final_samples = samples
		
print("Finished!")
print("Writing out results...")
for s in final_samples:
	proport = small_bins_dict[s]/total
	outfile.write(f"{s.decode('utf-8')}\t{proport}\n")

infile.close()
outfile.close()
print("Finished!")
		


