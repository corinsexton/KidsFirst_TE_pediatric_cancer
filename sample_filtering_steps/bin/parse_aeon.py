#!/usr/bin/env python


#Population,Superpopulation,PROCA888_PROCA888
#ACB,AFR,0.0
#ASW,AFR,0.0
#BEB,SAS,0.0
#CDX,EAS,0.0
#CEU,EUR,0.24

import sys
from collections import defaultdict


infile = open(sys.argv[1])
infile.readline()

pop_dict = defaultdict(int)

for line in infile:
	ll = line.split(',')
	pop_dict[ll[1]] += float(ll[2])


prev_value = 0
for i in pop_dict:
	value = pop_dict[i]
	if value > prev_value:
		max_pop = i
		prev_value = value

print(max_pop)
