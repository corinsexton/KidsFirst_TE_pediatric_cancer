#!/usr/bin/env python3

import csv
import sys
import glob

infile=sys.argv[1]
with open(infile, 'r', newline='', encoding='utf-8') as f:
    reader = csv.reader(f)
    header = next(reader)
    
    # Identify the data columns (all columns past the second)
    col_names = header[2:]
    
    # Dictionary to hold sums of columns by Superpopulation
    sums_by_superpop = {}
    
    # Read each row and accumulate sums
    for row in reader:
        if not row:
            continue
        superpop = row[1]   # 2nd column: Superpopulation
        values   = [float(x) for x in row[2:]]  # columns from 3rd onward

        if superpop not in sums_by_superpop:
            sums_by_superpop[superpop] = [0.0] * len(values)
        
        for i, val in enumerate(values):
            sums_by_superpop[superpop][i] += val
    
    # For each data column, find which superpop is highest
    for i, col_name in enumerate(col_names):
        best_superpop = None
        best_val = float('-inf')
        for sp, sums in sums_by_superpop.items():
            if sums[i] > best_val:
                best_val = sums[i]
                best_superpop = sp
        print(col_name, best_superpop)



