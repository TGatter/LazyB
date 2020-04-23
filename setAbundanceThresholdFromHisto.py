#!/usr/bin/env python
## Script reads the histofile of k-mer counting
## upper abundance threshold is set

import sys
import numpy as np

histo = sys.argv[1]
total_non_unique_k_mers = int(sys.argv[2])

# calculate quartiles 1 and 3 
q1_th_value = round((total_non_unique_k_mers + 1) * 0.25)
q1 = 0
q3_th_value = round((total_non_unique_k_mers + 1) * 0.75)
q3 = 0

current_value = 0
with open(histo, 'r') as file:
    for line in file:
        split_line = line.rstrip().split()
        abundance = int(split_line[0])
        frequency = int(split_line[1])
        if(abundance > 1):    # exclude k-mers with frequency 1 == non_unique k-mers
            current_value += frequency
            if(q1 == 0 and current_value >= q1_th_value): 
                q1 = abundance
            elif(q3 == 0 and current_value >= q3_th_value):
                q3 = abundance
                break

iqr = q3 - q1
upper_outlier = q3 + 2 * iqr
print(str(upper_outlier))

