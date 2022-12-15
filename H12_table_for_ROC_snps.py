#script takes SF results file and slim fixed file and creates table for plotting ROC curve

#python3  H12_table_for_ROC.py \
#-hFile 'garud_H12/results/stationary/rr_fixed_mu_fixed/benProp_0.0005_2Nes_100/1Mb_rep1_2kb_win.txt' \
#-fixedFile "../results/stationary/rr_fixed_mu_fixed/benProp_0.0005_2Nes_100/1Mb_rep1.fixed" \
#-outFile  "../ROC/tables/stationary/rr_fixed_mu_fixed/benProp_0.0005_2Nes_100/1Mb_rep1_200_fuzz.txt"

from __future__ import print_function
import sys
import numpy as np
import pandas as pd
import math
import argparse


#parsing user given constants
parser = argparse.ArgumentParser(description='Information about length of region and sample size')
parser.add_argument('-hFile', dest = 'hFile', action='store', nargs = 1, type = str, help = 'path to H12 results file')
parser.add_argument('-fixedFile', dest = 'fixedFile', action='store', nargs = 1, type = str, help = 'path to file containing fixed mutation data (.fixed format)')
parser.add_argument('-outFile', dest = 'outFile', action='store', nargs = 1, type = str, help = 'path to output file')

args = parser.parse_args()
h_file = args.hFile[0]
f_file = args.fixedFile[0]
out_file = args.outFile[0]


#Read in H12 results
res = pd.read_csv(h_file, sep='\t', header=0)


#Read in .fixed file
fixed = pd.read_csv(f_file, skiprows=2, sep=' ',
                   names=['tempID', 'permID', 'mutType', 'location', 's', 'h', 'initial_subpop', 'origin_gen',
                         'fix_gen'])
#Keep only beneficial mutations
sweeps = fixed[fixed.mutType=='m0']
#Keep only fixations that occured up to 0.5Ne generations ago
sweeps = sweeps[sweeps.fix_gen >= (17*7000)-(0.5*7000)]
#Reset index for sweeps df
sweeps = sweeps.reset_index(drop=True)

binary = [0 for x in range(0, len(res))]

for i in range(0, len(res)):
    for j in sweeps.location:
        if((j >= res.win_start[i]) & (j <= res.win_end[i])):
            binary[i] = 1
            
res['sweep_detected'] = binary

df = res[['sweep_detected', 'H12']]

#Output df to file
df.to_csv(out_file, sep='\t', index=False, header=True)