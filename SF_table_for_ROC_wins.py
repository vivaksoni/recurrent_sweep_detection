#script takes SF results file and slim fixed file and creates table for plotting ROC curve

#python3  SF_table_for_ROC.py \
#-sfFile 'SF_inputFiles/results/stationary/rr_fixed_mu_fixed/benProp_0.0005_2Nes_100/1Mb_rep1_2kb_win.out' \
#-fixedFile "../results/stationary/rr_fixed_mu_fixed/benProp_0.0005_2Nes_100/1Mb_rep1.fixed" \
#-outFile  "../ROC/tables/stationary/rr_fixed_mu_fixed/benProp_0.0005_2Nes_100/1Mb_rep1_200_fuzz.txt" \
#-winSize 500 

from __future__ import print_function
import sys
import numpy as np
import pandas as pd
import math
import argparse


#parsing user given constants
parser = argparse.ArgumentParser(description='Information about length of region and sample size')
parser.add_argument('-sfFile', dest = 'sfFile', action='store', nargs = 1, type = str, help = 'path to sweepfinder results file')
parser.add_argument('-fixedFile', dest = 'fixedFile', action='store', nargs = 1, type = str, help = 'path to file containing fixed mutation data (.fixed format)')
parser.add_argument('-outFile', dest = 'outFile', action='store', nargs = 1, type = str, help = 'path to output file')
parser.add_argument('-regionLen', dest = 'regionLen', action='store', nargs = 1, type = int, help = 'size of region')
parser.add_argument('-winSize', dest = 'winSize', action='store', nargs = 1, type = int, help = 'size of window')

args = parser.parse_args()
regionLen = args.regionLen[0]
winSize = args.winSize[0]
sf_file = args.sfFile[0]
f_file = args.fixedFile[0]
out_file = args.outFile[0]


#Read in SF results
res = pd.read_csv(sf_file, sep='\t', header=0)

#Read in .fixed file
fixed = pd.read_csv(f_file, 
                    skiprows=2, 
                    sep=' ',
                   names=['tempID', 'permID', 'mutType', 'location', 's', 'h', 'initial_subpop', 'origin_gen',
                         'fix_gen'])
#Keep only beneficial mutations
sweeps = fixed[fixed.mutType=='m0']
#Keep only fixations that occured up to 4Ne generations ago
sweeps = sweeps[sweeps.fix_gen >= (17*7000)-(0.5*7000)]
#Reset index for sweeps df
sweeps = sweeps.reset_index(drop=True)

sweeps['win_start'] = sweeps.location - (winSize/2)
sweeps['win_start'] = np.where(sweeps.win_start < 0, 0, sweeps.win_start)
sweeps['win_end'] = sweeps.location + (winSize/2)
sweeps['win_end'] = np.where(sweeps.win_end > regionLen, regionLen, sweeps.win_end)

wins = [i for i in range(1, regionLen + winSize, winSize)]
res['genomic_window'] = pd.cut(res['location'], wins, labels=wins[:-1])

res['sweep'] = 0

for i in range(0, len(sweeps)):
    res['sweep'] = np.where((res.location >= sweeps.win_start[i]) & 
                                (res.location <= sweeps.win_end[i]), 1, res.sweep)


l1 = list(res.groupby('genomic_window')['LR'].max())    
l2 = list(res.groupby('genomic_window')['sweep'].sum())

df = pd.DataFrame([l2,l1]).T
df.columns = ['sweep', 'LR']
df['sweep'] = np.where(df.sweep > 1, 1, df.sweep)
#Output df to file
df.to_csv(out_file, sep='\t', index=False, header=True)