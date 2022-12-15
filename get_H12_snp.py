#Script calculates Garud's H12.

#Example usage:
#python3  get_H12_snp.py \
#-msFile "sweep_detection/results/stationary/benProp_0.05_2Nes_10/1Mb_rep1.ms" \
#-outFile "sweep_detection/scripts/garud_H12/stationary/benProp_0.05_2Nes_10/1Mb_rep1.txt" \
#-regionLen 997204 -samples 100 -winSize 1000

import argparse
import sys
import pandas as pd
import math
import os
import numpy as np
from scipy import stats
import allel

#parsing user given constants
parser = argparse.ArgumentParser(description='Information about length of region and sample size')
parser.add_argument('-msFile', dest = 'msFile', action='store', nargs = 1, type = str, help = 'path to ms file')
parser.add_argument('-outFile', dest = 'outFile', action='store', nargs = 1, type = str, help = 'path to store output file')
parser.add_argument('-winSize', dest = 'winSize', action='store', nargs = 1, type = int, help = 'size of sliding windows')
parser.add_argument('-regionLen', dest = 'regionLen', action='store', nargs = 1, type = int, help = 'size of simulated region')
parser.add_argument('-samples', dest = 'samples', action='store', nargs = 1, type = int, help = 'no. of samples in .ms file')

args = parser.parse_args()
regionLen =  args.regionLen[0]
winSize = int(args.winSize[0]/2)
samples = args.samples[0]
msFile = args.msFile[0]
outFile = args.outFile[0]


#Function to parse .ms file data into nested list of positions and genotypes
#Function takes as input open reading file object and no. of samples
def get_nested_data_list(f_ms, samples):
    l_Pos = [] #list of positions of SNPs
    l_Genos = [] #list of alleles
    d_tmp = {} #dict to store individual allele info for each individual (values) at each site (keys)

    #positions on line 2
    pos_lines = [2]
    for position, line in enumerate(f_ms):
        if position in pos_lines:
            #store positions in list
            pos_list  = line.split()
    #Set file pointer to start of file
    f_ms.seek(0)

    i = 0
    #Loop through positions, storing in list
    for pos in pos_list[1:]:
        #Append position to l_Pos (after converting to float)
        l_Pos.append(float(pos))
        #Add dictionary key for each position, with empty value
        d_tmp[str(i)] = ""
        i += 1 
        
    
    #genotypes on line 3 onwards (use samples argument to determine length of file)
    g_lines = [x for x in range(3, samples + 4)]
    #Loop through lines (ie individuals)
    for position, line in enumerate(f_ms):
        if position in g_lines:
            #Remove newline character
            line1 = line.strip('\n')
            i = 0
            #For each individual, loop through each site, appending allele information for that individual to 
            #the site number in dict
            while i < len(line1):
                d_tmp[str(i)] = d_tmp[str(i)] + line1[i]
                i = i + 1

    f_ms.seek(0)

    #Create nested list of positions and genotypes
    l_data = [[j, d_tmp[str(i)]] for i,j in enumerate(l_Pos)]
    return(l_data)


#Read in .ms file
f_ms = open(msFile, 'r')
l_data = get_nested_data_list(f_ms, samples)

snps = [x[1] for x in l_data]
snps = [[float(x) for x in y] for y in snps]
snps_pos = [int(np.round(x[0]*regionLen))-1 for x in l_data]

hapArr = np.zeros(shape=[regionLen,samples])

for i,j in enumerate(snps_pos):
    hapArr[j] = snps[i]

hapArr = allel.HaplotypeArray(hapArr, dtype='i1')

#Create df of snp postions and windows around snps
df = pd.DataFrame([[x, x-winSize, x+winSize] for x in snps_pos], columns=['snp_position', 'win_start', 'win_end'])
df['win_start'] = np.where(df.win_start < 0, 0, df.win_start)
df['win_end'] = np.where(df.win_end > regionLen, regionLen, df.win_end)

lst = []
#Loop through df, subsetting haploytype array and estimating H12 on subsetted array
for x in range(0, len(df)):
    subhap = hapArr[df.win_start[x]:df.win_end[x]]
    lst.append(allel.garud_h(subhap)[1])

df['H12'] = lst

df.to_csv(outFile, sep='\t', header=True, index=False)

