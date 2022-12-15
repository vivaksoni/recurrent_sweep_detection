#script takes .ms file and creates allele frequency file for SweepFinder2
#python3  get_aff_for_SF2.py -msFile "sweep_detection/results/stationary/benProp_0.05_2Nes_10/1Mb_rep1.ms" \
#-fixedFile "sweep_detection/results/stationary/benProp_0.05_2Nes_10/1Mb_rep1.fixed" \
#-outFile  "sweep_detection/scripts/SF_inputFiles/results/stationary/benProp_0.05_2Nes_10/1Mb_rep1.aff" \
#-regionLen 997204 -samples 100 

from __future__ import print_function
import sys
import pandas as pd
import math
import os
import argparse
import numpy as np


#parsing user given constants
parser = argparse.ArgumentParser(description='Information about length of region and sample size')
parser.add_argument('-msFile', dest = 'msFile', action='store', nargs = 1, type = str, help = 'path to ms file (.ms format)')
parser.add_argument('-fixedFile', dest = 'fixedFile', action='store', nargs = 1, type = str, help = 'path to fixed file (.fixed format)')
parser.add_argument('-outFile', dest = 'outFile', action='store', nargs = 1, type = str, help = 'path to output file')
parser.add_argument('-regionLen', dest = 'regionLen', action='store', nargs = 1, type = int, help = 'size of simulated region')
parser.add_argument('-samples', dest = 'samples', action='store', nargs = 1, type = int, help = 'no. of samples in .ms file')

args = parser.parse_args()
chr_len =  args.regionLen[0]
samples = args.samples[0]
ms_file = args.msFile[0]
f_fixed = args.fixedFile[0]
out_file = args.outFile[0]


#Function to parse .ms file data into nested list of positions and genotypes
#Function takes as input open reading file object and no. of samples
def get_nested_data_list(f_ms, f_fixed, samples, chr_len):
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
        x = int(np.round(float(pos)*chr_len,0))
        l_Pos.append(int(x))
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
    
    #Sum genotype information to get DACs
    for i,j in enumerate(l_data):
        l_data[i].append(sum([int(x) for x in l_data[i][1]]))
    
    #Convert dictionary to dataframe, removing genotype column
    df = pd.DataFrame(l_data)[[0,2]]
    #Rename columns
    df.columns = ['position', 'x']
    #Add samples info
    df['n'] = samples
    #Add fold info
    df['folded'] = 0

    #Read in .fixed file
    fixed = pd.read_csv(f_fixed, skiprows=2, sep=' ',
                   names=['tempID', 'permID', 'mutType', 'position', 's', 'h', 'initial_subpop', 'origin_gen',
                         'fix_gen'])

    #Create dataframe of fixed mutations
    f_df = pd.DataFrame([[x,samples,samples,0] for x in fixed.position])
    f_df.columns = df.columns
    #Concatenate dfs and sort by position
    aff = pd.concat([df, f_df]).sort_values('position').reset_index(drop='True')
    #Remove duplicates, keeping the higher value (ie in case the beneficial has overlayed a previous mutation)
    aff = aff.groupby('position', group_keys=False).apply(lambda x: x.loc[x['x'].idxmax()]).reset_index(drop=True)
    
    return(aff)


f_ms = open(ms_file, 'r')
df = get_nested_data_list(f_ms, f_fixed, samples, chr_len)

df.to_csv(out_file, sep='\t', header=True, index=False)

print ("Allele frequency file output to " + out_file)