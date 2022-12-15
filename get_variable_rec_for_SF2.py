from __future__ import print_function
import sys
import pandas as pd
import os
import argparse
import numpy as np
import pandas as pd

parser = argparse.ArgumentParser(description='Information about length of region and sample size')
parser.add_argument('-affFile', dest = 'affFile', action='store', nargs = 1, type = str, help = 'path to aff file')
parser.add_argument('-recFile', dest = 'recFile', action='store', nargs = 1, type = str, help = 'rec rate file')
parser.add_argument('-outFile', dest = 'outFile', action='store', nargs = 1, type = str, help = 'path to output file')

args = parser.parse_args()
f_aff = args.affFile[0]
f_rec = args.recFile[0]
out_file = args.outFile[0]


#Read in aff file
df = pd.read_csv(f_aff, header=0, sep='\t')
#Read in rr file
rr = pd.read_csv(f_rec, names=['pos','ro'], sep='\t')

#Create nested list of ro at each position
l = [[x*10**-6] * 10000 for x in rr.ro]
#Flattern nested list
lst = [item for sublist in l for item in sublist]

#Create list for results, starting with a 0
res = [0]
#Loop through positions
for i in range(0, len(df.position)-1):
    #Set indexes to calculate sum of rates between polymorphisms
    p2 = df.position[i]
    p1 = df.position[i+1] 
    res.append(np.sum(lst[p2:p1]))

df['rate'] = res
df = df[['position', 'rate']]

#Output to file
df.to_csv(out_file, sep='\t', header=True, index=False)