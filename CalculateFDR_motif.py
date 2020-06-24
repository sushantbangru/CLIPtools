#this script takes the z-score of the k-mers and then calculates FDR of this discovery using method described in mCROSS pub
from __future__ import division #needed if using python2 instead pf python3
from Bio import SeqIO #for parsing the fasta file
import sys
import os 
import argparse
import pandas as pd
import csv
from itertools import izip

def fdr_calculator(test_value):
	a = 1
	b = 0
	anti_value = 0 - test_value
	for x in range(num_motif):
		if float(S1[x+1]) > test_value:
			a = a + 1
		elif float(S1[x+1]) < anti_value:
			b = b + 1
	q =  b/a
	return q

parser=argparse.ArgumentParser(
    description='''How to use FDR finder for calculating FDR of z-scores''',
    epilog="""Happy to help!!""")
parser.add_argument('--nrows', type=int, help='number of samples for which z-score lists exist in the table')
parser.add_argument('--filename', type=str, help='file containing the z-scores; should be csv; please list columns as motif,S1,S2,S3,...')
parser.add_argument('--filetosave', type=str, help='filename to save the output in; give it as filename.csv')
args=parser.parse_args()
writfile = args.filetosave
print writfile
if args.nrows:
	print("number of individual z-score rows are " + str(args.nrows))
if args.filename:
	print("reading z-scores from " + args.filename)

col_list = ["motif","S1","S2","S3","S4"]
df = pd.read_csv(args.filename, names = col_list)
S1 = df.S1.tolist()
S2 = df.S2.tolist()
S3 = df.S3.tolist()
S4 = df.S4.tolist()
motif = df.motif.tolist()

S1_FDR = ["S1-FDR"]
S2_FDR = ["S2-FDR"]
S3_FDR = ["S3-FDR"]
S4_FDR = ["S4-FDR"]	
num_motif = len(S1) - 1
#print num_motif

for x in range(num_motif):
	if float(S1[x+1]) > 0:
		FDR = fdr_calculator(float(S1[x+1]))
		S1_FDR.append(str(FDR))
	else:
		S1_FDR.append("1")

for x in range(num_motif):
	if float(S2[x+1]) > 0:
		FDR = fdr_calculator(float(S2[x+1]))
		S2_FDR.append(str(FDR))
	else:
		S2_FDR.append("1")

for x in range(num_motif):
	if float(S3[x+1]) > 0:
		FDR = fdr_calculator(float(S3[x+1]))
		S3_FDR.append(str(FDR))
	else:
		S3_FDR.append("1")

for x in range(num_motif):
	if float(S4[x+1]) > 0:
		FDR = fdr_calculator(float(S4[x+1]))
		S4_FDR.append(str(FDR))
	else:
		S4_FDR.append("1")

with open(writfile, 'wb') as f:
    writer = csv.writer(f)
    writer.writerows(izip(motif, S1_FDR, S2_FDR, S3_FDR, S4_FDR))

print("job complete! Thanks!")
