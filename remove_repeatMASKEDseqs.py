#use with python2.7
#this script searches through fasta sequences and can be used filter out fasta sequences with excessive number of repeat masked sequences
from __future__ import division #needed if using python2 instead pf python3
from Bio import SeqIO #for parsing the fasta file
import sys
import os 
import argparse

def repeat_counter(total_fasta_string):
    score = 0.0
    z = len(total_fasta_string)
    sequence = []
    for x in range(z):
        if total_fasta_string[x] == "a":
            score = score + 1.0
        elif total_fasta_string[x] == "t":
            score = score + 1.0
        elif total_fasta_string[x] == "g":
            score = score + 1.0
        elif total_fasta_string[x] == "c":
            score = score + 1.0
        else:
            continue
    score = score/z
    return score	
    
parser=argparse.ArgumentParser(
    description='''Hello. Please see help below to use this script to remove sequences with excessive repeat MASKED sequences''',
    epilog="""Happy to help!!""")
parser.add_argument('--cutoff', type=float, default=0.2, help='sequences with the given ratio or higher of masked sequnces will be removed')
parser.add_argument('--filename', type=str, help='file containing the sequences; should be FASTA')
args=parser.parse_args()
if args.cutoff:
	print("cutoff being used is " + str(args.cutoff))
if args.filename:
	print("reading sequences from " + args.filename)
CO = args.cutoff

total = 0
count = 0
output_filename = "filtered.fa"
output_handle = open(output_filename, "w")
#this part reads out the sequences from the fast file and stores in a array
for record in SeqIO.parse(args.filename,"fasta"):
    total = total + 1
    keep_it_or_not = repeat_counter(record)
    if keep_it_or_not <= CO:
	count = count + 1
        SeqIO.write(record, output_handle, "fasta")
output_handle.close()
print(str(count) + " records selected out of " + str(total))

    
