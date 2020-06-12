#use with python2.7
from __future__ import division #needed if using python2 instead pf python3
from Bio import SeqIO #for parsing the fasta file
import sys
import os 
import argparse

def sequence_searcher(DNA):
	length = len(DNA) - kmer_length + 1
	for x in range(length):
		kmer = DNA[x:x+kmer_length]
		if kmer == search_motif:
            		keep_score[x] = keep_score[x] + 1
		else:
			continue
			
def fasta_extractor(total_fasta_string):
    z = len(total_fasta_string)
    sequence = []
    for x in range(z):
        if total_fasta_string[x] == "A":
            sequence.append("A")
        elif total_fasta_string[x] == "T":
            sequence.append("T")
        elif total_fasta_string[x] == "G":
            sequence.append("G")
        elif total_fasta_string[x] == "C":
            sequence.append("C")
        else:
            continue
    concatenate_bb = ''.join(map(str,sequence))
    return concatenate_bb

parser=argparse.ArgumentParser(
    description='''How to use seq search to find instances of your motif''',
    epilog="""Happy to help!!""")
parser.add_argument('--motif', type=str, default="TTGGG", help='motif to search; should be UPPRERCASE')
parser.add_argument('--filename', type=str, help='file containing the sequences; should be FASTA')
args=parser.parse_args()
if args.motif:
	print("motif being used is " + args.motif)
if args.filename:
	print("reading sequences from " + args.filename)
search_motif = args.motif
kmer_length = len(search_motif)

nucleotide_seq_list = []
#this part reads out the sequences from the fast file and stores in a array
for seq_record in SeqIO.parse(args.filename,"fasta"):
    c = repr(seq_record.seq)
    nucleotide_seq_list.append(c)

#this part initializes a list with 0s; list length is same as length of the input sequences
keep_score = []
ab = nucleotide_seq_list[1]
cd = fasta_extractor(ab)
denominator = len(cd)
for x in range(denominator):
	keep_score.append(0)

#this part takes that array of sequences from previous part and now throws them into the seq counter
alpha = len(nucleotide_seq_list)
for x in range(alpha):
    q = nucleotide_seq_list[x]
    nuc_seq_holder = fasta_extractor(q)
    sequence_searcher(nuc_seq_holder)

for x in range(denominator):
	keep_score[x] = keep_score[x]/alpha
print keep_score

print("My work is done here! Your welcome!! ")	


