import sys
import argparse
import tempfile
from Bio.Seq import Seq

#This script takes a Jellyfish dump file and adds the reverse compliment of this in a new column

parser=argparse.ArgumentParser(description="Add reverse compliment column to jelly dump -d dump -o output")
parser.add_argument('-d', '--dump', type=str, required=True, help="Provide a Jelly_dump file")
parser.add_argument('-o', '--output', action='store', help="Directs the output to a name of your choice")
args = parser.parse_args()

with open(args.output, 'w') as out_file:
	with open(args.dump, 'r') as d_file:
		for line in d_file:
			line = line.rstrip('\n').rstrip('\t').split( )
			kmer = Seq(line[0])
			revcomp = kmer.reverse_complement()
			line_out = ' '.join([str(kmer), str(revcomp)]) + '\n'
			out_file.write(line_out)
