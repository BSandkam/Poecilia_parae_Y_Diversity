#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
USAGE: kmers-BAM2BED.py SAM/BAM
-----
'''

import sys
import re
import pysam


global CIGAR_PAT
CIGAR_PAT = re.compile(r'\d+[MIDNSHP=XB]{1}')


def usage():
	print(__doc__)
	sys.exit(1)


def cigar2tuple(cigar):
	counts = {}
	for centry in CIGAR_PAT.findall(cigar):
		try:
			counts[centry[-1]] += int(centry[:-1])
		except KeyError:
			counts[centry[-1]] = int(centry[:-1])
	return tuple(counts.items())


if __name__ == '__main__':

	if len(sys.argv) != 2:
		usage()

	mode = 'rb'
	infile = sys.argv[1]
	if infile.endswith('.sam'):
		mode = 'rt'
	elif infile == '-':  # assume is SAM (also needs header)
		mode = 'rt'

	sam = pysam.AlignmentFile(infile, mode)
	for read in sam:

		if read.is_unmapped:
			continue

		strand = '+'
		seq = read.query_sequence
		if read.is_reverse:
			strand = '-'
			seq = read.get_forward_sequence()

		t_cigar = read.cigartuples
		t_cigar_len = sum(x[1] for x in t_cigar)

		# clipped bases do not count for edit distance (NM)
		# so we need to add these to compare against the alternative alignments
		# in pysam, operations 4 and 5 are SOFT and HARD clip respectively
		clipped = sum(x[1] for x in t_cigar if x[0] in (4, 5))

		best_nm = read.get_tag('NM') + clipped
		sys.stdout.write(
			'{}\t{}\t{}\t{}\t{}\t{}\n'.format(
				read.reference_name, read.reference_start, read.reference_start + t_cigar_len,
				seq, best_nm, strand))

		# assume k-mer is mapped non-uniquely, so by default we get all alignments
		# kinda deals with N and M
		unique = False
		try:
			xt = read.get_tag('XT')
			if xt == 'U':
				unique = True
		except KeyError:
			pass

		try:
			xa = read.get_tag('XA').split(';')
			if not unique:
				for xal in xa:
					if xal:
						rname, rpos, cigar, nm = xal.split(',')
						alt_cigar = cigar2tuple(cigar)
						alt_cigar_len = sum(x[1] for x in alt_cigar)
						alt_clipped = sum(x[1] for x in alt_cigar if x[0] in ('S', 'H'))

						nm = int(nm) + alt_clipped

						if nm <= best_nm:
							if rpos.startswith('+') or rpos.startswith('-'):
								rpos = rpos[1:]
							rpos = int(rpos) - 1  # XA positions are 1-based
							sys.stdout.write(
								'{}\t{}\t{}\t{}\t{}\t{}\n'.format(
									rname, rpos, rpos + alt_cigar_len,
									seq, nm, strand))
		except KeyError:
			pass
		except AttributeError:
			pass
