#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Shorten fasta headers/ids to a given length

Author: Pedro Almeida
-----
'''

import argparse
import sys, gzip


def fasShortHeads(fasta, maxlength=10):

	# read file 1st time to determine possible names
	seen = {}
	with fopen(fasta, 'rt') as f:
		viter = (l.rstrip("\r\n") for l in f)
		for line in viter:
			if line.startswith('>'):
				name = line.split()[0][1:maxlength + 1]
				if name not in seen:
					print('>{}'.format(name))
					seen[name] = 1
				else:
					print('# ERROR:')
					print('Found conflicting names: {} -> {}'.format(line, name))
					print('Correct these before continuing or increase --length')
					sys.exit(1)
			else:
				print(line)


def fopen(file, mode):
	''' Opens a (gzip) file in open mode '''
	if file == '-':
		fobj = sys.stdin
	elif file.lower().endswith('.gz'):
		fobj = gzip.open(file, mode)
	else:
		fobj = open(file, mode)

	return fobj


def get_parser():
	''' get parser object for script '''
	parser = argparse.ArgumentParser(
			description=__doc__,
			formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument(
		'fasta', metavar='FASTA',
		help='input fastA file (reads from pipe with -)')
	parser.add_argument(
		'-l', '--length', metavar='INT', type=int,
		default=10,
		help='maximum length of fasta headers [%(default)s]')

	args = parser.parse_args()
	return args


if __name__ == '__main__':
	args = get_parser()
	fasShortHeads(args.fasta, args.length)
