"""

Author: Benjamin Furman
This script is designed to take a fasta file and create a mate pair library for each of the sequences in the file. Take a file of sequences, digest them into a sets of mate pairs with some sort of insert size, and make a lot of those for each sequence.

At minimum you need to specify a fasta file (-f), insert size (-in), length of sequences for ends of mate pairs (-seq), and the step size betwen start points of mate pairs (-step). You can also choose to exclud scaffolds < a certain length (-exculde), although any scaffold less than the length of 1 mate pair (inset + seqLen*2) is automatically excluded. A minimum of two mate pairs per scaffold will be made (to ensure both ends are captured), unless a scaffold is exactly the length of the desired matepair.

You can also have a fastq output of a made up quality (-fq quality), you can have the output gzipped (-gz), you can not have the mate pair made for the end of the scaffold so generated mate pairs are exactly evenly spaced (-noend). By default the first read in the pair will be reverse complimented, but you can set '-no_revcomp' to stop that action.

examples:

python3.7 mate_pair_maker.py -f test.fa -in 100 -seq 25 -step 50

# change output file name to myRun_mp_1.fa and myRun_mp_2.fa
python3.7 mate_pair_maker.py -f test.fa -in 100 -seq 25 -step 50 -out myRun

# allows 10% of the seqeuences from the tips of the mate pairs to be Ns
python3.7 mate_pair_maker.py -f test.fa -in 100 -seq 25 -step 50 -tip_N 0.1

# allows 10% of the seqeuences from the tips of the mate pairs to be Ns and 15% of the middle of the mate pair
python3.7 mate_pair_maker.py -f test.fa -in 100 -seq 25 -step 50 -tip_N 0.1 -mid_N 0.15

# adds an individual identifier to the fasta headers e.g., changes >1 to >mySeq_1
python3.7 mate_pair_maker.py -f test.fa -in 100 -seq 25 -step 50 -tip_N 0.1 -indv mySeq

# gzip output
python3.7 mate_pair_maker.py -f test.fa -in 100 -seq 25 -step 50 -gz

# gzip and fastq output with madeup "~" quality
python3.7 mate_pair_maker.py -f test.fa -in 100 -seq 25 -step 50 -gz -fq "~"

"""

def get_cmdline_args():
	import argparse
	parser = argparse.ArgumentParser()
	parser.add_argument("-f", "--fasta_file", help = "the fasta file of sequences you want digested into mate pairs (gzipped or not)", required = True)

	parser.add_argument("-in", "--insert_size", help = "insert size of mate pairs", required = True, type = int)

	parser.add_argument("-seq", "--seq_length", help = "length of sequencing data you want from the ends of the mate pairs", default = 0, required = True, type = int)

	parser.add_argument("-step", "--step_size", help = "how far apart to put the starts of each mate pair. If this is less than the mate pair inset_size+seq_length, then you will get overlapping kmers", required = True, type = int)

	parser.add_argument("-exclude", "--exclude_sequence_length", help = "ignore input sequences less than this length", type = int)

	parser.add_argument("-fq", "--fastq_output", help = "can export fastq with a made up quality as specified by this parameter (e.g., -fq ~)")

	parser.add_argument("-indv", "--individual_ID", help = "individual identifier to be added to fasta sequence header")

	parser.add_argument('-noend', '--no_fromEnd_matepair', help = "by default, the script ensures that the end of a sequence is turned into a mate pair, regardless of how it fits into the step size. This option turns that off.", action = 'store_true')

	parser.add_argument('-gz', '--gzip_output', action = 'store_true', help = "gzips the matepair output files.")

	parser.add_argument("-tip_N", "--prop_tip_Ns", default = 1.0, help = "proportion of Ns in the sequences from the ends of the mate pairs", type = float)

	parser.add_argument("-mid_N", "--prop_Ns_middle", default = 1.0, help = "proportion of Ns in the middle of the mate pair", type = float)

	parser.add_argument("-no_revcomp", "--no_revcomp_p1", action = 'store_true', default = 'store_false', help = "default the script will revcomp the first seq in the read pair. Add this flag to stop that.")

	parser.add_argument("-rev2", "--revcom_read2", action = 'store_true', default = 'store_false', help = "default the script will revcomp the first seq in the read pair. Add this flag to revcomp the second read in the pair. Add -no_revcomp and -rev2 to revcomp the second and not the first (so a normal library construction).")

	parser.add_argument("-out", "--outfile_basename", help = "will be added to the output matepair file names.", default = "matepairs")

	# parser.add_argument('-run_N', "--run_of_Ns", help = "if set, this will remove SCAFFOLDS that have a run of Ns in the middle larger than this value.")

	args = parser.parse_args()
	return(args)

def fasta_to_dict(file):
	"""reads fasta file to a dictionary"""
	Dict = {}
	collect = 0
	import re
	import gzip
	from os.path import exists
	if not exists(file) :
		exit("{} not found".format(file))
	else:
		if ".gz" in file:
			IN = gzip.open(file, 'rt')
		else:
			IN = open(file)

	for line in IN.readlines():
		line = line.strip()
		if ">" in line and collect == 1:
			collect = 0
		if not ">" in line and collect == 1:
			Dict[name] = Dict[name]+line
		if ">" in line and collect == 0:
			collect = 1
			name = re.sub('>','',line)
			name = name.split()[0]
			Dict[name] = ''
	IN.close()
	return(Dict)


def remove_short_contigs(dict, min_lens):
	"""removes short contigs"""
	for seq in list(dict.keys()):
		if sum([True for min in min_lens if len(dict[seq]) < min]) > 0:
			del dict[seq]


def check_Ns(seqs, max_threshold):
	"""Takes list of sequences and determines if they have too many Ns in them. Only total, not sequential (though they could also be sequential)."""

	for seq in seqs:
		if seq.count('N')/len(seq) > max_threshold:
			return(False)

	return(True)

# def check_mid_Ns(seq, max_threshold):
# 	"""Takes the whole make pair seq and checks for the number of Ns in that seq.""" ## would be more accurate to look for runs of N of a certain size.
# 	pass

def make_a_mate_pair(seq, seq_len, no_revcomp = False, revcomp_read2 = False):
	""" take a chunk of sequence, spits out each end of seq_len size"""
	p1 = seq[0:seq_len]
	p2 = seq[-seq_len:]

	if no_revcomp is True:
		if revcomp_read2 is True:
			p2 = rev_comp(p2)
			return([p1,p2])
		else:
			return([p1,p2])
	else:
		p1 = rev_comp(p1)
		if revcomp_read2 is True:
			p2 = rev_comp(p2)
			return([p1,p2])
		else:
			return([p1,p2])

def make_mate_pairs(seq, seq_len, insert_size, step_size, prop_tip_Ns, prop_middle_Ns, no_revcomp = False, revcomp_read2 = False):
	"""takes a sequence and rolls through and makes sets of make pairs"""
	total_mate_length = seq_len*2+insert_size
	p1_seqs = list()
	p2_seqs = list()
	for i in range(0, len(seq)-total_mate_length+1, step_size):
		pairs = make_a_mate_pair(seq[i:i+total_mate_length], seq_len, no_revcomp, revcomp_read2)

		if ( check_Ns(pairs, prop_tip_Ns) and
				check_Ns([seq[i:i+total_mate_length]], prop_middle_Ns)
			):
			p1_seqs.append(pairs[0])
			p2_seqs.append(pairs[1])

	return(p1_seqs, p2_seqs)

def check_if_seq_end_covered(target_seq_length, mate_pair_len, step_size):
	"""confirms if the end of the sequence was made into a mate pair"""
	if (target_seq_length - mate_pair_len) % step_size > 0:
		return(False)
	else:
		return(True)

def make_pair_of_end(seq, seq_len, insert_size, prop_tip_Ns, prop_middle_Ns, no_revcomp = False, revcomp_read2 = False):
	"""captures the end of a sequences as a mate pair."""
	total_mate_length = seq_len*2+insert_size
	pairs  = make_a_mate_pair(seq[-total_mate_length:], seq_len, no_revcomp, revcomp_read2)
	if ( check_Ns(pairs, prop_tip_Ns) and
			check_Ns([seq[-total_mate_length:]], prop_middle_Ns)
		):
		return(pairs)
	else:
		return([])


def rev_comp(seq):
	"""Will take a seq and spit back the rev comp of it"""
	complements = {
		"A" : "T", "T" : "A", "C" : "G", "G" : "C", "." : ".", "-" : "-",
		"Y" : "R", "R" : "Y", "M" : "K", "K" : "M", "W" : "W", "S": "S", "N" : "N"
	}
	new_seq = list()
	for i in seq[::-1]:
		try:
			new_seq.append(complements[i.upper()])
		except:
			exit("unknown base pair '{}' in your input sequences".format(i))
			# new_seq.append("N")

	return("".join(new_seq))


def open_output_files(fastq_output, gzip_file, base_name):
	"""get normal or gzip files ready"""
	import gzip
	from os.path import exists
	base_name = base_name + "_mp"
	if fastq_output:
		file_ext = ".fq"
	else:
		file_ext = ".fa"

	if gzip_file: # arguments online suggest that zipping this was is not efficient, but zipping after is better (of course). So probably want to repeat at the end the zipping, but this will help with lots of mate pairs per seq.
		file_ext = file_ext + ".gz"

	p1_file_name = '{}_1{}'.format(base_name, file_ext)
	p2_file_name = '{}_2{}'.format(base_name, file_ext)
	if exists(p1_file_name) or exists(p2_file_name):
		exit("mate pair files exits in this directory, I won't run here and overwrite them")

	if gzip_file:
		p1_fileH = gzip.open(p1_file_name,'wt')
		p2_fileH = gzip.open(p2_file_name,'wt')
	else:
		p1_fileH = open(p1_file_name,'w')
		p2_fileH = open(p2_file_name,'w')

	return(p1_fileH, p2_fileH)


def rezip_file(file_name):
	"""rezips a file handle for best compression"""
	import os
	os.system("gzip -q {} > file_name".format(file_name))

def write_to_output(fh, seqs, name, pair, fq):
	"""output seqs with names to file"""

	if fq is None:
		start = ">"

		names = [start + name + "_" + str(val+1) + "_" + pair for val in range(0, len(seqs))]

		for x in zip(names, seqs):
			fh.write(x[0] + "\n" + x[1] + "\n")
	else:
		start = "@"

		names = [start + name + "_" + str(val+1) + "_" + pair for val in range(0, len(seqs))]

		qual_names = ["+" + name + "_" + str(val+1) + "_" + pair for val in range(0, len(seqs))]

		qualities = [fq*len(seq) for seq in seqs]

		for x in zip(names, seqs, qual_names, qualities):
			fh.write(x[0] + "\n" + x[1] + "\n" + x[2]+ "\n" + x[3] + "\n")

def message(scaff = 0, total = 0, message = 0):
	import sys
	import getpass
	username = getpass.getuser()

	import random
	num = random.randint(0,3)
	lyrics = (
	"\n\n\t\t...hold tight {}, I'm still working on it\n\n".format(username),
	"\n\n\t\tI think I can, I think I can, I think I can... generate your pseudomate pairs\n\n",
	"\n\n\t\tI dont know what you heard about me, \n\n\t\tBut I cant get a mate pair assembly,\n\n\t\tNo scaffolds, no consensus, you can't see,\n\n\t\tSo lets infer them with P-M-P \n\n\t\t--Dave 'Cp-G MC' Metzger\n\n"
	)

	if message == 0 or message == 1:
		print(lyrics[num], file = sys.stderr)

	if message is 1:
		perc = round((scaff / total)*100)
		print("\n\n\t\t...working on scaffold {} of {} ({}% done)...\n\n".format(scaff, total, perc), file = sys.stderr)

	if message is 2:
		print("\n\n\t\t...loading the scaffolds...\n\n", file = sys.stderr)
	if message is 3:
		print("\t\t...done loading {} scaffolds...\n\n".format(total), file = sys.stderr)
	if message is 4:
		print("\n\n\t\trezipping the files for optimal compression\n\n")


def check_py_version():
	"""this is pointless as it dies from syntax errors anyway if it's python version 2"""
	import sys
	vers = sys.version_info
	if sys.version_info<(3,0):
		exit("sorry, you need to run this with some version of python3")

def main():
	import warnings
	from datetime import datetime

	check_py_version()

	message()
	silly_num = 0
	scaff_count = 0

	args = get_cmdline_args()

	message(message = 2)
	seqs_dict = fasta_to_dict(args.fasta_file)
	total_scaffolds = len(seqs_dict.keys())
	message(message = 3, total = total_scaffolds)

	mate_pair_len = args.insert_size + (args.seq_length * 2)

	if args.seq_length > args.insert_size:
		exit("end seq longer than inset size, no longer a mate pair. Dying now.")

	if args.step_size > mate_pair_len:
		warnings.warn("given your step size and desired mate pair length, you'll have non-overlapping mate pairs")

	if args.step_size == mate_pair_len:
		warnings.warn("given your step size and desired mate pair length, you'll have mate pairs that cover the same sequence twice (read 2 of first mate pair is read 1 of second mate pair)")

	if args.exclude_sequence_length:
		remove_short_contigs(seqs_dict, [args.exclude_sequence_length, mate_pair_len])

		if len(seqs_dict.keys()) == 0:
			exit("there were no input sequences longer than your 'exclude' length. No output made.")

	p1_fileH, p2_fileH = open_output_files(args.fastq_output, args.gzip_output, args.outfile_basename)

	for name,seq in seqs_dict.items():
		scaff_count += 1
		### just printing a silly message
		if datetime.now().second % 13 == 0 and silly_num == 1:
				pass
		elif datetime.now().second % 13 == 0 and silly_num == 0:
			message(scaff = scaff_count, total = total_scaffolds, message = 1)
			silly_num = 1
		else:
			silly_num = 0

		if args.individual_ID:
			new_name = name + "_" + args.individual_ID
		else:
			new_name = name

		p1_seqs, p2_seqs = make_mate_pairs(seq, args.seq_length,
		 									args.insert_size, args.step_size,
											args.prop_tip_Ns,
											args.prop_Ns_middle,
											args.no_revcomp_p1,
											args.revcom_read2
											)

		if ( check_if_seq_end_covered(len(seq),
									mate_pair_len, args.step_size) == False
									and
			 args.no_fromEnd_matepair == False
			):

			end_pairs = make_pair_of_end(seq, args.seq_length,
										args.insert_size,
										args.prop_tip_Ns,
										args.prop_Ns_middle,
										args.no_revcomp_p1,
										args.revcom_read2
										)

			if len(end_pairs) > 0:
				p1_seqs.append(end_pairs[0])

				p2_seqs.append(end_pairs[1])

		write_to_output(p1_fileH, p1_seqs, new_name, "1", args.fastq_output)
		write_to_output(p2_fileH, p2_seqs, new_name, "2", args.fastq_output)

	## rezip the final files for better compression
	p1_fileH.close()
	p2_fileH.close()
	if args.gzip_output:
		message(message = 4)
		rezip_file(p1_fileH.name)
		rezip_file(p2_fileH.name)
	message(scaff = total_scaffolds, total = total_scaffolds, message = 1)

if __name__ == '__main__':
	main()
