import sys

# python3 seq_lines_to_fasta.py myfile.txt
# python3 seq_lines_to_fasta.py myfile.txt 99 ### optional input to start numbering at 99

in_file = sys.argv[1]
outfile = "fasta_of_" + in_file + ".fa"
out_fh = open(outfile, "w")

try:
    int = int(sys.argv[2])
except:
    int = 1

with open(in_file) as f:
    for line in f:
        out_fh.write(">tempname_{}\n".format(int) + line)
        int += 1

out_fh.close()
