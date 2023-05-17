#!/usr/bin/python2.6
# -*- coding: utf-8 -*-
#==============================================================================
import argparse
import sys
import os
from collections import defaultdict
import numpy as np
import collections
#==============================================================================
#Command line options==========================================================
#==============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("ptarget", type=str,
                    help="name used for the target species when running RACA, e.g. piwn, ppic, ppar")
parser.add_argument("conserved_segments", type=str,
                    help="Conserved segments output file from RACA")
parser.add_argument("scaff_counter", type=str,
                    help="")
parser.add_argument("target_segs", type=str,
                    help="target species segments.refined.txt output file from RACA")
parser.add_argument("outfile", type=str,
                    help="")
# This checks if the user supplied any arguments. If not, help is printed.
if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)
# Shorten the access to command line arguments.
args = parser.parse_args()
#==============================================================================
#Functions=====================================================================
#==============================================================================
           
#==============================================================================
#Main==========================================================================
#==============================================================================
def main():

    target = args.ptarget

    unique_scaffs = []
    dupl_scaffs = []
    # extracts uniquely mapping scaffolds from the RACA output
    with open(args.scaff_counter, "r") as counter:
        for line in counter:
            line = line.rstrip()
            scaffold = line.split(",")[0]
            counts = line.split(",")[1]
            if counts == "1":
                unique_scaffs.append(scaffold)
            else:
                dupl_scaffs.append(scaffold)
    print "Number of uniquely mapping scaffolds = ", len(unique_scaffs)
    print "Number of nonuniquely mapping scaffolds = ", len(dupl_scaffs)

    # this section reads the target species segments.refined.txt output file from RACA and makes a dictionary of 
    # scaffold name, scaffold start (bp in scaffold total length), scaffold end (bp in scaffold total length), pcf for that scaffold, 
    # scaffold start in pcf, scaffold end in pcf
    scaff_length = defaultdict(list)
    with open(args.target_segs, "r") as target_segs:
        for line in target_segs:
            if line.startswith("RACA"):
                pcf = str(line.split("RACA.")[1].split(":")[0])
                start_pcf = line.split(":")[1].split("-")[0]
                end_pcf = line.split("-")[1].split()[0]
            elif line.startswith(target):
                scaffold = line.split(target + ".")[1].split(":")[0]
                if scaffold in unique_scaffs:
                    start_scaff = line.split(":")[1].split("-")[0]
                    end_scaff = line.split("-")[1].split()[0]
                    if scaffold not in scaff_length:
                        scaff_length[scaffold].append(start_scaff)
                        scaff_length[scaffold].append(end_scaff)
                        scaff_length[scaffold].append(pcf)
                        scaff_length[scaffold].append(start_pcf)
                        scaff_length[scaffold].append(end_pcf)
                    else:
                        print "ERROR - scaffold wrongly placed in the uniquely mapping scaffolds list", scaffold
                        sys.exit()


    # this section find the start and stop positions of uniquley mapping scaffolds in P reticulata chromosomes
    scaff_chromo = defaultdict(list)
    with open(args.conserved_segments , "r") as segs:
        for line in segs:
            if line.startswith("xiph"):
                xiph_chromo = line.split("xiph.")[1].split(":")[0]
                xiph_chromo_start = line.split(":")[1].split("-")[0]
                xiph_chromo_end = line.split("-")[1].split()[0]
            elif line.startswith(target):
                pwin_scaffold = line.split(target + ".")[1].split(":")[0]
                pwin_scaffold_start = line.split(":")[1].split("-")[0]
                pwin_scaffold_end = line.split("-")[1].split()[0]
                if pwin_scaffold in scaff_length:
                    if scaff_length[pwin_scaffold][0] == pwin_scaffold_start and scaff_length[pwin_scaffold][1] == pwin_scaffold_end:
                        scaff_chromo[pwin_scaffold].append(pwin_scaffold_start)
                        scaff_chromo[pwin_scaffold].append(pwin_scaffold_end)
                        scaff_chromo[pwin_scaffold].append(scaff_length[pwin_scaffold][2])
                        scaff_chromo[pwin_scaffold].append(scaff_length[pwin_scaffold][3])
                        scaff_chromo[pwin_scaffold].append(scaff_length[pwin_scaffold][4])    
                        scaff_chromo[pwin_scaffold].append(xiph_chromo)
                        scaff_chromo[pwin_scaffold].append(xiph_chromo_start)
                        scaff_chromo[pwin_scaffold].append(xiph_chromo_end)

    with open(args.outfile, "w") as out:
        out.write("scaffold,scaffold_start,scaffold_end,pcf,pcf_start,pcf_end,xiph_chromo,xiph_chromo_start,xiph_chromo_end\n")
        for scaff in scaff_chromo:
            out.write(scaff)
            out.write(",")
            out.write(str(scaff_chromo[scaff][0]))
            out.write(",")
            out.write(str(scaff_chromo[scaff][1]))
            out.write(",")
            out.write(scaff_chromo[scaff][2])
            out.write(",")
            out.write(str(scaff_chromo[scaff][3]))
            out.write(",")
            out.write(str(scaff_chromo[scaff][4]))
            out.write(",")
            out.write(scaff_chromo[scaff][5])
            out.write(",")
            out.write(str(scaff_chromo[scaff][6]))
            out.write(",")
            out.write(str(scaff_chromo[scaff][7]))
            out.write("\n")
 


if __name__ == '__main__':
    main()
