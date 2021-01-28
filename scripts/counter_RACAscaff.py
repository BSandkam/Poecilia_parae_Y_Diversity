#!/usr/bin/python2.6
# -*- coding: utf-8 -*-
#==============================================================================
import argparse
import sys
import os
from collections import defaultdict
import collections
#==============================================================================
#Command line options==========================================================
#==============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("target_segs", type=str,
                    help="target species segments.refined.txt output file from RACA")
parser.add_argument("ptarget", type=str,
                    help="name used for the target species when running RACA, e.g. pwin, ppic, ppar")
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

    scaff_dicty = []
    # reads in the target species segments.refined.txt output file from RACA, extracts the scaffold names and counts how many times each scaffold is listed
    with open(args.target_segs, "r") as target_segs:
        for line in target_segs:
            if line.startswith(target):
                scaff = line.split(target + ".")[1].split(":")[0]
                scaff_dicty.append(scaff)
    counter_scaffs = collections.Counter(scaff_dicty)
    
    # writes each scaffold and number of times it appeared in the target species segments.refined.txt file
    with open(args.outfile, "w") as out:
        for scaff in counter_scaffs:
            out.write(scaff)
            out.write(",")
            out.write(str(counter_scaffs[scaff]))
            out.write("\n")


if __name__ == '__main__':
    main()
