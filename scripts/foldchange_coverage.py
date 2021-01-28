#!/usr/bin/python2.6
# -*- coding: utf-8 -*-
#==============================================================================
import argparse
import sys
import os
from collections import defaultdict
import numpy as np
#==============================================================================
#Command line options==========================================================
#==============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("females", type=str,
                    help="A SOAPcov file of scaffolds and cov for females")
parser.add_argument("males", type=str,
                    help="A SOAPcov file of scaffolds and cov for males")
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
def extract_depth(source):
    depthdict = defaultdict(list)
    with open(source, "r") as infile:
        for line in infile:
            line = line.rstrip()
            if not line.startswith("Scaffold"):
                line = line.split(",")
                scaffold = line[0]
                sumdepth = line[2]
                averagedepth = line[3]
                logaveragedepth = line[4]
                depthdict[scaffold].append(sumdepth)
                depthdict[scaffold].append(averagedepth)
                depthdict[scaffold].append(logaveragedepth)
    return depthdict

def combine_depth(females_depth, males_depth):
    combined_dict = defaultdict(list)
    for scaffold in females_depth: 
        fem_sumdepth = float(females_depth[scaffold][0])
        fem_averagedepth = float(females_depth[scaffold][1])
        fem_logaveragedepth = float(females_depth[scaffold][2])
        
        mal_sumdepth = float(males_depth[scaffold][0])
        mal_averagedepth = float(males_depth[scaffold][1])
        mal_logaveragedepth = float(males_depth[scaffold][2])

        #combined sum and averages are identical
        combined_logaveragedepth = mal_logaveragedepth - fem_logaveragedepth

        combined_dict[scaffold].append(combined_logaveragedepth)
        combined_dict[scaffold].append(mal_averagedepth)
        combined_dict[scaffold].append(mal_logaveragedepth)
        combined_dict[scaffold].append(fem_averagedepth)
        combined_dict[scaffold].append(fem_logaveragedepth)
    
    return combined_dict
           
#==============================================================================
#Main==========================================================================
#==============================================================================
def main():
    #extract coverage in males and females
    females_depth = extract_depth(args.females)
    print "Number of scaffolds with coverage in females =", len(females_depth)
    males_depth = extract_depth(args.males)
    print "Number of scaffolds with coverage = males", len(males_depth)

    #check identical scaffolds
    for females in females_depth:
        if females not in males_depth:
            print "error"
    for males in males_depth:
        if males not in females_depth:
            print "error"

    #combine coverage
    combined_depth = combine_depth(females_depth, males_depth)
    print "Number of scaffolds with average coverage =", len(combined_depth)

    with open(args.outfile, "w") as outfile:
        header = "Scaffold, MFLogaverage, Maverage, Mlogaverage, Faverage, Flogaverage"
        outfile.write(header)
        outfile.write("\n")
        for scaffold in combined_depth:
            outfile.write(scaffold)
            outfile.write(",")
            outfile.write(str(combined_depth[scaffold][0]))
            outfile.write(",")
            outfile.write(str(combined_depth[scaffold][1]))
            outfile.write(",")
            outfile.write(str(combined_depth[scaffold][2]))
            outfile.write(",")
            outfile.write(str(combined_depth[scaffold][3]))
            outfile.write(",")
            outfile.write(str(combined_depth[scaffold][4]))
            outfile.write("\n")

if __name__ == '__main__':
    main()
