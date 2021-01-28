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
parser.add_argument("scaffolds_raca", type=str,
                    help="A file (output from 2.scaff_posinfo_RACA.py) with target scaffolds and their position in PCFs and P. ret chromosomes")
parser.add_argument("cov_fc", type=str,
                    help="A file with coverage fold change for each scaffold (output from 07.cov_fold_change.py)")
parser.add_argument("outfile", type=str,
                    help="A file with target scaffolds with RACA position and associated coverage estimates")
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

    # makes dictionary of scaffolda and their positional information in PCFs and P. ret chromosomes
    scaff_pos_info = defaultdict(list)
    with open(args.scaffolds_raca, "r") as scaffolds:
        next(scaffolds)
        for line in scaffolds:
            line = line.rstrip()
            scaffold = line.split(",")[0]
            scaff_start = line.split(",")[1]
            scaff_end = line.split(",")[2]
            pcf = line.split(",")[3]
            pcf_start = line.split(",")[4]
            pcf_end = line.split(",")[5]
            pret_chromo = line.split(",")[6]
            pret_chromo_start = line.split(",")[7]
            pret_chromo_end = line.split(",")[8]
            scaff_pos_info[scaffold].append(scaff_start)
            scaff_pos_info[scaffold].append(scaff_end)
            scaff_pos_info[scaffold].append(pcf)
            scaff_pos_info[scaffold].append(pcf_start)
            scaff_pos_info[scaffold].append(pcf_end)
            scaff_pos_info[scaffold].append(pret_chromo)
            scaff_pos_info[scaffold].append(pret_chromo_start)
            scaff_pos_info[scaffold].append(pret_chromo_end)
    print "Number of target species' scaffolds with RACA positional information = ", len(scaff_pos_info)
    
    # outputs a file containing only the scaffolds in the dictionary above, together with their coverage and positional information
    count = 0
    with open(args.outfile, "w") as out:
        with open(args.cov_fc, "r") as cov:
            for line in cov:
                if line.startswith("Scaffold,"):
                    line = line.rstrip()
                    out.write(line)
                    out.write(",")
                    out.write("PCF")
                    out.write(",")
                    out.write("PCF_start")
                    out.write(",")
                    out.write("PCF_end")
                    out.write(",")
                    out.write("Pret_LG")
                    out.write(",")
                    out.write("Pret_LG_start")
                    out.write(",")
                    out.write("Pret_LG_end")
                    out.write("\n")
                else:
                    scaff = line.split(",")[0]
                    if scaff in scaff_pos_info:
                        count += 1
                        line = line.rstrip()
                        MFLogav = line.split(",")[1]
                        Mav = line.split(",")[2]
                        Mlogav = line.split(",")[3]
                        Fav = line.split(",")[4]
                        Flogav = line.split(",")[5]
                        PCF = scaff_pos_info[scaff][2]
                        PCF_start = scaff_pos_info[scaff][3]
                        PCF_end = scaff_pos_info[scaff][4]
                        LG = scaff_pos_info[scaff][5]
                        LG_start = scaff_pos_info[scaff][6]
                        LG_end = scaff_pos_info[scaff][7]
                        out.write(scaff)
                        out.write(",")
                        out.write(str(MFLogav))
                        out.write(",")
                        out.write(str(Mav))
                        out.write(",")
                        out.write(str(Mlogav))
                        out.write(",")
                        out.write(str(Fav))
                        out.write(",")
                        out.write(str(Flogav))
                        out.write(",")
                        out.write(str(PCF))
                        out.write(",")
                        out.write(str(PCF_start))
                        out.write(",")
                        out.write(str(PCF_end))
                        out.write(",")
                        out.write(str(LG))
                        out.write(",")
                        out.write(str(LG_start))
                        out.write(",")
                        out.write(str(LG_end))
                        out.write("\n")
    print "Number of scaffolds with RACA positional information and coverage information = ", count



if __name__ == '__main__':
    main()