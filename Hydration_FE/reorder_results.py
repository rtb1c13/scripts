#!/usr/bin/env python

# Script to print out ligand files in format ready for
# R anova & Tukey HSD analysis
# 
# Author: Richard Bradshaw, R.T.Bradshaw@soton.ac.uk

# Requirements: Numpy

import numpy as np
import sys,argparse

### Argparser
def parse():
   parser = argparse.ArgumentParser()
   parser.add_argument("-f","--files",help="File names for analysis (start with file used for t-tests)",nargs='+',required=True)
   args = parser.parse_args()
   return args


def main():
   global args
   args = parse()
   # create list of file prefixes
   prefixes=[]
   for _ in args.files:
      prefixes.append(_.split(".")[0])
   for idx,fn in enumerate(args.files):
      with open(fn,'r') as readresult:
         for lig,line in enumerate(readresult,1):
            vals = line.split()
            with open("lig"+str(lig)+".txt",'a') as outfile:
               outfile.write("%s %s\n" % (vals[1],prefixes[idx])) # Don't need expt result, start at 1
               outfile.write("%s %s\n" % (vals[2],prefixes[idx]))
               outfile.write("%s %s\n" % (vals[3],prefixes[idx]))


### Main below here
if __name__ == "__main__":
   main()
