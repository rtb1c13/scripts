#!/usr/bin/env python

# Matplotlib histogram of fields

import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-f","--files",help="Files for fields at first atom",nargs='+',type=str,required=True)
parser.add_argument("-f2","--filestwo",help="Files for fields at second atom",nargs='+',type=str)

args = parser.parse_args()



def combine_files(files):
   """Reads in a list of files and returns concatenated
      Numpy array."""
   cat = np.array([])
   for fn in files:
      x = np.loadtxt(fn)
      cat = np.concatenate((cat,x))
   return cat

def histogram(values):
   """Shows histogram of results in values"""
   n, bins, patches = plt.hist(values, 50, normed=1, color='blue')
   plt.show()

def desc_stats(values,desc):
   """Writes out some descriptive statistics of the values
      read in. Prepends output with 'desc'"""
   xbar = np.mean(values)
   sigma = np.std(values)
   error = stats.sem(values)
   with open("stats.txt",'a') as f:
#      f.write("Description    Mean      Std. Dev.  Std. Err.\n")
      f.write("%s %7.3f %7.3f    %7.3f\n" % (desc.ljust(14),xbar,sigma,error))

#def combined_stats(values1,prefix):
#   """Writes out stats for combination of values read in
#      e.g. average field between two, 

def combine_atoms(val1,val2):
   """Combines two arrays of fields at given atoms
      and returns the average field and field drop
      (average and difference of values in arrays)"""
   ave = (val1 + val2) / 2
   diff = val2 - val1
   return ave,diff

def main():
   atm1 = combine_files(args.files)
   with open("stats.txt",'w') as f:
      f.write("Description    Mean      Std. Dev.  Std. Err.\n")
   desc_stats(atm1,"atm1")
   if args.filestwo is not None:
      atm2 = combine_files(args.filestwo)
      desc_stats(atm2,"atm2")
      Fvib,dFvib = combine_atoms(atm1,atm2)
      desc_stats(Fvib,"Ave field")
      desc_stats(dFvib,"Field drop")
      histogram(Fvib)
   else:
      histogram(atm1)

main()

#y = mlab.normpdf( bins, mu, sigma)
#l = plt.plot(bins, y, 'r--', linewidth=1)

#From tuts
#plt.xlabel('Smarts')
#plt.ylabel('Probability')
#plt.title(r'$\mathrm{Histogram\ of\ IQ:}\ \mu=100,\ \sigma=15$')
#plt.axis([40, 160, 0, 0.03])
#plt.grid(True)

#plt.show()
