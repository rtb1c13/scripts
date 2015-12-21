#!/usr/bin/env python

# Script to perform analyses for file of HFE data.
# Will calculate mean error, mean absolute deviation (MUE),
# Kendall's tau, R**2 and print out to file
# Author: Richard Bradshaw, R.T.Bradshaw@soton.ac.uk

import numpy as np
import scipy.stats as stats
import scikits.bootstrap as boot
import sys,argparse

### Argparser
def parse():
   parser = argparse.ArgumentParser()
   parser.add_argument("-f","--files",help="File names for analysis (start with file used for t-tests)",nargs='+',required=True)
   args = parser.parse_args()
   return args

### Functions below here
def MUE(data):
   """Calculates mean absolute deviation
      (mean unsigned error) for entries in array.
      Returns MUE and std error in MUE"""

   absdev = []
   for i in data:
      absdev.append(np.absolute(i[1] - i[0]))
   MUE = np.mean(absdev)
   MUEse = stats.sem(absdev)
   return MUE, MUEse

def MSE(data):
   """Calculates mean signed error for entries in array.
      Returns MSE and std error in MSE"""

   signdev = []
   for i in data:
      signdev.append(i[1] - i[0])
   MSE = np.mean(signdev)
   MSEse = stats.sem(signdev)
   return MSE, MSEse

def correls(data):
   """Calculates Pearson's R and Kendall's tau
      for entries in array.
      
      Returns:
      slope, intercept, r_value, p_value, std_err
      (= values from linear regression)
      tau, taup
      (= values from Kendall's tau)
      as entries in dictionary correldict"""

   slope, intercept, r_value, p_value, std_err = stats.linregress(data[:,0],data[:,1])
   correldict = {'slope': slope,'intercept': intercept,'r_value': r_value, 'p_value': p_value,'std_err': std_err}
   tau, taup = stats.kendalltau(data[:,0],data[:,1])
   correldict['tau'] = tau
   correldict['taup'] = taup
   return correldict 

def correls_for_bootstrap(data):
   """Identical to correls but suitable for bootstrapping
      
      Returns:
      slope, intercept, r_value, p_value, std_err
      (= values from linear regression)
      tau, taup
      (= values from Kendall's tau)
      NOT as dictionary"""

   slope, intercept, r_value, p_value, std_err = stats.linregress(data[:,0],data[:,1])
   r2 = r_value**2
   tau, taup = stats.kendalltau(data[:,0],data[:,1])
   return slope,intercept,r_value,r2,p_value,std_err,tau,taup 

def write_data(fn,data):
   """Performs descriptive stats and writes stats to output file"""

   f = open(fn,'w')
   mue,muese = MUE(data)
   f.write("Errors are 95% CIs\n")
   f.write("MUE = %5.3f +/- %5.3f\n" % (mue,muese*1.96))
   mse,msese = MSE(data)
   f.write("MSE = %5.3f +/- %5.3f\n" % (mse,msese*1.96))
   correldict = correls(data)
   f.write("R^2 = %3.2f\n" % correldict['r_value']**2)
   f.write("K-Tau = %3.2f\n\n" % correldict['tau'])
   f.write("BOOTSTRAPPED RESULTS (10k resamples, 95% CIs)\n")
   CIs = boot.ci(data,MUE)
   f.write("MUE = %5.3f < %5.3f < %5.3f\n" % (CIs[0][0],mue,CIs[1][0]))
   CIs = boot.ci(data,MSE)
   f.write("MSE = %5.3f < %5.3f < %5.3f\n" % (CIs[0][0],mse,CIs[1][0]))
   CIs = boot.ci(data,correls_for_bootstrap)
   f.write("Pearson's R = %3.2f < %3.2f < %3.2f\n" % (CIs[0][2],correldict['r_value'],CIs[1][2]))
   f.write("R^2 = %3.2f < %3.2f < %3.2f\n" % (CIs[0][3],correldict['r_value']**2,CIs[1][3]))
   f.write("K-Tau = %3.2f < %3.2f < %3.2f\n\n" % (CIs[0][6],correldict['tau'],CIs[1][6]))
   f.close()

def t_test_diffs(data,data2,outfile):
   """Will calc p values between MSE and MUE distributions
      for 2 datasets"""

   f = open(outfile,'a')
   f.write("Student's t-test between error distributions of %s and %s\n" %(args.files[0].split(".")[0],outfile.split(".")[0]))
   signdev = []
   for i in data:
      signdev.append(i[1] - i[0])
   signdev2 = []
   for i in data2:
      signdev2.append(i[1] - i[0])
   tstat,p = stats.shapiro(signdev)
   tstat,p2 = stats.shapiro(signdev2)
   f.write("Signed err Shapiro-Wilks (normality) p-orig, p-current = %6.4f, %6.4f\n" % (p,p2))
   tstat,p = stats.ttest_ind(signdev,signdev2)
   f.write("Signed err Independent p = %6.4f\n" % p)
   tstat,p = stats.ttest_rel(signdev,signdev2)
   f.write("Signed err Related p = %6.4f\n" % p)
   tstat,p = stats.wilcoxon(signdev,signdev2)
   f.write("Signed err Wilcoxon p = %6.4f\n" % p)
   tstat,p = stats.mannwhitneyu(signdev,signdev2)
   f.write("Signed err Mann-Whitney U p = %6.4f\n" % p)
   absdev = []
   for i in data:
      absdev.append(np.absolute(i[1] - i[0]))
   absdev2 = []
   for i in data2:
      absdev2.append(np.absolute(i[1] - i[0]))
   tstat,p = stats.shapiro(absdev)
   tstat,p2 = stats.shapiro(absdev2)
   f.write("Abs err Shapiro-Wilks (normality) p-orig, p-current = %6.4f, %6.4f\n" % (p,p2))
   tstat,p = stats.ttest_ind(absdev,absdev2)
   f.write("Abs err Independent p = %6.4f\n" % p)
   tstat,p = stats.ttest_rel(absdev,absdev2)
   f.write("Abs err Related p = %6.4f\n" % p)
   tstat,p = stats.wilcoxon(absdev,absdev2)
   f.write("Abs err Wilcoxon p = %6.4f\n" % p)
   tstat,p = stats.mannwhitneyu(absdev,absdev2)
   f.write("Abs err Mann-Whitney U p = %6.4f\n" % p)
   f.close()

def calc_anova(anova,outfile):
   """Will calc One-way Anova for each entry
      in array 'Anova', divided up into groups of
      3 repeats"""
   f = open(outfile,'w')
   f.write("Ligand   p-value\n")
   for lig in enumerate(anova,start=1):
      groups = np.reshape(lig[1],(len(args.files),3))
      anova_args = []
      [ anova_args.append(i) for i in groups ]
      f_stat,p = stats.f_oneway(*anova_args)
      f.write("%8s %8.4f\n" %(lig[0],p))
   f.close()

### Main below here

def main():
   global args
   args = parse()
   # Run through desc stats for file 1
   prefix = args.files[0].split(".")[0] 
   fn = prefix+".txt"
   outfile = prefix+"_results.txt"
   data = np.genfromtxt(fn)
#   ligindex = [4,5,6,7,8,27,28,29,30,31,32,33,34] # Can be used to analyse a subset of ligands
#   data = data[ligindex]
#   print data
   data = np.hsplit(data,[1]) # Split expt values apart
   expt = data[0]
   comput = data[1]
   anova = np.copy(comput) # Beginning of array for anova
   origdata = np.hstack((expt,np.reshape(np.mean(comput,axis=1),(47,1)))) # Reshape then stack, now each ligand has an entry of length 2
   write_data(outfile,origdata)
   # Now do stats and t test for others
   for argfile in args.files[1:]:
      prefix = argfile.split(".")[0] 
      fn = prefix+".txt"
      outfile = prefix+"_results_vs_file1.txt"
      data = np.genfromtxt(fn)
#      data = data[ligindex]
      data = np.hsplit(data,[1]) # Split expt values apart
      expt = data[0]
      comput = data[1]
      anova = np.hstack((anova,comput)) # Add data to Anova array
      currdata = np.hstack((expt,np.reshape(np.mean(comput,axis=1),(47,1)))) # Reshape then stack, now each ligand has an entry of length 2
      write_data(outfile,currdata)
      t_test_diffs(origdata,currdata,outfile)
   # Now do Anova      
   calc_anova(anova,"Anovas.txt")

   
### Main below here
if __name__ == "__main__":
   main()
