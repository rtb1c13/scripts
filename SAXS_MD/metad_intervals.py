#!/usr/bin/env python

# Script to calculate widths of SAXS metadynamics CVs based on
# relative errors in expt measurements. Currently, CV widths
# are 2*expt error, centred at 0 (perfect agreement with expt)

import numpy as np
import argparse,sys

### Argparser ###
def parse():
   parser = argparse.ArgumentParser()
   parser.add_argument("-p","--plumed",help="Path to Plumed SAXS profile",type=str,required=True)
   parser.add_argument("-e","--expt",help="Path to experimental SAXS profile",type=str,required=True)
   parser.add_argument("-l","--length",help="Number of Q-values to read in from Plumed profile (i.e. length). Default=50",type=int,default=50)
   parser.add_argument("-q","--qvalues",help="List of Q-values for analysis. Defaults to 0.1",nargs='*',type=float,default=[0.1])
   parser.add_argument("-o","--output",help="Filename for output. Defaults to 'plumed_outputs.txt'",type=str,default="plumed_outputs.txt")
   parser.add_argument("-s","--scale",help="Scale factor for width of MetaD interval. Multiplied by the relative experimental error at each Q-value. Defaults to 3.0",type=float,default=3.0)
   parser.add_argument("-a","--adjust",help="Adjust MetaD intervals using a reference set of intensities (the second entry in the plumed intensity files. Defaults to False",type=bool,default=False)

   if len(sys.argv)==1:
      parser.print_help()
      sys.exit(1)

   args = parser.parse_args()
   return args

### Class for SAXS profile ###
class Profile:
   """Class to handle reading/writing (TBC) of 
      plumed-calculated SAXS profiles"""

   def __init__(self,fn,qvals=50,ref=False):
      """Reads in plumed SAXS profile of given filename
         and number of Q-values (default 50). Usage:
         Profile(filename,[qvals=50])"""
      self._inpfile = fn
      self.plumed = np.genfromtxt(fn,skip_header=1,max_rows=qvals,dtype="float32")
      if ref:
         self.xtal = np.genfromtxt(fn,skip_header=3+qvals,max_rows=qvals,dtype="float32")
         

   def get_expt(self,fn):
      """Reads in experimental profile for given Plumed profile
         ready for calculation of relative errors etc. Usage:
         get_expt(filename)"""
      self.expt = np.genfromtxt(fn, dtype="float32")

   def _closest(self,arr,val):
      """ Function to return index of closet value in array.
          Currently unused."""
      diffs = np.abs(arr-val)
      return diffs.argmin()

# Interpolate a relative error from the two nearest values in a  Q/intensity/error array 
   def _itplerr(self,val):
      """Interpolate a relative error from the two nearest values in
         experimental Q/intensity/error array"""
      qdiffs = self.expt[:,0]-val
      idx1 = np.argmax(qdiffs[qdiffs<0]) # -ve value closest to 0
      idx2 = idx1 + 1
      relatives = self.expt[:,2] / self.expt[:,1]
      y = relatives[idx2] - relatives[idx1]
      x = qdiffs[idx2] - qdiffs[idx1]
      grad = y/x
      return relatives[idx2] - (grad*qdiffs[idx2])

# Interpolate an intensity from the two nearest values in a Q/intensity array 
   def _itplint(self,arr,val):
      """Interpolate a relative error from the two nearest values in
         experimental Q/intensity/error array"""
      qdiffs = arr[:,0]-val
      try:
         idx1 = np.argmax(qdiffs[qdiffs<0]) # -ve value closest to 0
      except ValueError: # No negative values
         if np.min(qdiffs) == 0:
            idx1 = np.argmin(qdiffs)
         else:
            raise ValueError("No intensity values at Q < %4.3f. Interpolation not possible" % val)
            sys.exit(1) 
      idx2 = idx1 + 1
      y = arr[idx2,1] - arr[idx1,1]
      x = qdiffs[idx2] - qdiffs[idx1]
      grad = y/x
      return arr[idx2,1] - (grad*qdiffs[idx2])

   def get_rel_err(self,qvals,adj_intervals=False):
      """Calculates the relative error in a plumed intensity at a
         given Q value from an associated experimental dataset"""
      self.relatives = []
      for q in qvals:
         r_int = self._itplint(self.plumed,q)
         r_err = self._itplerr(q)
         self.relatives.append([q,r_int*r_err])
      self.relatives.sort()

   def write_plumed(self,fn="plumed_outputs.txt",sc=3,adj_intervals=False):
      """Writes plumed inputs at desired Q values, scaled
         according to scale factor""" 
      for _ in self.relatives:
         if adj_intervals:
            ref_int = self._itplint(self.xtal,_[0])
            r_int = self._itplint(self.plumed,_[0])
            with open(fn,'a') as wf:
               wf.write("Q=%5.3f INTERVAL=%5.3f,%5.3f GRID_MIN=%5.3f GRID_MAX=%5.3f GRID_BIN=1000\n"
                        % (_[0],-sc*_[1] - (r_int-ref_int),
                        sc*_[1] - (r_int-ref_int),
                        -sc*1.2*_[1] - (r_int-ref_int),
                        sc*1.2*_[1] - (r_int-ref_int)))
         else:
            with open(fn,'a') as wf:
               wf.write("Q=%5.3f INTERVAL=%5.3f,%5.3f GRID_MIN=%5.3f GRID_MAX=%5.3f GRID_BIN=1000\n"
                        % (_[0],-sc*_[1],sc*_[1],-sc*1.2*_[1],sc*1.2*_[1]))
            

### Write out a series of Plumed inputs ###
def generate_plumed_input():
   prof = Profile(args.plumed,args.length,ref=args.adjust)
   prof.get_expt(args.expt)
   prof.get_rel_err(args.qvalues)
   prof.write_plumed(args.output,args.scale,adj_intervals=args.adjust)


### Main below here ###
if __name__ == "__main__":
   global args
   args = parse()
   generate_plumed_input()
