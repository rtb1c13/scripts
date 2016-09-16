#!/usr/bin/env python

# Script to calculate widths of SAXS metadynamics CVs based on
# relative errors in expt measurements. Currently, CV widths
# are 2*expt error, centred at 0 (perfect agreement with expt)

import numpy as np

# Function to return index of closest value in an array
def get_closest(arr,val):
   diffs = np.abs(arr-val)
   return diffs.argmin()

# Interpolate a relative error from the two nearest values in a  Q/intensity/error array 
def calc_relative(arr,val):
   qdiffs = arr[:,0]-val
   idx1 = np.argmax(qdiffs[qdiffs<0]) # -ve value closest to 0
   idx2 = idx1 + 1
   relatives = arr[:,2] / arr[:,1]
   y = relatives[idx2] - relatives[idx1]
   x = qdiffs[idx2] - qdiffs[idx1]
   grad = y/x
   return relatives[idx2] - (grad*qdiffs[idx2])

#def write_relatives( 


# Read in the expt datafile and the expt profile calculated with plumed

expt = np.genfromtxt("expt/6G11_newunits_merged.dat",dtype="float32")
with open("single_metad_fixedCV/Q0015/profile/profile-cv1.out",'r') as f:
   profile = f.readlines()[1:51]

for line in profile:
   q,r_int = [float(x) for x in line.split()]
   r_err = calc_relative(expt,q)
   int_err = r_int*r_err
   with open("plumed_inputs.txt",'a') as wf:
      wf.write("Q=%5.3f INTERVAL=%5.3f,%5.3f GRID_MIN=%5.3f GRID_MAX=%5.3f GRID_BIN=1000\n"
               % (q,-3*int_err,3*int_err,-3.6*int_err,3.6*int_err))


