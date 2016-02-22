#!/usr/bin/env python

# Author: Richard Bradshaw, R.T.Bradshaw@soton.ac.uk

# Script to create and calculate SAXS scattering intensities due to individual
# atom-atom scatters. Can then create a difference matrix of atom-atom
# contributions to intensity at a given Q-value, for given reference & test
# conformations.

# imports - MDAnalysis distances, traj readin, numpy, matplotlib
# Requirements: MDAnalysis, numpy, matplotlib
# netcdf4-python for netcdf support (dependency of mdanalysis)

import MDAnalysis.coordinates.TRJ as readamber 
import MDAnalysis.analysis.distances as distanal

from MDAnalysis import Universe
import numpy as np
import sys,argparse
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# argparser


# Exception

class QCalc_Error(Exception):
   """Exception in Debye differences module"""



# Define MDanalysis universe

def setup_universe(parm,traj,fformat='NCDF'):
   """Defines an MDAnalysis Universe object with the desired topology
      and coordinates.
      Usage: setup_universe(parm,traj,[fformat])
      fformat defaults to NetCDF for Amber files. If you're using a
      different format you'll have to set fformat yourself.
      Returns MDAnalysis Universe object"""
   return Universe(parm,traj,format=fformat)

def atoms_slice(univ,atoms='protein'): # Could probably do this more cleverly to allow index based slicing options. Ok for now.
   """Sets up a group of atoms on which to perform later analysis.
      Usage: atoms_slice(universe,[atoms])
      Requires a universe object, and an optional string in CHARMM selection format
      (this is the same as you'd use in VMD). Default is 'protein'.
      Returns MDAnalsyis AtomGroup"""
   return univ.select_atoms(atoms)

def get_trajslice(univ,start=0,end=None,stride=1): # For looping through all conformations
   """Usage: get_trajslice(universe,[start],[end],[stride]).
      Requires a universe object and optional start, end & stride
      Returns a generator for trajectory Timesteps with stride n
      Start = 0,Stride = 1 by default and will return a generator for the whole
      trajectory. Trajectories are zero-indexed."""
   if end is None:
      trajslice = univ.trajectory[start::stride]
   else:
      trajslice = univ.trajectory[start:end:stride]
   return trajslice   
   
def get_coords(atoms,timestep=None):
   """Usage: get_coord(atoms,[timestep]).
      Requires an AtomGroup and optional trajectory timestep object
      Returns a numpy array of shape (n,3) with x,y,z coordinates of atoms at current
      timestep."""

   if timestep is None:
      try:
         return atoms.get_positions()
      except AttributeError:
         raise QCalc_Error("Your atoms don't have a method associated to get their positions. You're probably trying to get the intensity due to a single atom. This won't work with this script")
   else:
      try:
         return atoms.get_positions(ts=timestep)
      except AttributeError:
         raise QCalc_Error("Your atoms don't have a method associated to get their positions. You're probably trying to get the intensity due to a single atom. This won't work with this script")

# Matrix generation
def gen_matrix(i1,atoms):
   """Creates symmetrical matrix from input array of length N(N-1/2).
      Atoms read from atom group provided."""

   addedvals = len(atoms)-1
   curridx = 0
   matx = np.zeros((len(atoms),len(atoms)))
   for i in xrange(len(atoms)):
      matx[i][i+1:] = i1[curridx:curridx+addedvals]
      curridx += addedvals
      addedvals -= 1
   matx = matx + matx.T # Add its own transpose to get symmetry
   return matx



# Object for a single Q-value
class QPoint:

   def __init__(self,q):
      self.q = q

   def get_ff(self,ffnm,atms):
      """Reads in form factor file and deposits into numpy array"""
      
      self._ffs = []
      self._ffarray = []
      with open(ffnm,'r') as f:
         for line in f:
            if line.startswith('#'):
               continue
            else:
               if line.split()[1] == '1': #NOT hydrogens
                  continue
               else:
                  self._ffs.append(float(line.split()[2])) # FF in column 2. All atoms read in.
      if len(self._ffs) != len(atms):
         raise QCalc_Error("Mismatch in the length of form factor file and atoms selected for analysis.")

      for i in self._ffs:
         for j in self._ffs:                   # Double loop over i, j
            self._ffarray.append(i*j)
      self._ffarray = np.reshape(np.array(self._ffarray),(len(atms),len(atms)))

   def apply_debye(self,rij):
      """Calculates Debye factor for spherical scatterers
         sin(Q*rij)/(Q*rij)
         Evaluates to 1 if rij = 0"""
      olderrs = np.seterr(invalid='ignore') # Avoids complaints from nan
      df = np.sin(self.q*rij)/(self.q*rij)
      df[np.isnan(df)] = 1. # Above calc returns nan where rij=0 (i.e. for atomic scattering)
      np.seterr(invalid=olderrs['invalid'])
      return df
   
   def calc_intensities(self,dist):
      """Calculates and returns intensities at this Q value from an array of rij distances."""
 # Need to put tests for size/shape consistency etc in here
      debyearr = np.apply_along_axis(self.apply_debye,1,dist)
      try:
         debyearr = debyearr * np.exp(-0.23*(self.q**2)) * self._ffarray # Modulation parameter
      except ValueError, err:
          print err
          raise QCalc_Error("Error in calculating intensities. Your form factors or distance matrix may be calculated incorrectly. See message above from Numpy")
      return debyearr

def gen_diffmatrix(i1,i2,atoms):
   """Creates difference matrix from two intensity arrays"""

   diffs = i1 - i2
#   addedvals = len(atoms)-1
#   curridx = 0
#   diffmatrix = np.zeros((len(atoms),len(atoms)))
#   for i in xrange(len(atoms)):
#      diffmatrix[i][i+1:] = diffs[curridx:curridx+addedvals]
#      curridx += addedvals
#      addedvals -= 1
   return diffs

def gen_matrix(i1,atoms):
   """Creates symmetrical matrix from input array of length N(N-1/2).
      Atoms read from atom group provided."""

   addedvals = len(atoms)-1
   curridx = 0
   matx = np.zeros((len(atoms),len(atoms)))
   for i in xrange(len(atoms)):
      matx[i][i+1:] = i1[curridx:curridx+addedvals]
      curridx += addedvals
      addedvals -= 1
   matx = matx + matx.T # Add its own transpose to get symmetry
   return matx



def absmax(nparr):
   """Returns the absolute maximum value from a numpy array"""
   _ = np.absolute(nparr)
   return np.max(_)

#### main below here

u = setup_universe('6G08_box.prmtop','6G08_allruns.mdcrd')
prot = atoms_slice(u,atoms='protein and not name H*')
#prot = prot.CA # Edit to only do C-alphas
conf1 = get_coords(prot,u.trajectory[94]) # edit for 95ns frame
conf2 = get_coords(prot,u.trajectory[576]) # edit for 577ns frame
d1 = distanal.self_distance_array(conf1) # No PBC necessary if trajectory already wrapped
d1 = gen_matrix(d1,prot)
d2 = distanal.self_distance_array(conf2)
d2 = gen_matrix(d2,prot)

for qval in [0.1,0.12,0.14,0.16,0.18,0.2,0.22,0.24,0.26,0.28,0.3]:
#   qval = qval*0.01
   Q1 = QPoint(qval)
   Q1.get_ff('6G08_run1_stripped.ff',prot)
   intensity1 = Q1.calc_intensities(d1)
   intensity2 = Q1.calc_intensities(d2)
   print 'Q = %s, intensity at frame 95 = %s' % (qval,np.sum(intensity1))
#   with open('95_plot_mdcrd.dat','a') as f:
#      f.write('%s %s\n' % (qval,np.sum(intensity1)))
   print 'Q = %s, intensity at frame 577 = %s' % (qval,np.sum(intensity2))
#   with open('577_plot.dat','a') as f2:
#      f2.write('%s %s\n' % (qval,np.sum(intensity2)))

   matx = gen_diffmatrix(intensity1,intensity2,prot)
#   matx = matx + matx.T # Add two semi-populated halves together (matrix plus its transpose)

   # Matplotlib fluff
   arrmax = absmax(matx)
   arrmin = arrmax*-1
   plt.figure()
   plt.imshow(matx,interpolation='None',vmin=arrmin/2,vmax=arrmax/2,cmap=cm.RdBu,extent=[1,len(prot),len(prot),1])
   plt.title('NoH_diffs frames 95-577 Q=%s' % qval)
   plt.colorbar()
   plt.savefig('NoH_diffs_95-577_%s.png' % qval,bbox_inches='tight',dpi=300)
   plt.show()
