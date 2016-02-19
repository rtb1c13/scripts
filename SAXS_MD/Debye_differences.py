#!/usr/bin/env python

# Script to create and calculate SAXS scattering intensities due to individual
# atom-atom scatters. Can then create a difference matrix of atom-atom
# contributions to intensity at a given Q-value, for two given conformations.

# imports - MDAnalysis distances, traj readin, numpy, matplotlib
# Requires netcdf4-python for netcdf support (dependency of mdanalysis)
import MDAnalysis.coordinates.TRJ as readamber 
import MDAnalysis.analysis.distances as distanal

from MDAnalysis import *
import numpy as np
import argparse
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# argparser



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

def get_trajslice(univ,start=0,end=None,stride=1):
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
      return atoms.get_positions()
   else:
      return atoms.get_positions(ts=timestep)





# Object for a single Q-value
class QPoint:

   def __init__(self,q):
      self.q = q

   def get_ff(self,ffnm):
      """Reads in form factor file and deposits into numpy array"""
      
      self._ffs = []
      self._ffarray = []
      with open(ffnm,'r') as f:
         for line in f:
            if line.startswith('#'):
               continue
            else:
#               if line.split()[1] == '1': #NOT hydrogens
#                  continue
#               else:
                  self._ffs.append(float(line.split()[2])) # FF in column 2. All atoms read in.
      for idx,val in enumerate(self._ffs):
         for j in self._ffs[idx+1:]:                   # Double loop over i, j>i
            self._ffarray.append(val*j)
      self._ffarray = np.array(self._ffarray)
      self._ffarray = self._ffarray*2  # Correct for sum over i=1, j=1

#   def read_conformers(self,conf1,conf2=None):


   # meths to read trajectory, extract a frame
   # Get formfactors
   # Write PDB?
   # Get distmatrix

# class difference

   # Given 2 conformers, calc distmatrix
   # Write formatted distmatrix
   # Plot (somehow? Matplotlib?)

   def apply_debye(self,rij):
      """Performs first part of Debye formula for spherical scatterers
         sin(Q*rij)/(Q*rij)"""

      return np.sin(self.q*rij)/(self.q*rij)
   
   def calc_intensities(self,dist):
      """Calculates and returns intensities at this Q value from an array of rij distances."""
 # Need to put tests for size/shape consistency etc in here
     
      debyedist = np.apply_along_axis(self.apply_debye,0,dist)
      return debyedist * self._ffarray      

def gen_diffmatrix(i1,i2,atoms):
   """Creates difference matrix from two intensity arrays"""

   diffs = i1 - i2
   addedvals = len(atoms)-1
   curridx = 0
   diffmatrix = np.zeros((len(atoms),len(atoms)))
   for i in xrange(len(atoms)):
      diffmatrix[i][i+1:] = diffs[curridx:curridx+addedvals]
      curridx += addedvals
      addedvals -= 1
   return diffmatrix

def absmax(nparr):
   """Returns the absolute maximum value from a numpy array"""
   _ = np.absolute(nparr)
   return np.max(_)

#### main below here

u = setup_universe('6G08_box.prmtop','6G08_allruns.mdcrd')
prot = atoms_slice(u,atoms='resid 1-10')
prot = prot.CA # Edit to only do C-alphas
conf1 = get_coords(prot,u.trajectory[94]) # edit for 95ns frame
conf2 = get_coords(prot,u.trajectory[576]) # edit for 577ns frame
d1 = distanal.self_distance_array(conf1,box=u.trajectory[94].dimensions)
d2 = distanal.self_distance_array(conf2,box=u.trajectory[576].dimensions)

for qval in [0.10,0.12,0.14,0.16,0.18,0.20,0.22,0.24,0.26,0.28,0.30]:
   Q1 = QPoint(qval)
   Q1.get_ff('6G08_run1_CA_1-10.ff')
   intensity1 = Q1.calc_intensities(d1)
   intensity2 = Q1.calc_intensities(d2)
   print 'Q = %s, intensity at frame 95 = %s' % (qval,sum(intensity1))
   print 'Q = %s, intensity at frame 577 = %s' % (qval,sum(intensity2))

   matx = gen_diffmatrix(intensity1,intensity2,prot)
   matx = matx + matx.T # Add two semi-populated halves together (matrix plus its transpose)

   # Matplotlib fluff
#   arrmax = absmax(matx)
#   arrmin = arrmax*-1
#   plt.figure()
#   plt.imshow(matx,interpolation='None',vmin=arrmin,vmax=arrmax,cmap=cm.RdBu,extent=[1,len(prot),len(prot),1])
#   plt.title('NoH_diffs frames 95-577 Q=%s' % qval)
#   plt.colorbar()
#   plt.savefig('NoH_diffs_95-577_%s.png' % qval,bbox_inches='tight',dpi=600)
#   plt.show()
