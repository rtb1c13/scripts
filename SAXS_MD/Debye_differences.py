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

def parse():
   parser = argparse.ArgumentParser()
   parser.add_argument("-p","--parm",help="Parameter/topology file for desired structure",type=str,required=True)
   parser.add_argument("-t","--traj",help="Coordinate file for desired structure",type=str,required=True)
   parser.add_argument("-ff","--formfactor",help="Formfactor file for desired structure. Note that formfactors for Hydrogens should be included in their attached heavy atoms. This script does not currently support explicit hydrogens",type=str,required=True)
   parser.add_argument("-f","--frames",help="Frames to analyse. If you supply one frame here, a SAXS profile will be generated. Two frames will generate a difference plot. Three frames will be used as the start, stop & stride values through a trajectory, and generate a profile for each one (WARNING, TIME CONSUMING!). Defaults to calculating the profile of the first frame",nargs='+',type=int,default=1)
   parser.add_argument("-q","--qvalues",help="List of Q-values for analysis. Defaults to 0.1",nargs='*',type=float,default=[0.1])
   parser.add_argument("-s","--select",help="Selection string for atoms to analyse from trajectory (in CHARMM/VMD format). Formfactor file should be pre-processed to only include these atoms.",type=str,default='protein and not name H*')
 

   if len(sys.argv)==1:
      parser.print_help()
      sys.exit(1)

   args = parser.parse_args()
   return args


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

def atoms_slice(univ,atoms='protein and not name H*'): # Could probably do this more cleverly to allow index based slicing options. Ok for now.
   """Sets up a group of atoms on which to perform later analysis.
      Usage: atoms_slice(universe,[atoms])
      Requires a universe object, and an optional string in CHARMM selection format
      (this is the same as you'd use in VMD). Default is 'protein and not name H*'.
      Returns MDAnalysis AtomGroup"""
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
#def gen_matrix(i1,atoms):
#   """Creates symmetrical matrix from input array of length N(N-1/2).
#      Atoms read from atom group provided."""
#
#   addedvals = len(atoms)-1
#   curridx = 0
#   matx = np.zeros((len(atoms),len(atoms)))
#   for i in xrange(len(atoms)):
#      matx[i][i+1:] = i1[curridx:curridx+addedvals]
#      curridx += addedvals
#      addedvals -= 1
#   matx = matx + matx.T # Add its own transpose to get symmetry
#   return matx



# QPoint class for a single Q-value
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
      """Calculates and returns intensities at this Q value from an array of
         rij distances."""
      debyearr = np.apply_along_axis(self.apply_debye,1,dist)
      try:
         debyearr = debyearr * np.exp(-0.23*(self.q**2)) * self._ffarray # Modulation parameter
      except ValueError, err:
          print err
          raise QCalc_Error("""Error in calculating intensities. Your 
                               form factors or distance matrix may be 
                               calculated incorrectly. See message above
                               from Numpy""")
      return debyearr

# Conformer class for single frames of a trajectory
class Conformer:
   """Class to hold information for SAXS analysis of single structure"""
   def __init__(self,univ,frame=None,**kwargs):
      self.atms = atoms_slice(univ,**kwargs)
      if frame is not None:
         self.coords = get_coords(self.atms,univ.trajectory[frame]) 
      else:
         self.coords = get_coords(self.atms)
      self.intensities={}
      self.label = univ.trajectory.frame + 1

   def _gen_matrix(self):
      """Creates symmetrical matrix from input array of length N(N-1/2).
         Atoms read from atom group of current object."""

      addedvals = len(self.atms)-1
      curridx = 0
      matx = np.zeros((len(self.atms),len(self.atms)))
      for i in xrange(len(self.atms)):
         matx[i][i+1:] = self._dists[curridx:curridx+addedvals]
         curridx += addedvals
         addedvals -= 1
      matx = matx + matx.T # Add its own transpose to get symmetry
      return matx

   def get_distmatrix(self,**kwargs):
      """Generates a distance matrix using the self_distance_array function
         of MDAnalysis."""

      self._dists = distanal.self_distance_array(self.coords,**kwargs) # No PBC necessary if 
                                                                       # trajectory already wrapped
      self.distmatx = self._gen_matrix()


# Functions for plotting
def absmax(nparr):
   """Returns the absolute maximum value from a numpy array"""
   _ = np.absolute(nparr)
   return np.max(_)

def write_profile(conf,qmin=0.01,qmax=0.50,step=0.01):
   """Writes a SAXS profile of a given conformer to file.
      Usage: write_profile(Conformer,[qmin],[qmax],[step])
      Requires Conformer object as input, by default prints out
      profile between Q=0.01 and Q=0.50"""

   for q in np.arange(qmin,qmax+step,step):
      try:
         conf.intensities[q]
         continue
      except KeyError:
         Q1 = QPoint(q)
         Q1.get_ff(args.formfactor,conf.atms)
         conf.intensities[q] = Q1.calc_intensities(conf.distmatx)
   with open(args.traj+'.f%s.dat' % conf.label,'a') as f:
      for q in np.arange(qmin,qmax+step,step):
         f.write('%5.3f %12.6f\n' % (q,np.sum(conf.intensities[q])))

def plot_diffmatrix(ref,probe,q=0.1):
   """Plots 2D difference matrix for given Conformer objects
      at desired Q-value (default 0.1).
      Usage: plot_diffmatrix(Conformer1,Conformer2,[q])"""
   try:
      diffs = ref.intensities[q] - probe.intensities[q]
   except ValueError, err:
      print err
      raise QCalc_Error("""Your intensity arrays are different sizes - 
                           are you sure your structures/atom selections
                           are the same size? See Numpy error above.""")
   except KeyError:
      raise QCalc_Error("""One or both of your conformers is missing an 
                           intensity matrix.""")
   arrmax = absmax(diffs)
   arrmin = arrmax*-1
   plt.figure()
   plt.imshow(diffs,interpolation='None',vmin=arrmin/2,vmax=arrmax/2,cmap=cm.RdBu)
#   plt.gca().set_xticks(range(1,len(ref.atms)+1,200))
#   plt.gca().set_yticks(range(len(ref.atms),0,-200))
   plt.suptitle('Intensity difference matrix, frames %s-%s Q=%s' % (ref.label,probe.label,q),y=0.98,fontsize=14)
   plt.title('Atoms selected = "%s", zero indexed' % args.select,fontsize=10)
   plt.colorbar()
   plt.savefig('Diffs_%s-%s_%s.png' % (ref.label,probe.label,q),bbox_inches='tight',dpi=300)
      

############ main below here ################
def main():
   global args
   args = parse()
   u = setup_universe(args.parm,args.traj)
   if len(args.frames) == 1:                   # Profile only
      conf1 = Conformer(u,frame=args.frames[0]-1,atoms=args.select)
      conf1.get_distmatrix()
      write_profile(conf1)
   elif len(args.frames) == 2:                 # Difference matrix
      conflist = []
      for conf in args.frames:
         currconf = Conformer(u,frame=conf-1,atoms=args.select) 
         currconf.get_distmatrix()
         for q in args.qvalues:
            Q1 = QPoint(q)
            Q1.get_ff(args.formfactor,currconf.atms)
            currconf.intensities[q] = Q1.calc_intensities(currconf.distmatx)
         conflist.append(currconf)
      for q in args.qvalues:
         plot_diffmatrix(conflist[0],conflist[1],q)
   else:                                       # Trajectory slice
      trajslice = get_trajslice(u,args.frames[0],args.frames[1],args.frames[2])
      for ts in trajslice:
         conf1 = Conformer(u,frame=ts.frame-1,atoms=args.select)
         conf1.get_distmatrix()
         write_profile(conf1)
###############################################################      
main()      
