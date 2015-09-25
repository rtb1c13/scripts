#!/usr/bin/env python

# Script to calculate the Amoeba-based electric field at a desired site
# Requirements: MDTraj, numpy

import mdtraj as md
import numpy as np
from subprocess import call

# Class for trajectory
class Anal_traj:
   """Class to define all the desired analyses for a
      Tinker trajectory read in using MDTraj."""   

   def __init__(self,arcname):
      """Reads in Tinker trajectory arcname."""
      self.traj = md.load_arc(arcname)

   def getcoords(self,atm1,atm2,atm3):
      """Gets coordinates of atom indices
         atm1 & atm2. Units in Angstrom."""
      self.coords1 = np.zeros((self.traj.n_frames,3))
      self.coords2 = np.zeros((self.traj.n_frames,3))
      for i in range(0,self.traj.n_frames):
         self.coords1[i] = self.traj[i].xyz[0][atm1] *10
         self.coords2[i] = self.traj[i].xyz[0][atm2] *10

   def midpoints(self,atm1,atm2):
      """Defines the Cartesian coords midpoint between
         two atom indices (atm1 & atm2) for every frame
         in a trajectory. Returns 3 x trajectory length 
         numpy array."""
      self.getcoords(atm1,atm2)
      self.midp = np.zeros((self.traj.n_frames,3))
      for i in range(0,self.traj.n_frames):
         self.midp[i] = (self.coords2[i] + self.coords1[i]) / 2

   def vectors(self,atm1,atm2):
      """Defines the interatomic vector, length and
         unit vector along that path for atom indices
         atm1, atm2. Returns three arrays:
          1) X/Y/Z vector components
          2) Interatomic vector length
          3) X/Y/Z unit vector components."""
      self.vec = np.zeros((self.traj.n_frames,3))
      self.leng = np.zeros(self.traj.n_frames)
      self.unitvec = np.zeros((self.traj.n_frames,3))
      for i in range(0,self.traj.n_frames):
         self.vec[i] = self.coords2[i] - self.coords1[i]
         self.leng[i] = np.linalg.norm(self.vec[i])
         self.unitvec[i] = self.vec[i] / self.leng[i]

   def cross(self,atm1,atm2,atm3):
      """Defines the vector cross product, length and
         unit vector along that path for atom indices
         atm1, atm2, atm3, centred at atm1.Returns three arrays:
          1) X/Y/Z vector components
          2) Interatomic vector length
          3) X/Y/Z unit vector components."""
      self.coords1 = np.zeros((self.traj.n_frames,3))
      self.coords2 = np.zeros((self.traj.n_frames,3))
      self.coords3 = np.zeros((self.traj.n_frames,3))
      for i in range(0,self.traj.n_frames):
         self.coords1[i] = self.traj[i].xyz[0][atm1] *10
         self.coords2[i] = self.traj[i].xyz[0][atm2] *10
         self.coords3[i] = self.traj[i].xyz[0][atm3] *10
      self.vec1 = np.zeros((self.traj.n_frames,3))
      self.vec2 = np.zeros((self.traj.n_frames,3))
      self.leng = np.zeros(self.traj.n_frames)
      self.unitvec = np.zeros((self.traj.n_frames,3))
      for i in range(0,self.traj.n_frames):
         self.vec1[i] = self.coords2[i] - self.coords1[i]
         self.vec2[i] = self.coords3[i] - self.coords1[i]
         cross = np.cross(self.vec1[i],self.vec2[i])
         self.leng[i] = np.linalg.norm(cross)
         self.unitvec[i] = cross / self.leng[i]
      
# Class for dipoles
class Anal_uind:
   """Class to define all the desired analyses for a
      Tinker induced dipoles file read in using MDTraj."""   

   def __init__(self,uindname):
      """Reads in Tinker uind file uindname."""
      self.traj = md.load_arc(uindname)
      # Convert units - mdtraj reads in as distance in nm!
      self.traj.xyz = self.traj.xyz *10


# Functions below here

def newxyz(prefix,trajobj):
   """Creates new xyzfiles with additional atom at 
      the midpoint between two other atoms. Currently
      hard coded to look for box coords in xyz and 
      print out new atom with type 257.
      
      Creates files as 'prefix_new.000'. No more than
      999 snapshots!"""

   for i in range(0,trajobj.traj.n_frames):
      xyzname=prefix + ".%03d" % int(i+1)
      modxyz=prefix + "_new.%03d" % int(i+1)
      xyzf = open(xyzname,'r')
      # Print out first two lines - number of atoms and box coords
      modf = open(modxyz,'w')
      line1=int(xyzf.readline())
      line1 +=1
      modf.write("%s \n" % line1)
      line2=xyzf.readline()
      modf.write("%s" % line2)
      # Now print out rest of the lines
      for line in range(line1-1):
         line=xyzf.readline()
         modf.write("%s" % line)
      # Insert the final line here
      modf.write("%6d  Du %12.6f %11.6f %11.6f   257 \n" % (line1,trajobj.midp[i][0],trajobj.midp[i][1],trajobj.midp[i][2]))
      modf.close()
      xyzf.close()
   
# Assign induced dipoles (hence field - unit polarisability) from analyze  
def analyze_dipl(prefix,trajobj):
   """Assigns induced dipole components for atom at
      midpoint using analyze. Note that the keyfile
      'test.key' must be prepared elsewhere, and only the
      dummy atom at the midpoint can be 'active'."""

   field = np.zeros((trajobj.traj.n_frames,3))
   for i in range(0,trajobj.traj.n_frames):
      modxyz=prefix + "_new.%03d" % int(i+1)
      call('analyze %s -k test.key m > analout.txt ' % modxyz,shell=True)
      f = open("analout.txt",'r')
      for line in f.readlines():
         if line.startswith(" Dipole X,Y,Z-Components"):
            s=line.split()
            field[i][0] = s[3]
            field[i][1] = s[4]
            field[i][2] = s[5]
      f.close()
   return field

# This one makes use of existing atoms - no dummy atms. Takes longer as prints all interactions.
# Can change analyze.f to improve this!
def analyze_dipl_detailed(prefix,trajobj,atm1,atm2):
   """Assigns induced dipole components for defined atoms
      using analyze. Note that the keyfile
      'test.key' must be prepared elsewhere, and should
      print out dipoles to 6dp."""
   # 3d array for field
   field = np.zeros((2,trajobj.traj.n_frames,3))
   for i in range(0,trajobj.traj.n_frames):
      xyz=prefix + ".%03d" % int(i+1)
      call('analyze %s -k test.key d > analout.txt ' % xyz,shell=True)
      call('grep -A%s "Induced Dipole Moments (Debyes) :" analout.txt > dipls.txt ' % (trajobj.traj.n_atoms+3),shell=True)
      dipls = np.loadtxt("dipls.txt",skiprows=4,usecols=(1,2,3))
      field[0][i][0] = dipls[atm1][0]
      field[0][i][1] = dipls[atm1][1]
      field[0][i][2] = dipls[atm1][2]
      field[1][i][0] = dipls[atm2][0]
      field[1][i][1] = dipls[atm2][1]
      field[1][i][2] = dipls[atm2][2]
   field[0] = field[0] / 1.3340 #Glycine C
   field[1] = field[1] / 1.0730 #Proline N
   return field


def writefield(fn,trajobj,field):
   """Writes out field projections in MV/cm 
      along unit vector at desired coordinate"""

   outfile = open(fn,'w')
   for i in range(0,trajobj.traj.n_frames):
      fieldproj = np.dot(field[i], trajobj.unitvec[i])
      fieldproj = fieldproj * 299.79
      outfile.write("%8.3f \n" % fieldproj)
   outfile.close()

# Projection along cross product of defined vectors
def writefield_cross(fn,trajobj,field):
   """Writes out field projections in MV/cm 
      along unit vector at desired coordinate"""

   outfile = open(fn,'w')
   for i in range(0,trajobj.traj.n_frames):
      fieldproj = np.dot(field[i], trajobj.unitvec[i])
      fieldproj = fieldproj * 299.79
      outfile.write("%8.3f \n" % fieldproj)
   outfile.close()

### MAIN BELOW HERE ###


# Read in tinker archive. INSERT FILENAME HERE.
# Currently needs to be an archive without box dims (due to bug)
arcname = "Run1_stripped.arc"

# Indices of atoms we want to measure the field between
atm1 = 2501
atm2 = 2502
atm3 = 2506

# Processing starts here
x = Anal_traj(arcname)
#x.getcoords(atm1,atm2)
#x.midpoints(atm1,atm2)
#x.vectors(atm1,atm2)
x.cross(atm1,atm2,atm3)

prefix="1M9C_run1"
#newxyz(prefix,x)
field = analyze_dipl_detailed(prefix,x,atm1,atm3)
writefield_cross("fields_atC_run1.txt",x,field[0])
writefield_cross("fields_atN_run1.txt",x,field[1])
#writefield("fields_atO.txt",x,field[1])

