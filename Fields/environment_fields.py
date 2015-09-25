#!/usr/bin/env python

# Author Richard Bradshaw, R.T.Bradshaw@soton.ac.uk
# Script to calculate the Amoeba-based electric field at a desired site
# Requirements: MDTraj, numpy, Tinker analyze

#Usage: xxxxxxxxxxxxx


import mdtraj as md
import numpy as np
import sys,argparse
from subprocess import call


#Argparser
def parse():
   parser = argparse.ArgumentParser()
   parser.add_argument("-a","--atoms",help="Atom numbers for analysis (Max 4)",nargs='+',type=int)
   parser.add_argument("-ga","--gasatoms",help="Atom numbers for analysis (Max 4)",nargs='+',type=int)
   parser.add_argument("-d","--dummy",help="Calculates fields by inserting a dummy atom instead of using existing atom. Use with caution!",action="store_true")
   parser.add_argument("-p","--polarise",help="Polarisabilities for atoms 1 & 2. Defaults to carbonyl C and O.",nargs=2,type=float,default=['1.3340','0.8370'])
   parser.add_argument("-st","--solvtraj",help="Solvated trajectory file (Tinker arc, NO BOX INFO DUE TO BUG IN MDTRAJ)",type=str)
   parser.add_argument("-gt","--gastraj",help="Gas phase trajectory file (Tinker arc, NO BOX INFO DUE TO BUG IN MDTRAJ)",type=str)
   parser.add_argument("-sp","--solvprefix",help="Prefix for solvated coordinate files (Tinker xyz, with box info)",type=str)
   parser.add_argument("-gp","--gasprefix",help="Prefix for gas-phase coordinate files (Tinker xyz, no box info)",type=str)

   if len(sys.argv)==1:
      parser.print_help()
      sys.exit(1)

   args = parser.parse_args()
   return args

#Class for atoms
class Atom:
   """Class to define properties for atom numbers
      read in from command line arguments"""

   def __init__(self,num):
      """Defines index for atom num"""
      self.idx = num-1

   def polarise(self,num):
      """Defines polarisability for atom num"""
      listidx = args.atoms.index(num)
      if listidx > 1:
         return
      else:
         self.alpha = args.polarise[listidx]
   def polarisegas(self,num):
      """Defines polarisability for atom num"""
      listidx = args.gasatoms.index(num)
      if listidx > 1:
         return
      else:
         self.alpha = args.polarise[listidx]

# Class for trajectory
class Anal_traj:
   """Class to define all the desired analyses for a
      Tinker trajectory read in using MDTraj."""   

   def __init__(self,arcname):
      """Reads in Tinker trajectory arcname."""
      self.traj = md.load_arc(arcname)

   def getcoords(self,atm1,atm2):
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

   def cross(self,atmlst):
      """Defines the vector cross product, length and
         unit vector along that path for atom indices
         atm1, atm2, atm3, centred at atm1.Returns three arrays:
          1) X/Y/Z vector components
          2) Interatomic vector length
          3) X/Y/Z unit vector components."""
      if len(atmlst) > 3:
         atm1 = atmlst[0].idx
         atm2 = atmlst[2].idx
         atm3 = atmlst[3].idx
      else:
         atm1 = atmlst[0].idx
         atm2 = atmlst[1].idx
         atm3 = atmlst[2].idx
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
def analyze_dipl_detailed(prefix,trajobj,atmlst):
   """Assigns induced dipole components for defined atoms
      using analyze. Note that the keyfile
      'test.key' must be prepared elsewhere, and should
      print out dipoles to 6dp."""
   atm1=atmlst[0].idx
   atm2=atmlst[1].idx
   alpha1=atmlst[0].alpha
   alpha2=atmlst[1].alpha
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
   field[0] = field[0] / alpha1 
   field[1] = field[1] / alpha2 
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
#def main()

# Read in tinker archive. INSERT FILENAME HERE.
# Currently needs to be an archive without box dims (due to bug)
args = parse()
arcname = args.solvtraj
arcname_gas = args.gastraj
atmlst=[]
for i in args.atoms:
   j = Atom(i)
   j.polarise(i)
   atmlst.append(j)
gatmlst=[]
for i in args.gasatoms:
   j = Atom(i)
   j.polarisegas(i)
   gatmlst.append(j)

# Indices of atoms we want to measure the field between, and calculate z axis with
#atm1 = 2506 # Proline N
#atm2 = 2507 # Proline CA
#atm3 = 2513 # Proline CD
#atm4 = 2501 # Glycine C

#atm5 = 36 # Pro N
#atm6 = 37 # Pro CA
#atm7 = 43 # Pro CD
#atm8 = 31 # Gly C

# Processing starts here
arc = Anal_traj(arcname)
gasarc = Anal_traj(arcname_gas)
#arc.getcoords(atm1,atm2)
#gasarc.getcoords(atm1,atm2)
#x.midpoints(atm1,atm2)
#arc.vectors(atm1,atm2)
#gasarc.vectors(atm1,atm2)

# This time we need to create a new Z axis using the proline N/CA/CD
arc.cross(atmlst)
gasarc.cross(gatmlst)

# These should be the same right?
print arc.unitvec - gasarc.unitvec


# Here we have to subtract the self field and average field str. across C=O
prefix=args.solvprefix
prefix2=args.gasprefix
#newxyz(prefix,x)
solfield = analyze_dipl_detailed(prefix,arc,atmlst)
gasfield = analyze_dipl_detailed(prefix2,gasarc,gatmlst)
field = solfield - gasfield

writefield_cross(("fields_at_atom_%d.txt" % args.atoms[0]),arc,field[0])
writefield_cross(("fields_at_atom_%d.txt" % args.atoms[1]),arc,field[1])
#writefield("fields_atO.txt",x,field[1])

