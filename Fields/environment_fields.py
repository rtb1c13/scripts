#!/usr/bin/env python

# Author Richard Bradshaw, R.T.Bradshaw@soton.ac.uk
# Script to calculate the Amoeba-based electric field at a desired site
# Requirements: MDTraj, numpy, Tinker analyze

#Usage: xxxxxxxxxxxxx


import mdtraj as md
import numpy as np
import sys,argparse
from subprocess import call


### Argparser ###
def parse():
   parser = argparse.ArgumentParser()
   parser.add_argument("-a","--atoms",help="Atom numbers for analysis (Max 4)",nargs='+',type=int,required=True)
   parser.add_argument("-ga","--gasatoms",help="Atom numbers for analysis (Max 4)",nargs='+',type=int)
   parser.add_argument("-t","--type",help="Type of analysis to run. 'dummy' will insert a dummy atom between atoms 1 & 2. 'standard' will evaluate fields at atoms 1 & 2. 'env' will take the difference between the fields in the solvated and gas (probe only) phase at atoms 1 & 2. Defaults to standard. ",choices=['dummy','standard','env'],default='standard')
   parser.add_argument("-p","--polarise",help="Polarisabilities for atoms 1 & 2. Defaults to carbonyl C and O.",nargs=2,type=float,default=[1.3340,0.8370])
   parser.add_argument("-st","--solvtraj",help="Solvated trajectory file (Tinker arc, NO BOX INFO DUE TO BUG IN MDTRAJ)",type=str,required=True)
   parser.add_argument("-gt","--gastraj",help="Gas phase trajectory file (Tinker arc, NO BOX INFO DUE TO BUG IN MDTRAJ)",type=str)
   parser.add_argument("-sp","--solvprefix",help="Prefix for solvated coordinate files (Tinker xyz, with box info)",type=str)
   parser.add_argument("-gp","--gasprefix",help="Prefix for gas-phase coordinate files (Tinker xyz, no box info)",type=str)

   if len(sys.argv)==1:
      parser.print_help()
      sys.exit(1)

   args = parser.parse_args()
   return args


### Class for atoms ###
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

### Class for trajectory ###
class Anal_traj:
   """Class to define all the desired analyses for a
      Tinker trajectory read in using MDTraj."""   

   def __init__(self,arcname):
      """Reads in Tinker trajectory arcname."""
      self.traj = md.load_arc(arcname)

   def getcoords(self,atmlst):
      """Gets coordinates of atom indices
         atm1 & atm2 given in a list. Units in Angstrom."""
      self.coords1 = np.zeros((self.traj.n_frames,3))
      self.coords2 = np.zeros((self.traj.n_frames,3))
      for i in range(0,self.traj.n_frames):
         self.coords1[i] = self.traj[i].xyz[0][atmlst[0].idx] *10
         self.coords2[i] = self.traj[i].xyz[0][atmlst[1].idx] *10

   def midpoints(self,atmlst):
      """Defines the Cartesian coords midpoint between
         two atom indices (atm1 & atm2) for every frame
         in a trajectory. Returns 3 x trajectory length 
         numpy array."""
      self.getcoords(atmlst)
      self.midp = np.zeros((self.traj.n_frames,3))
      for i in range(0,self.traj.n_frames):
         self.midp[i] = (self.coords2[i] + self.coords1[i]) / 2

   def vectors(self):
      """Defines the interatomic vector, length and
         unit vector along that path, assuming coordinates
         for atoms are already defined. Returns three arrays:
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
      call('analyze %s -k test.key m &> analout.txt ' % modxyz,shell=True)
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
def analyze_dipl_detailed(prefix,trajobj,atmlst,gas):
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
      if gas is True:
         call('analyze %s -k testgas.key d &> analout.txt ' % xyz,shell=True)
      else:
         call('analyze %s -k test.key d &> analout.txt ' % xyz,shell=True)
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
def writefield_components(fn,trajobj,field):
   """Writes out field xyz components in MV/cm 
      at desired coordinate, along with unit vector
      for projection"""

   currfield = field * 299.79
   output = np.concatenate((currfield,trajobj.unitvec),axis=1)
   np.savetxt(fn,output,fmt="%8.3f")

def dummy_field():
   """Analyses field between two atoms by inserting
      a dummy atom at midpoint between them. Field is
      projected along internuclear vector"""
   print "Reading in trajectory %s as the complete system. Ignoring any other flags" % args.solvtraj
   arcname = args.solvtraj
   arc = Anal_traj(arcname)
   print "Calculating field between atoms %d and %d. Ignoring any other atoms defined" % (args.atoms[0],args.atoms[1])
   atmlst=[]
   for i in range(0,2):
      j = Atom(args.atoms[i])
      j.polarise(args.atoms[i])
      atmlst.append(j)
   # Get coordinates, midpoints & interatomic vector
   arc.midpoints(atmlst)
   arc.vectors()
   # Get field
   newxyz(args.solvprefix,arc)
   field = analyze_dipl(args.solvprefix,arc)
   writefield("fields_at_dummy.txt",arc,field)

def standard_field():
   """Analyses field at two atoms No gas-phase correction is made for
      the self field of the probe. Field is projected along
      internuclear vector or cross product of two related bonds"""
   print "Reading in trajectory %s as the complete system. Ignoring any other flags" % args.solvtraj
   arcname = args.solvtraj
   arc = Anal_traj(arcname)
   print "Calculating field at atoms %d and %d." % (args.atoms[0],args.atoms[1])
   atmlst=[]
   for i in args.atoms:
      j = Atom(i)
      j.polarise(i)
      atmlst.append(j)
   # Decide on where to project according to length of atom list
   if len(atmlst) == 2:
      print "Projecting along internuclear vector %d to %d" % (args.atoms[0],args.atoms[1])
      arc.getcoords(atmlst)
      arc.vectors()
   elif len(atmlst) == 3:
      print "Projecting along cross product %d to %d x %d to %d" % (args.atoms[0],args.atoms[1],args.atoms[0],args.atoms[2])
      arc.cross(atmlst)
   elif len(atmlst) == 4:
      print "Projecting along cross product %d to %d x %d to %d" % (args.atoms[0],args.atoms[2],args.atoms[0],args.atoms[3])
      arc.cross(atmlst)
   # Get field
   gas = False
   field = analyze_dipl_detailed(args.solvprefix,arc,atmlst,gas)
   writefield(("field_at_atom_%i.txt" % args.atoms[0]),arc,field[0])
   writefield(("field_at_atom_%i.txt" % args.atoms[1]),arc,field[1])

def env_field():
   """Analyses field at two atoms Correction is made for
      the self field of the probe. Field is projected along
      internuclear vector or cross product of two related bonds"""
   print "Reading in trajectory %s as the solvated system and %s as the gas system" % (args.solvtraj, args.gastraj)
   arcname = args.solvtraj
   arcname_gas = args.gastraj
   arc = Anal_traj(arcname)
   gasarc = Anal_traj(arcname_gas)
   print "Calculating field at atoms %d and %d." % (args.atoms[0],args.atoms[1])
   atmlst=[]
   gatmlst=[]
   for i in args.atoms:
      j = Atom(i)
      j.polarise(i)
      atmlst.append(j)
   for i in args.gasatoms:
      j = Atom(i)
      j.polarisegas(i)
      gatmlst.append(j)
   # Decide on where to project according to length of atom list
   if len(atmlst) == 2:
      print "Projecting along internuclear vector %d to %d" % (args.atoms[0],args.atoms[1])
      arc.getcoords(atmlst)
      arc.vectors()
      gasarc.getcoords(gatmlst)
      gasarc.vectors()
   elif len(atmlst) == 3:
      print "Projecting along cross product %d to %d x %d to %d" % (args.atoms[0],args.atoms[1],args.atoms[0],args.atoms[2])
      arc.cross(atmlst)
      gasarc.cross(gatmlst)
   elif len(atmlst) == 4:
      print "Projecting along cross product %d to %d x %d to %d" % (args.atoms[0],args.atoms[2],args.atoms[0],args.atoms[3])
      arc.cross(atmlst)
      gasarc.cross(gatmlst)
   # Get field
   gas = False
   solfield = analyze_dipl_detailed(args.solvprefix,arc,atmlst,gas)
   gas = True
   gasfield = analyze_dipl_detailed(args.gasprefix,gasarc,gatmlst,gas)
   field = solfield - gasfield
   writefield(("field_at_atom_%i.txt" % args.atoms[0]),arc,field[0])
   writefield(("field_at_atom_%i.txt" % args.atoms[1]),arc,field[1])
   writefield_components(("field_components_at_atom_%i.txt" % args.atoms[0]),arc,field[0])
   writefield_components(("field_components_at_atom_%i.txt" % args.atoms[1]),arc,field[1])

### MAIN BELOW HERE ###
def main():

   global args
   global gas
   args = parse()
   if args.type == "dummy":
      dummy_field()
   elif args.type == "standard":
      standard_field()
   elif args.type == "env":
      env_field()

############################################
main()
