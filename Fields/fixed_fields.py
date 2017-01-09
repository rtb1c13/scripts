#!/usr/bin/env python

# Author Richard Bradshaw, R.T.Bradshaw@soton.ac.uk
# Script to calculate the Amber-based (fixed-charge)
# electric field at a desired site
# Requirements: AmberTools15 python APIs for sander & parmed, numpy

# ***NB*** currently set up for gas phase complex (pmeflag=False)
# Change this if you want to run e.g. arotein in PBC as the complex


#Usage: xxxxxxxxxxxxx

import sander
import sys,argparse
import numpy as np
from parmed import amber, unit as u, tools as ptools

### Argparser ###
def parse():
   parser = argparse.ArgumentParser()
   parser.add_argument("-p","--parm",help="Path to Amber parameter file for analysis",type=str,required=True)
   parser.add_argument("-c","--coords",help="Path to NetCDF trajectory file for analysis",type=str,required=True)
   parser.add_argument("-a","--atoms",help="Atom numbers for analysis (Max 4)",nargs='+',type=int,required=True)
   parser.add_argument("-m","--mask",help="Mask (in Amber atom mask format) for atoms not to discharge whan calculating electrostatic forces (i.e. the 'ligand')",type=str,required=True)

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

   def charge(self,prm):
      """Defines partial charge on atom num from prmtop"""
      listidx = args.atoms.index(self.idx+1)
      if listidx > 1:
         return
      else:
         self.chg = prm.orig.atoms[self.idx].charge

### Class for Parmfile ###
class Prmfile:
   """Class to define properties for prmtop file
      read in from command line arguments"""

   def __init__(self,parmname):
      """Reads in Amber prmtop parmname."""

      self.orig = amber.AmberParm(parmname)
      self.discharged = amber.AmberParm(parmname) # Bug in copy method = read in twice instead

   def discharge(self,atmmask):
      """Sets all charges of atoms NOT in provided
         atmmask to 0 in a provided parmtop object.
         Returns new, discharged, prmtop"""

      if atmmask=="":
         invmask = "*"
      else:
         invmask = "!" + atmmask
#      print "invmask
      act = ptools.change(self.discharged, 'charge', invmask, 0.000)
      act.execute()

###Class for trajectories###
class Anal_traj:
   """Class to define all the desired analyses for a
      NetCDF trajectory read in using ParmEd."""   

   def __init__(self,trajname):
      """Reads in NetCDF trajectory arcname."""
      self.traj = amber.NetCDFTraj.open_old(trajname)

   def cross(self,atmlst):
      """Defines the vector cross product, length and
         unit vector along that path for atom indices
         atm1, atm2, atm3, centred at atm1.Returns three arrays:
          1) X/Y/Z vector components
          2) Interatomic vector length
          3) X/Y/Z unit vector components."""

      if len(atmlst) > 3: # 4 atoms = cross product of 1,3,4
         atm1 = atmlst[0].idx
         atm2 = atmlst[2].idx
         atm3 = atmlst[3].idx
      else:               # 3 atoms = cross procuct of 1,2,3
         atm1 = atmlst[0].idx
         atm2 = atmlst[1].idx
         atm3 = atmlst[2].idx
      self.coords1 = np.zeros((self.traj.frame,3))
      self.coords2 = np.zeros((self.traj.frame,3))
      self.coords3 = np.zeros((self.traj.frame,3))
      for i in range(0,self.traj.frame):
         self.coords1[i] = self.traj.coordinates[i][atm1]
         self.coords2[i] = self.traj.coordinates[i][atm2]
         self.coords3[i] = self.traj.coordinates[i][atm3]
      self.vec1 = np.zeros((self.traj.frame,3))
      self.vec2 = np.zeros((self.traj.frame,3))
      self.leng = np.zeros(self.traj.frame)
      self.unitvec = np.zeros((self.traj.frame,3))
      for i in range(0,self.traj.frame):
         self.vec1[i] = self.coords2[i] - self.coords1[i]
         self.vec2[i] = self.coords3[i] - self.coords1[i]
         cross = np.cross(self.vec1[i],self.vec2[i])
         self.leng[i] = np.linalg.norm(cross)
         self.unitvec[i] = cross / self.leng[i]

   def vectors(self,atmlst):
      """Defines the interatomic vector, length and
         unit vector along that path between atoms
         1 & 2. Returns three arrays:
          1) X/Y/Z vector components
          2) Interatomic vector length
          3) X/Y/Z unit vector components."""
      self.vec = np.zeros((self.traj.frame,3))
      self.leng = np.zeros(self.traj.frame)
      self.unitvec = np.zeros((self.traj.frame,3))
      atm1 = atmlst[0].idx
      atm2 = atmlst[1].idx
      self.coords1 = np.zeros((self.traj.frame,3))
      self.coords2 = np.zeros((self.traj.frame,3))
      for i in range(0,self.traj.frame):
         self.coords1[i] = self.traj.coordinates[i][atm1]
         self.coords2[i] = self.traj.coordinates[i][atm2]
      for i in range(0,self.traj.frame):
         self.vec[i] = self.coords2[i] - self.coords1[i]
         self.leng[i] = np.linalg.norm(self.vec[i])
         self.unitvec[i] = self.vec[i] / self.leng[i]


   def sanderforce(self,parmstr,atmlst,boxflag=False,pmeflag=False):
      """Calculates energy and forces for a provided
         Amber parm string and set of coordinates using the Sander
         API. Also takes box dimensions and a flag for use of PME.
         Currently uses standard simulation options for gas phase or pme
         based on pmeflag.
         Only returns forces on first 2 atoms in atmlst
         Returns a sander energy object and a 3 x 2 x traj length array
         of forces with units"""

      self.forces=[]
      self.energies=[]
#      sander.APPLY_UNITS = True # More bugs/inconsistencies in pysander in amber15, may be fixed later
      indices = [a.idx for a in atmlst[:2]]
      if pmeflag is True:
#        These default options will be suitable in most cases
#        PME,cut=8.0,ntb=1,ntf=1,ntc=1
#        Should be adapted if shake is required (not important here for ele)
         inp = sander.pme_input()
      else:
         inp = sander.gas_input()
      for i in range(0,self.traj.frame):
         coord = self.traj.coordinates[i]
         if boxflag is True:
            box = self.traj.box[i]
#           print box
         else:
            box = None
         with sander.setup(parmstr,coord,box=box,mm_options=inp): #Pass a string, else a temp parmfile is created that fills up tmp directory!
            ene,frc = sander.energy_forces(as_numpy=False) # pysander __init__ bug/inconsistency, can't use as_numpy!
         frc = np.asarray(frc)
         frc = np.reshape(frc,((len(frc)/3.),3))
         frcslice = frc[indices]
         # Add units for forces: kcal mol-1 A-1
         frcslice = frcslice * u.kilocalorie / (u.mole * u.angstroms)
         self.forces.append(frcslice)
         self.energies.append(ene)
#      sander.APPLY_UNITS = False

#Functions below here

def subtract_frcs(frc1,frc2):
   """Subtracts one set of forces from another. E.g.
      totalfrc - (totalbutele + pepele) = envelefrc
      The envelefrc is the force due to the electrostatics
      of the environment (i.e. everything except peptide)
      Requires 2 lists of arrays. Returns 1 list of arrays."""

   if len(frc1) != len(frc2):
      print "The charged and discharged force sets are different lengths! Check your inputs"
      quit()
   newfrc = []
   for n in range(0,len(frc1)):
      newfrc.append(frc1[n] - frc2[n])
   return newfrc

def frc2field(frc,atmlst):
   """Turns a set of electrostatic forces on the desired
      peptide atoms into field vectors at those points.
      Requires atom list and the array of forces with
      units.
      Returns only force vectors at desired atoms (2)"""

   #Extract desired forces
   fields = np.zeros((2,len(frc),3))
   for n in range(0,len(frc)):
      for i in range(0,3):
#         vec = frc[n][atmlst[0].idx][i]
         vec = frc[n][0][i]
	 fields[0][n][i] = vec.value_in_unit(u.megajoules/(u.centimeter*u.mole))
   for n in range(0,len(frc)):
      for i in range(0,3):
#         vec = frc[n][atmlst[1].idx][i]
         vec = frc[n][1][i]
	 fields[1][n][i] = vec.value_in_unit(u.megajoules/(u.centimeter*u.mole))
   # Divide through by Faraday constant to get fields in MV/cm
   # AND divide through by partial charge
   fields[0] = fields[0]/(atmlst[0].chg * 96485.3329)
   fields[1] = fields[1]/(atmlst[1].chg * 96485.3329)
   return fields

def writefield(fn,trajobj,field):
   """Writes out field projections in MV/cm 
      along unit vector at desired coordinate"""

   outfile = open(fn,'w')
   for i in range(0,trajobj.traj.frame):
      fieldproj = np.dot(field[i], trajobj.unitvec[i])
      outfile.write("%8.3f \n" % fieldproj)
   outfile.close()

def writefield_components(fn,field):
   """Writes out components of field in MV/cm 
      at desired coordinate"""

   np.savetxt(fn,field,fmt='%8.3f')

def debug_energies(fn,traj):
   """Prints out energies calculated internally
      using pysander API for comparison with those
      from sander standalone.

      Note: These will NOT agree with energies in 
      simulation output files - output file energies
      are a half-timestep removed from coordinates!
      Rerun sander with imin=5 to obtain correct energies
      for frames in coordinate file."""

   outfile = open(fn,'w')
   keylist = ["Tot","Bond","Angle","Dihedral","Elec","Elec 1-4","VDW","VDW 1-4"]
   row_format = "{:>12}" * len(keylist)
   outfile.write(row_format.format(*keylist)+"\n")
   row_format = "{:>12.4f}" * len(keylist)
   for e in traj.energies:
      enelist = [e.tot,e.bond,e.angle,e.dihedral,e.elec,e.elec_14,e.vdw,e.vdw_14]
      outfile.write(row_format.format(*enelist)+"\n") 
   outfile.close()

### MAIN BELOW HERE ###
def main():

   global args
   args = parse()
   prmfile = args.parm
   trajfile = args.coords

   # Assign params and trajectories
   parameters = Prmfile(prmfile)
   parameters.discharge(args.mask)
   parameters.discharged.write_parm("discharged.parm7")
   orig_trajectory = Anal_traj(trajfile)
   new_trajectory = Anal_traj(trajfile)
   atoms = [Atom(i) for i in args.atoms]
   [i.charge(parameters) for i in atoms]
   # Work out original forces
   orig_trajectory.sanderforce(args.parm,atoms,boxflag=False,pmeflag=False)
   new_trajectory.sanderforce("discharged.parm7",atoms,boxflag=False,pmeflag=False)
   debug_energies("Energies.txt",orig_trajectory)
   elefrc = subtract_frcs(orig_trajectory.forces,new_trajectory.forces)
   # Fields and projection
   fields = frc2field(elefrc,atoms)
   if len(atoms) > 2:
      orig_trajectory.cross(atoms)
   else:
      orig_trajectory.vectors(atoms)
   writefield("field_at_%i.txt" % args.atoms[0],orig_trajectory,fields[0])
   writefield("field_at_%i.txt" % args.atoms[1],orig_trajectory,fields[1])
   writefield_components("field_components_at_%i.txt" % args.atoms[0],fields[0])
   writefield_components("field_components_at_%i.txt" % args.atoms[1],fields[1])

############################################
main()
