#!/usr/bin/env python

# Author Richard Bradshaw, R.T.Bradshaw@soton.ac.uk
# Script to calculate the Amber-based (fixed-charge)
# electric field at a desired site
# Requirements: AmberTools15 python APIs for sander & parmed, numpy

#Usage: xxxxxxxxxxxxx

import sander
import ParmedTools as parmed
import numpy as np
from chemistry import amber, unit as u


### Class for atoms ###
class Atom:
   """Class to define properties for atom numbers
      read in from command line arguments"""

   def __init__(self,num):
      """Defines index for atom num"""
      self.idx = num-1

   def charge(self,prm):
      """Defines partial charge on atom num from prmtop"""
#      listidx = args.atoms.index(num)
#      if listidx > 1:
#         return
#      else:
      self.chg = prm.orig.atoms[self.idx].charge

### Class for Parmfile ###
class Prmfile:
   """Class to define properties for prmtop file
      read in from command line arguments"""

   def __init__(self,parmname):
      """Reads in Amber prmtop parmname."""
      self.orig = amber.AmberParm(parmname)
      self.discharged = amber.AmberParm(parmname)

   def discharge(self,atmmask):
      """Sets all charges of atoms NOT in provided
         atmmask to 0 in a provided parmtop object.
         Returns new, discharged, prmtop"""

      invmask = "!" + atmmask
#      print invmask
#      newparm = self.orig.copy(amber.AmberParm)
      act = parmed.change(self.discharged, 'charge', invmask, 0.000)
      act.execute()
#      print newparm
#      self.discharged = newparm

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
      if len(atmlst) > 3:
         atm1 = atmlst[0].idx
         atm2 = atmlst[2].idx
         atm3 = atmlst[3].idx
      else:
         atm1 = atmlst[0].idx
         atm2 = atmlst[1].idx
         atm3 = atmlst[2].idx
      self.coords1 = np.zeros((self.traj.frame,3))
      self.coords2 = np.zeros((self.traj.frame,3))
      self.coords3 = np.zeros((self.traj.frame,3))
      for i in range(0,self.traj.frame):
         self.coords1[i] = self.traj.coordinates(i)[(atm1*3):(atm1*3)+3]
         self.coords2[i] = self.traj.coordinates(i)[(atm2*3):(atm2*3)+3]
         self.coords3[i] = self.traj.coordinates(i)[(atm3*3):(atm3*3)+3]
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

   def sanderforce(self,parm,boxflag=False,pmeflag=False):
      """Calculates energy and forces for a provided
         Amber parm object and set of coordinates using the Sander
         API. Also takes box dimensions and a flag for use of PME.
         Currently uses standard simulation options for gas phase or pme
         based on pmeflag.
         Returns a sander energy object and a 3 x n array of forces with units"""

      self.forces=[]
      self.energies=[]
      if pmeflag is True:
         inp = sander.pme_input()
      else:
         inp = sander.gas_input()
      for i in range(0,self.traj.frame):
         coord = self.traj.coordinates(i)
         if boxflag is True:
            box = self.traj.box(i)
#            print box
         else:
            box = None
         with sander.setup(parm,coord,box=box,mm_options=inp):
            ene,frc = sander.energy_forces()
#         print i
         frc = np.asarray(frc)
         frc = np.reshape(frc,((len(frc)/3.),3))
         # Add units for forces: kcal mol-1 A-1
         frc = frc * u.kilocalorie / (u.mole * u.angstroms)
         self.forces.append(frc)
         self.energies.append(ene)

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
         vec = frc[n][atmlst[0].idx][i]
	 fields[0][n][i] = vec.value_in_unit(u.megajoules/(u.centimeter*u.mole))
   for n in range(0,len(frc)):
      for i in range(0,3):
         vec = frc[n][atmlst[1].idx][i]
	 fields[1][n][i] = vec.value_in_unit(u.megajoules/(u.centimeter*u.mole))
   # Divide through by Faraday constant to get fields in MV/cm
   # AND divide through by partial charge
   fields[0] = fields[0]/(atmlst[0].chg * 96485.3329)
   fields[1] = fields[1]/(atmlst[1].chg * 96485.3329)
   return fields





#import mdtraj as md
#import numpy as np
#import sys,argparse
#from subprocess import call


### Argparser ###
#def parse():
#   parser = argparse.ArgumentParser()
#   parser.add_argument("-a","--atoms",help="Atom numbers for analysis (Max 4)",nargs='+',type=int,required=True)
#   parser.add_argument("-ga","--gasatoms",help="Atom numbers for analysis (Max 4)",nargs='+',type=int)
#   parser.add_argument("-t","--type",help="Type of analysis to run. 'dummy' will insert a dummy atom between atoms 1 & 2. 'standard' will evaluate fields at atoms 1 & 2. 'env' will take the difference between the fields in the solvated and gas (probe only) phase at atoms 1 & 2. Defaults to standard. ",choices=['dummy','standard','env'],default='standard')
#   parser.add_argument("-p","--polarise",help="Polarisabilities for atoms 1 & 2. Defaults to carbonyl C and O.",nargs=2,type=float,default=['1.3340','0.8370'])
#   parser.add_argument("-st","--solvtraj",help="Solvated trajectory file (Tinker arc, NO BOX INFO DUE TO BUG IN MDTRAJ)",type=str,required=True)
#   parser.add_argument("-gt","--gastraj",help="Gas phase trajectory file (Tinker arc, NO BOX INFO DUE TO BUG IN MDTRAJ)",type=str)
#   parser.add_argument("-sp","--solvprefix",help="Prefix for solvated coordinate files (Tinker xyz, with box info)",type=str)
#   parser.add_argument("-gp","--gasprefix",help="Prefix for gas-phase coordinate files (Tinker xyz, no box info)",type=str)
#
#   if len(sys.argv)==1:
#      parser.print_help()
#      sys.exit(1)
#
#   args = parser.parse_args()
#   return args


# Projection along cross product of defined vectors
def writefield_cross(fn,trajobj,field):
   """Writes out field projections in MV/cm 
      along unit vector at desired coordinate"""

   outfile = open(fn,'w')
   for i in range(0,trajobj.traj.frame):
      fieldproj = np.dot(field[i], trajobj.unitvec[i])
      outfile.write("%8.3f \n" % fieldproj)
   outfile.close()

#def env_field():
#   """Analyses field at two atoms Correction is made for
#      the self field of the probe. Field is projected along
#      internuclear vector or cross product of two related bonds"""
#   print "Reading in trajectory %s as the solvated system and %s as the gas system" % (args.solvtraj, args.gastraj)
#   arcname = args.solvtraj
#   arcname_gas = args.gastraj
#   arc = Anal_traj(arcname)
#   gasarc = Anal_traj(arcname_gas)
#   print "Calculating field at atoms %d and %d." % (args.atoms[0],args.atoms[1])
#   atmlst=[]
#   gatmlst=[]
#   for i in args.atoms:
#      j = Atom(i)
#      j.polarise(i)
#      atmlst.append(j)
#   for i in args.gasatoms:
#      j = Atom(i)
#      j.polarisegas(i)
#      gatmlst.append(j)
#   # Decide on where to project according to length of atom list
#   if len(atmlst) == 2:
#      print "Projecting along internuclear vector %d to %d" % (args.atoms[0],args.atoms[1])
#      arc.getcoords(atmlst)
#      arc.vectors()
#      gasarc.getcoords(gatmlst)
#      gasarc.vectors()
#   elif len(atmlst) == 3:
#      print "Projecting along cross product %d to %d x %d to %d" % (args.atoms[0],args.atoms[1],args.atoms[0],args.atoms[2])
#      arc.cross(atmlst)
#      gasarc.cross(gatmlst)
#   elif len(atmlst) == 4:
#      print "Projecting along cross product %d to %d x %d to %d" % (args.atoms[0],args.atoms[2],args.atoms[0],args.atoms[3])
#      arc.cross(atmlst)
#      gasarc.cross(gatmlst)
#   # Get field
#   gas = False
#   solfield = analyze_dipl_detailed(args.solvprefix,arc,atmlst,gas)
#   gas = True
#   gasfield = analyze_dipl_detailed(args.gasprefix,gasarc,gatmlst,gas)
#   field = solfield - gasfield
#   writefield(("field_at_atom_%i.txt" % args.atoms[0]),arc,field[0])
#   writefield(("field_at_atom_%i.txt" % args.atoms[1]),arc,field[1])

### MAIN BELOW HERE ###
#def main():
#
#   global args
#   global gas
#   args = parse()
#   if args.type == "dummy":
#      dummy_field()
#   elif args.type == "standard":
#      standard_field()
#   elif args.type == "env":
#      env_field()


prmfile = "../1M9C_solv_ff14SB.prmtop"
trajfile = "../1M9C_run1.mdcrd"

parameters = Prmfile(prmfile)
parameters.discharge(":164-169")
orig_trajectory = Anal_traj(trajfile)
new_trajectory = Anal_traj(trajfile)
atoms = [2507,2505,2517,2508]
for i in range(0,len(atoms)):
  atoms[i] = Atom(atoms[i])
  atoms[i].charge(parameters)
# Work out original forces
orig_trajectory.sanderforce(parameters.orig,boxflag=True,pmeflag=True)
new_trajectory.sanderforce(parameters.discharged,boxflag=True,pmeflag=True)
elefrc = subtract_frcs(orig_trajectory.forces,new_trajectory.forces)
#print elefrc[0][atoms[0].idx]
fields = frc2field(elefrc,atoms)
orig_trajectory.cross(atoms)
writefield_cross("atom_2507.out",orig_trajectory,fields[0])
writefield_cross("atom_2505.out",orig_trajectory,fields[1])

############################################
#main()
