#!/usr/bin/env python

# Script to print out a list of atoms within a defined radius
# of an atom/bond midpoint in a Tinker trajectory.
# DOES NOT TAKE PBC INTO ACCOUNT
# Requirements:  MDTraj, numpy, argparse

import mdtraj as md
import numpy as np
import sys,argparse


### Argparser ###
def parse():
   parser = argparse.ArgumentParser()
   parser.add_argument("-t","--trajectory",help="Trajectory file (Tinker arc, NO BOX INFO DUE TO BUG IN MDTRAJ)",type=str,required=True)
   parser.add_argument("-a","--atoms",help="Atom numbers for analysis (Max 2)",nargs='+',type=int,required=True)
   parser.add_argument("-r","--radius",help="Radius outside which atoms defined as inactive. Defaults to 9.0A",type=float,default=9.0)

   if len(sys.argv)==1:
      parser.print_help()
      sys.exit(1)

   args = parser.parse_args()
   return args


### Class for trajectory (taken from environment_fields) ###
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
         self.coords1[i] = self.traj[i].xyz[0][atmlst[0]-1]
         self.coords2[i] = self.traj[i].xyz[0][atmlst[1]-1]

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

### Functions below here ###

# Calc atoms in radius

def filter_atoms(trajectory,radius):
   """Returns list of atoms within a certain radius of
      a defined trajectory midpoint"""
   actives,inactives = [],[]
   for fr_idx,frame in enumerate(trajectory.traj):
      tmp_actives,tmp_inactives = [],[]
      for num,xyz in enumerate(frame.xyz[0],start=1):
         dist = np.linalg.norm(xyz-trajectory.midp[fr_idx])
         if dist > radius/10: # Radius in nm, as coords in nm
            tmp_actives.append(num)
         else:
            tmp_inactives.append(num)
      actives.append(tmp_actives)
      inactives.append(tmp_inactives)
   return actives,inactives   




# Main below here
if __name__ == "__main__":
   args = parse()
   archive = Anal_traj(args.trajectory)
   if len(args.atoms) > 1:
      archive.midpoints([args.atoms[0],args.atoms[1]])
   else:
      archive.midpoints([args.atoms[0],args.atoms[0]])
   actlst,inactlst=filter_atoms(archive,args.radius)
   for idx,frame in enumerate(inactlst):
      with open("restraints_{}.txt".format(idx+1),'w') as f:
         for atom in frame:
            f.write("restrain-position %s %10.6f %10.6f %10.6f 500.0\n" % 
                    (atom,archive.traj[0].xyz[idx][atom-1][0]*10,
                    archive.traj[0].xyz[idx][atom-1][1]*10,
                    archive.traj[0].xyz[idx][atom-1][2]*10))


