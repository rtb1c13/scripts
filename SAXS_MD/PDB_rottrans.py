#!/usr/bin/env python

# Script to make transform coordinates of pdb. Makes use of biopython (only python 2.6 on this machine?!)
# and python array module.

from Bio.PDB import *
import array as arr

#Name of PDB files
infile = "../6G08_fab.pdb"
outfile = "../6G08_fab_trans.pdb"

#Rotation matrix & translation vector
rot = rotmat(Vector(1,0,0),Vector(1,0,0)) # i.e. identity matrix
trans = arr.array('f',(10,10,10)) # +10 A in all directions

#Biopython below here

parser = PDBParser()
oldstruct = parser.get_structure('6g08',infile)
newstruct = oldstruct.copy
for i in newstruct.get_atoms():
   i.transform(rot,trans)

io = PDBIO()
io.set_structure(newstruct)
io.save(outfile)
