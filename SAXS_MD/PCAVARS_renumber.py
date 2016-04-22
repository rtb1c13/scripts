#!/usr/bin/env python

# Script to renumber atoms in a gromacs eigenvector PDB. Designed
# to make sure atoms subjected to PCAVARS/PCARMSD in plumed
# are correctly identified

# Requirements: Biopython
# Inputs: Average PDB structure (with correct atom names)
#         Eigenvector PDB (all models will be renumbered)
# Outputs: Renumbered Eigenvector PDB


from Bio import PDB
import itertools

# Inputs
avename = '6G08_average.pdb'
eigtraj = '6G08_eigenvec.pdb'
outname = 'test.pdb'

par = PDB.PDBParser()
avepdb = par.get_structure('ave',avename)
eigpdb = par.get_structure('eig',eigtraj)


# Extract atoms from average
ave_atoms = PDB.Selection.unfold_entities(avepdb,target_level='A') 

# Iterate over atoms in eigpdb
for model in eigpdb:
   for chain in model:
      eigatms = chain.get_atoms() # atoms generator
      i = itertools.izip(ave_atoms,eigatms)
      for old,new in i:
         new.set_serial_number(old.get_serial_number) # Replace new with old
         print new


io = PDB.PDBIO()
io.set_structure(eigpdb)
io.save(outname)
