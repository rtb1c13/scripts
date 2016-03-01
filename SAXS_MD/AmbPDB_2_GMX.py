#!/usr/bin/env python

# Script to convert Amber protein hydrogens (created by ambpdb) into Gromacs
# names. For some reason, eg. HB2 + HB3 are expected to be HB1 and HB2 in the
# Gromacs aminoacids.rtp parameter files.

# Requirements: Biopython

# Desired PDB file should be the first argument. Output is renamed PDB with
# 'GMX_' prepended.

import sys
from Bio.PDB import PDBIO
from Bio.PDB.PDBParser import PDBParser

proteins = ['ALA','ARG','ASH','ASP','ASN','CYS','CYX','GLH','GLN','GLU','GLY','HID','HIE','HIP','HIS','ILE','LEU','LYN','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL']


# Argparse this...
try:
   pdbname = sys.argv[1]
except IndexError:
   print "PDB file required as first argument. "
   sys.exit()

try:
   sys.argv[2]
   print "Renaming termini (providing any character as a second argument will turn this on)"
   termini = True
except IndexError:
   termini = False
   print "Not renaming termini (add a second argument if you want to turn this on)"


parser = PDBParser(QUIET=True) # Lack of atom symbols = warnings if verbose
structure = parser.get_structure('ambpdb',pdbname)

# Convert first and last residues to termini
def rename_termini(struct):
   """Rename termini of chains to be NXXX and CXXX. TAKE CARE - RESULTS IN
      NON-STANDARD PDB FORMAT (residue names should be 3 characters).
      Input = Biopython structure object, returns structure object with
      modified first/last resnames in each chain."""

   for chain in struct.get_chains():
      first = chain.child_list[0]
      last = chain.child_list[-1]
      if first.resname in proteins:
         first.resname = 'N'+first.resname
      else:
         print """The first residue is either already an N-terminus, or doesn't
                  appear to be a natural amino acid. Continuing anyway."""
      if last.resname in proteins:
         last.resname = 'C'+last.resname
      else:   
         for idx,res in enumerate(chain.child_list,-1):
            if res.resname in proteins:
               continue
            else:
               last = chain.child_list[idx] # idx starts from -1
               print "You seem to have non-protein residues in a chain. Renaming residue %s %s as C-terminus" % (last.resname,last.id[1])
               last.resname = 'C'+last.resname
   return struct

def amb2gmx(struct):
   """Rename Amber atoms to those expected in the Gromacs amber parameter files.
      Mainly e.g. HB2/HB3 -> HB1/HB2. Cycles over all residues. Input =
      Biopython structure object, returns structure object with modified atom 
      names"""

   for residue in struct.get_residues():
      if residue.resname in ['GLY','NGLY','CGLY']:
         residue['HA2'].fullname = ' HA1'
         residue['HA3'].fullname = ' HA2'
      elif residue.resname in ['SER','NSER','CSER']:
         residue['HB2'].fullname = ' HB1'
         residue['HB3'].fullname = ' HB2'
      elif residue.resname in ['LEU','NLEU','CLEU']:
         residue['HB2'].fullname = ' HB1'
         residue['HB3'].fullname = ' HB2'
      elif residue.resname in ['ILE','NILE','CILE']:
         residue['HG12'].fullname = 'HG11'
         residue['HG13'].fullname = 'HG12'
      elif residue.resname in ['ASN','NASN','CASN']:
         residue['HB2'].fullname = ' HB1'
         residue['HB3'].fullname = ' HB2'
      elif residue.resname in ['GLN','NGLN','CGLN']:
         residue['HB2'].fullname = ' HB1'
         residue['HB3'].fullname = ' HB2'
         residue['HG2'].fullname = ' HG1'
         residue['HG3'].fullname = ' HG2'
      elif residue.resname in ['ARG','NARG','CARG']:
         residue['HB2'].fullname = ' HB1'
         residue['HB3'].fullname = ' HB2'
         residue['HG2'].fullname = ' HG1'
         residue['HG3'].fullname = ' HG2'
         residue['HD2'].fullname = ' HD1'
         residue['HD3'].fullname = ' HD2'
      elif residue.resname in ['HID','NHID','CHID','HIE','NHIE','CHIE','HIP','NHIP','CHIP','HIS','NHIS','CHIS']:
         residue['HB2'].fullname = ' HB1'
         residue['HB3'].fullname = ' HB2'
      elif residue.resname in ['TRP','NTRP','CTRP']:
         residue['HB2'].fullname = ' HB1'
         residue['HB3'].fullname = ' HB2'
      elif residue.resname in ['PHE','NPHE','CPHE']:
         residue['HB2'].fullname = ' HB1'
         residue['HB3'].fullname = ' HB2'
      elif residue.resname in ['TYR','NTYR','CTYR']:
         residue['HB2'].fullname = ' HB1'
         residue['HB3'].fullname = ' HB2'
      elif residue.resname in ['GLU','NGLU','CGLU','GLH','NGLH','CGLH']:
         residue['HB2'].fullname = ' HB1'
         residue['HB3'].fullname = ' HB2'
         residue['HG2'].fullname = ' HG1'
         residue['HG3'].fullname = ' HG2'
      elif residue.resname in ['ASP','NASP','CASP','ASH','NASH','CASH']:
         residue['HB2'].fullname = ' HB1'
         residue['HB3'].fullname = ' HB2'
      elif residue.resname in ['LYS','NLYS','CLYS','LYN','NLYN','CLYN']:
         residue['HB2'].fullname = ' HB1'
         residue['HB3'].fullname = ' HB2'
         residue['HG2'].fullname = ' HG1'
         residue['HG3'].fullname = ' HG2'
         residue['HD2'].fullname = ' HD1'
         residue['HD3'].fullname = ' HD2'
         residue['HE2'].fullname = ' HE1'
         residue['HE3'].fullname = ' HE2'
      elif residue.resname in ['PRO','NPRO','CPRO']:
         residue['HB2'].fullname = ' HB1'
         residue['HB3'].fullname = ' HB2'
         residue['HG2'].fullname = ' HG1'
         residue['HG3'].fullname = ' HG2'
         residue['HD2'].fullname = ' HD1'
         residue['HD3'].fullname = ' HD2'
      elif residue.resname in ['CYS','NCYS','CCYS','CYX','NCYX','CCYX']:
         residue['HB2'].fullname = ' HB1'
         residue['HB3'].fullname = ' HB2'
      elif residue.resname in ['MET','NMET','CMET']:
         residue['HB2'].fullname = ' HB1'
         residue['HB3'].fullname = ' HB2'
         residue['HG2'].fullname = ' HG1'
         residue['HG3'].fullname = ' HG2'
   return struct
  
def write_struct(struct,name):
   """Writes a corrected PDB file to disk with 'GMX_' prepended"""

   w = PDBIO()
   w.set_structure(struct)
   w.save('GMX_'+name)

#### Main below here ####

if termini is True:
   structure = rename_termini(structure) 
print """If you have separate chains in your PDB (separated by a TER card),
         make sure they have a chain identifier in column 22, in accordance
         with the PDB version 3 format. Biopython is expecting this!"""
structure = amb2gmx(structure)
write_struct(structure,pdbname)
print "Done."
