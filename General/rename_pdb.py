#!/usr/bin/env python

# Script to rename atoms in PDB based on 'original' MOE fit conformation,
# with coordinates substituted for 'new' conformation
# Uses Parmed to print PDB in Amber compatible format
# Usage: rename_pdb.py old_file new_file 

import parmed as pmd
import sys

old = pmd.load_file(sys.argv[1])
new = pmd.load_file(sys.argv[2])
old.write_pdb("Renamed_"+sys.argv[2],coordinates=new.coordinates)
