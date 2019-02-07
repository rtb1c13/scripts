#!/usr/bin/env python

# Script to convert Tinker PDB atom namings to Amber
# Note this won't change coordinates/atom ordering. To 
# get atom ordering consistent with an Amber prmtop, the 
# resulting PDB should be read into Leap & re-saved.
# Usage: tinker_pdb_to_amber.py input_pdb.pdb output_pdb.pdb 

import mdtraj as md
import sys

inpdb = sys.argv[1]
outpdb = sys.argv[2]

p = md.load(inpdb)
# Rename water H atoms
for r in p.top.residues:
    if r.name =='HOH':
        r.atom(1).name = 'H1'
        r.atom(2).name = 'H2'

# Rename N-term H atoms
for c in p.top.chains:
    for a in c.residue(0).atoms:
        if a.name == 'H':
            a.name = 'H1'

# Rename His to Hid/Hie/Hip
for r in p.top.residues:
    if r.name =='HIS':
        try:
            r.atom('HD1')
            r.atom('HE2')
            r.name = 'HIP'
            continue
        except KeyError:
            pass
        try:
            r.atom('HD1')
            r.name = 'HID'
            continue
        except KeyError:
            pass
        try:
            r.atom('HE2')
            r.name = 'HIE'
            continue
        except KeyError:
            print "Residue %s should be a histidine but doesn't have HD1 or HE2" % r

p.save(outpdb)
