#/usr/bin/python2.6


# Short script to print out SAXS profile calculated using the saxs module in IMP
# Note this should be the same as the intensities calculated with a Debye double sum
# script, give or take some error due to the use of distance distribution function
# in IMP instead of pairwise distances as in Debye formula.

# Output should be the same as running foxs with -e 1.0 -w 0.0

# NOTES: Only runs with python 2.6 (IMP compiled this way locally), make sure that the IMP
# site-packages are in your PYTHONPATH. First command line argument should be PDB input


# mostly taken from profile.py in imp/modules/saxs/examples
import IMP
import IMP.atom
import IMP.core
import IMP.saxs
import sys,os

m = IMP.kernel.Model()

# read PDB (can adapt PDB selector to include H specifically)
mp = IMP.atom.read_pdb(sys.argv[1], m, IMP.atom.NonWaterNonHydrogenPDBSelector())

# select particles from the model
particles = IMP.atom.get_by_type(mp, IMP.atom.ATOM_TYPE)

# calculate SAXS profile
model_profile = IMP.saxs.Profile()
model_profile.calculate_profile(particles)
model_profile.write_SAXS_file(sys.argv[1]+'.dat')


