#!/usr/bin/env python

# Script to convert & translate atm positions from Jacek's xyz to Amber rst
# Can then be converted further to netcdf eg. cpptraj

import parmed as pmd
import sys


# Read in. Arg 1 = Amber inpcrd, arg 2 = Tinker xyz
ambfile = pmd.amber.Rst7(sys.argv[1])
tinkfile = pmd.tinker.tinkerfiles.XyzFile(sys.argv[2])

# Copy coords across for all atoms
ambfile.coordinates = tinkfile.coordinates

# Dict for translation (manually curated)
# Key = amber atom, Val = tinker atom (only valid for Pro_dipep)
amb_2_tink = {1:4,2:1,3:5,4:6,5:2,6:3,7:7,8:14,9:19,10:20,11:13,12:17,13:18,14:12,15:15,16:16,17:8,18:11,19:9,20:10,21:21,22:23,23:22,24:24,25:25,26:26}

# Translate & write
for key,val in amb_2_tink.items():
   ambfile.coordinates[key-1] = tinkfile.coordinates[val-1]
ambfile.write("Parmedout.rst")


