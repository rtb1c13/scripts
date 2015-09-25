#!/usr/bin/env python

# Script to print out scaled quadrupole entries for a list of atom types (user-dfeined)

import sys

if len(sys.argv) < 3:
  print "Usage = scale_OH.py [FILENAME] [ATOM TYPES]"

# Read in input file and atom types
infile = sys.argv[1]
attypes = []
for i in sys.argv[2:]:
  attypes.append(int(i))

# I/O files
fi = open(infile,'r')
fo = open("scaled_OH.txt",'w')

# Read in inpult file sequentially
while True:
  line = fi.readline()
  if not line: break

  # Search for multipole lines and append to list
  if line.find("multipole") > -1:
    multipole_lines = [line]
    i=1
    while i < 5:
      line = fi.readline()
      multipole_lines.append(line)
      i = i+1
#   Match up multipoles to scale
    if int(multipole_lines[0].split()[1]) in attypes:
      fo.write(multipole_lines[0])
      fo.write(multipole_lines[1])
#     Scale OH quadrupoles
      xx = float(multipole_lines[2].split()[0]) * 0.6
      xy = float(multipole_lines[3].split()[0]) * 0.6
      yy = float(multipole_lines[3].split()[1]) * 0.6
      xz = float(multipole_lines[4].split()[0]) * 0.6
      yz = float(multipole_lines[4].split()[1]) * 0.6
      zz = float(multipole_lines[4].split()[2]) * 0.6
#     Print OH quadrupoles
      fo.write("%47.5f\n" % xx)
      fo.write("%47.5f %10.5f\n" % (xy,yy))
      fo.write("%47.5f %10.5f %10.5f\n" % (xz,yz,zz))


fi.close()
fo.close()
print "Scaling complete, printed out in scaled_OH.txt"
print "Substitute this for your existing multipoles"
