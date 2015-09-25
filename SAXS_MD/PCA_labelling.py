#!/usr/bin/env python

# Script to create entries to label PCA plots in R


import sys

# This reads in the file PDB_list.txt, which has a separate PDB ID on each line
# Change the file name as desired

pdbfile = open("PDB_list.txt",'r')
pdblist = pdbfile.read().splitlines()

# This defines the principal components you want on the x and y axes
# Change as desired
x = 1
y = 2

# Below here shouldn't need changing!

outfile = open("R_PCA_entries.txt",'w')

for pdb in pdblist:
   outfile.write("x%s_raw<-scan(\"../Structures/C-alpha_only/PCA_%s_1.txt\")\n" % (pdb,pdb) )
   outfile.write("x%s_proj<-project.pca(x%s_raw,pcadata)\n" % (pdb,pdb) )
   outfile.write("points(x%s_proj[%d],x%s_proj[%d],pch=20)\n" % (pdb,x,pdb,y) )
   outfile.write("text(x%s_proj[%d],x%s_proj[%d],pos=4,label=\"%s\",col=\"purple\")\n" % (pdb,x,pdb,y,pdb) )
   outfile.write("\n")

outfile.close()

