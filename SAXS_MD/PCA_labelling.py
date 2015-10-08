#!/usr/bin/env python

# Script to create entries to label PCA plots in R

import sys,csv


# This reads in the file PDB_list.txt, which has a separate PDB ID on each line
# Change the file name as desired

pdbfile = open("PDB_list_MD_T.txt",'r')
pdblist = pdbfile.read().splitlines()
pdbfile.close()

# Reads in a csv file as a csv object

csvfile = open("MD_CRYSOL_VALUES_ordered.csv", 'rb')
csvread = csv.reader(csvfile)

# This defines the principal components you want on the x and y axes
# Change as desired
x = 1
y = 3

# Defines endpoints of ranges for colour definitions
r1 = 3.0 
r2 = 6.0
r3 = 9.0

# Below here shouldn't need changing, except if you want to define new colours

outfile = open("R_PCA_entries.txt",'w')
i = 1

for pdb in pdblist:
   row = csvread.next()
   crysol = float(row[7])
   if int(row[1]) != i:
      print "The MD frame at position %s is actually labelled %s ! Check the inputs!" % (i,row[1])
      outfile.close()
      csvfile.close()
      sys.exit(1)
   outfile.write("x6G08_MD_%s_raw<-scan(\"/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_%s.txt\")\n" % (pdb,pdb) )
   outfile.write("x6G08_MD_%s_proj<-project.pca(x6G08_MD_%s_raw,pcadata)\n" % (pdb,pdb) )
   outfile.write("points(x6G08_MD_%s_proj[%d],x6G08_MD_%s_proj[%d],pch=20)\n" % (pdb,x,pdb,y) )
   if crysol < r1:
      outfile.write("text(x6G08_MD_%s_proj[%d],x6G08_MD_%s_proj[%d],pos=4,label=\"%s\",col=\"black\")\n" % (pdb,x,pdb,y,pdb) )
   elif r1 <= crysol < r2:
      outfile.write("text(x6G08_MD_%s_proj[%d],x6G08_MD_%s_proj[%d],pos=4,label=\"%s\",col=\"blue\")\n" % (pdb,x,pdb,y,pdb) )
   elif r2 <= crysol < r3:
      outfile.write("text(x6G08_MD_%s_proj[%d],x6G08_MD_%s_proj[%d],pos=4,label=\"%s\",col=\"green\")\n" % (pdb,x,pdb,y,pdb) )
   elif crysol >= r3:
      outfile.write("text(x6G08_MD_%s_proj[%d],x6G08_MD_%s_proj[%d],pos=4,label=\"%s\",col=\"red\")\n" % (pdb,x,pdb,y,pdb) )
   else:
      print "Something wrong in the csv file - column 8 doesn't see to correspond to crysol values"
      outfile.close()
      csvfile.close()
      sys.exit(1)
   outfile.write("\n")
   i += 1

outfile.close()
csvfile.close()
