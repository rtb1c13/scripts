#!/usr/bin/env python

# Script to split up fields byframe into 4 subsets
# 1) Doubly H-bonded (<3.0A)
# 2) Tyr16 H-bonded
# 3) Ash103 H-bonded
# 4) No H-bonds

import numpy as np

# Input data should be ordered: Fields (1 col), Tyr16 distances (2 cols), Ash103 distances (2 cols)
raw_data = np.genfromtxt("fields_Tyr_Ash.out",usecols=(0,2,4))
# Doubly H-bonded
tmp = raw_data[raw_data[:,1]<3.0]
doubly = tmp[tmp[:,2]<3.0]
np.savetxt("Doubly_H-bonded_amoeba.txt",doubly,fmt='%8.5f')
# Tyr16 only
tyr16 = tmp[tmp[:,2]>=3.0]
np.savetxt("Tyr16_H-bonded_amoeba.txt",tyr16,fmt='%8.5f')
# Ash103 only
tmp = raw_data[raw_data[:,2]<3.0]
ash103 = tmp[tmp[:,1]>=3.0]
np.savetxt("Ash103_H-bonded_amoeba.txt",ash103,fmt='%8.5f')
# Neither
tmp = raw_data[raw_data[:,2]>=3.0]
neither = tmp[tmp[:,1]>=3.0]
np.savetxt("No_H-bonded_amoeba.txt",neither,fmt='%8.5f')

arrs = [doubly,tyr16,ash103,neither]
with open("Hbond_stats_amoeba.txt",'a') as f:
   f.write("Double, Tyr16, Ash103, Neither\n")
   f.write("No. entries, Mean field, Std Dev, Std Err\n")
   for i in arrs:
      f.write("%5d %8.3f %5.3f %5.3f\n" % (len(i[:,0]),np.mean(i[:,0]),np.std(i[:,0]),np.std(i[:,0])/np.sqrt(len(i[:,0]))))
