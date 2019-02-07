#!/usr/bin/env python

# Script to split up fields byframe into 4 subsets
# 1) Doubly H-bonded (<3.0A)
# 2) Tyr16 H-bonded
# 3) Ash103 H-bonded
# 4) No H-bonds

import numpy as np

times = range(0,30,5)

for tp in times:
    fields = np.genfromtxt("%sns_fields.txt" % tp)
    tyrd = np.genfromtxt("TyrHO-NthO_%sns.out" % tp, usecols=1)
    ashd = np.genfromtxt("AshHO-NthO_%sns.out" % tp, usecols=1)
    tyrd = np.repeat(tyrd, 5000, axis=0) # Turn struct snapshots every 10ps into every 2fs
    ashd = np.repeat(ashd, 5000, axis=0) # Turn struct snapshots every 10ps into every 2fs
    fields,tyrd,ashd = [ np.reshape(i,(500000,1)) for i in [fields,tyrd,ashd] ] # Reshape
    raw_data = np.concatenate((fields,tyrd,ashd),axis=1)
    # Doubly H-bonded
    tmp = raw_data[raw_data[:,1]<3.0]
    try:
        doubly = tmp[tmp[:,2]<3.0]
    except IndexError:
        doubly = np.array([0,0,0],ndmin=2)
        pass
    np.savetxt("Doubly_H-bonded_%sns.txt" % tp,doubly,fmt='%8.5f')
    # Tyr16 only
    try:
        tyr16 = tmp[tmp[:,2]>=3.0]
    except IndexError:
        tyr16 = np.array([0,0,0],ndmin=2)
        pass
    np.savetxt("Tyr16_H-bonded_%sns.txt" % tp,tyr16,fmt='%8.5f')
    # Ash103 only
    tmp = raw_data[raw_data[:,2]<3.0]
    try:
        ash103 = tmp[tmp[:,1]>=3.0]
    except IndexError:
        ash103 = np.array([0,0,0],ndmin=2)
        pass
    np.savetxt("Ash103_H-bonded_%sns.txt" % tp,ash103,fmt='%8.5f')
    # Neither
    tmp = raw_data[raw_data[:,2]>=3.0]
    try:                                 # In case of no entries
        neither = tmp[tmp[:,1]>=3.0]
    except IndexError:
        neither = np.array([0,0,0],ndmin=2)
        pass
    np.savetxt("No_H-bonded_%sns.txt" % tp,neither,fmt='%8.5f')

    arrs = [doubly,tyr16,ash103,neither]
    with open("Hbond_stats_%sns.txt" % tp,'a') as f:
        f.write("Double, Tyr16, Ash103, Neither\n")
        f.write("No. entries, Mean field, Std Dev, Std Err\n")
        for i in arrs:
            f.write("%5d %8.3f %5.3f %5.3f\n" % (len(i[:,0]),np.mean(i[:,0]),np.std(i[:,0]),np.std(i[:,0])/np.sqrt(len(i[:,0]))))
