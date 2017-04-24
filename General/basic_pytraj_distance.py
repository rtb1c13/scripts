#!/usr/bin/env

# Basic script using pytraj to read & plot distance between COM of two residue masks
# Uses matplotlib, pytraj (Amber16 version)
# Usage: basic_pytraj_distance.py trajectory prmtop res1 res2

import pytraj as pt
import matplotlib.pyplot as plt
import sys

trajname, topname = (sys.argv[1], sys.argv[2])
mask1, mask2 = (sys.argv[3], sys.argv[4])

traj = pt.load(trajname, top=topname)
dist = pt.calc_distance(traj, mask=mask1+' '+mask2)


# Show plot
plt.plot(dist)
plt.show()

# Below here is manually inserted



