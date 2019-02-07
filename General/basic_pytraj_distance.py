#!/usr/bin/env python

# Basic script using pytraj to read & plot distance between COM of two residue masks
# Uses matplotlib, pytraj (Amber16 version)
# Usage: basic_pytraj_distance.py trajectory prmtop res1 res2 [res3 res4]

import pytraj as pt
import matplotlib.pyplot as plt
import numpy as np
import sys

trajname = sys.argv[1]
topname = sys.argv[2]
mask1 = sys.argv[3]
mask2 = sys.argv[4]

traj = pt.load(trajname, top=topname)
dist = pt.calc_distance(traj, mask=mask1+' '+mask2)


# Show plot
#plt.plot(dist)
#plt.show()

# Below here is manually inserted
cryst = pt.load('../../3OV6_DOPC_solv.inpcrd',top=topname)
cryst_dist = pt.calc_distance(cryst, mask=mask1+' '+mask2)

try:
 mask3, mask4 = (sys.argv[5], sys.argv[6])
 dist2 = pt.calc_distance(traj, mask=mask3+' '+mask4)
 cryst_dist2 = pt.calc_distance(cryst, mask=mask3+' '+mask4)
except IndexError:
 print 'You forgot the second pair of residues...'
x = np.arange(0, len(dist)/2, 0.5) # For snapshots every 0.5 ns
plt.plot(x, dist, color='black')
plt.plot(x, dist2, color='blue')
line1 = [ cryst_dist[0] ] * len(dist)
line2 = [ cryst_dist2[0] ] * len(dist)
plt.plot(x, line1, ls='--', color='black')
plt.plot(x, line2, ls='--', color='blue')

plt.gca().set_xlabel('Time / ns')
plt.gca().set_xlim((0,len(dist)/2)) # For snapshots every 0.5 ns
plt.gca().set_ylabel(r'Distance / $\AA$')
plt.gca().legend([mask1+'-'+mask2, mask3+'-'+mask4])
plt.gca().set_title('Distance between roof residues (initial values = %5.2f, %5.2f)' % (cryst_dist[0], cryst_dist2[0]))
plt.savefig('roof_res_distances.png',dpi=300)

# Save text files
np.savetxt(mask1+'-'+mask2+'_dist.txt',dist)
np.savetxt(mask3+'-'+mask4+'_dist.txt',dist2)

