#!/usr/bin/env python

# Script to sum segment averages and take diffs Out - In 

import recipy
import numpy as np
### For headless nodes, recipy has already imported matplotlib, so backend is already set
import matplotlib 
matplotlib.use('agg')
###
import matplotlib.pyplot as plt
### More for headless nodes with recipy
plt.switch_backend('agg')
###
from scipy.stats import sem
from scipy.stats import pearsonr as correl
from cycler import cycler

### Define defaults for matplotlib plots
plt.rc('lines', linewidth=1.5, markersize=4)
plt.rc('axes', prop_cycle=(cycler('color', ['k','b','r','orange','c','m','y','g'])), # Color cycle defaults to black
       labelweight='heavy', labelsize=22, titlesize=22) # Default fontsizes for printing
plt.rc('axes.spines', top=False, right=False) # Switch off top/right axes
plt.rc('legend', fontsize=16) # Default fontsizes for printing
plt.rc('xtick', labelsize=16) # Default fontsizes for printing
plt.rc('ytick', labelsize=16) # Default fontsizes for printing
plt.rc('figure', figsize=(11,8.5), titlesize=22, titleweight='heavy') # Default fontsizes for printing
#plt.rc('text', usetex=True)

times = [0.167, 1., 10., 120.]
#
out_ca = np.loadtxt("segfracs_1.dat", dtype=[('segs', np.int32, (2,)), ('fracs', np.float64, (4,))])
out_cb = np.loadtxt("segfracs_2.dat", dtype=[('segs', np.int32, (2,)), ('fracs', np.float64, (4,))])
out_cc = np.loadtxt("segfracs_3.dat", dtype=[('segs', np.int32, (2,)), ('fracs', np.float64, (4,))])
out_expt = np.loadtxt("segfracs_expt.dat", dtype=[('segs', np.int32, (2,)), ('fracs', np.float64, (4,))])

print "Read-in segs equal?", np.array_equal(out_ca['segs'][:,0:2], out_cb['segs'][:,0:2])
out_mean = np.mean((out_ca['fracs'], out_cb['fracs'], out_cc['fracs']), axis=0)
out_err = sem((out_ca['fracs'], out_cb['fracs'], out_cc['fracs']), axis=0)

# Plot
fig, axs = plt.subplots(2,2, figsize=(22,17))
xs = range(1,len(out_ca['segs'])+1)
for time, ax, j in zip(times, axs.flatten(), range(len(times))):
    xs = range(1,len(out_ca['segs'])+1)
    labels = [ str(i[0])+"-"+str(i[1]) for i in out_ca['segs'][:,0:2] ]
    ax.set_title("Time = %s min" % time)
    ax.set_xticks(xs)
    ax.set_xticklabels(labels, rotation='vertical')
    ax.set_ylabel("Deuterated fraction")
    ax.set_xlabel("Peptide segment")
    ax.set_ylim(0.0, 1.0)
    ax.plot(xs, out_mean[:,j], label="Predicted fraction, R = %3.2f" % correl(out_mean[:,j], out_expt['fracs'][:,j])[0])
    ax.plot(xs, out_expt['fracs'][:,j], label="Experimental fraction", linestyle=':')
    ax.fill_between(xs, out_mean[:,j]+out_err[:,j], out_mean[:,j]-out_err[:,j], color='gray', alpha=0.2, label="Std. err. across repeats")
    ax.legend()

plt.tight_layout()
np.savetxt("mean_segfracs.dat", out_mean)
plt.savefig("seg_curves_combined.png", dpi=120)

