#!/usr/bin/env python

# Reading in and converting TukeyHSDs to heatmaps
# Requirements: numpy,matplotlib

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.ticker import AutoMinorLocator


paramdict = {
"New_params":0,
"New_polgps":1,
"Multifit":2,
"Aug_cc":3,
"Wat14":4,
"OH_scale":5,
"Gaff":6
}

def mapsort(col):
      """Splits column into list of parameter sets, maps
      to a set of indices in an array and sorts (x index is largest)"""
      setlist = col.split("-")
      mappedlist = [ paramdict[i] for i in setlist ]
      mappedlist.sort()
      return np.array(mappedlist)

def map_ci(val,conf):
      """Returns 1 if p value < confidence interval, 0 otherwise"""
      if val < conf:
         return 1
      else:
         return 0      

def tukey_array(idx,fn,p=0.05):
      """Generates and returns a half 7x7 array of significant pairwise
      p-values as boolean true/false. Requires file index name, function to 
      convert into booleans, and desired p-value for significance (default
      0.05)"""
      # Read in columns 1 and 5, convert 1 to a coordinate array
      tukeyres = np.genfromtxt("TukeyHSD_%d.txt" % idx,skip_header=7,skip_footer=1,usecols=(0,4),converters={0: lambda col: mapsort(col)},dtype=[object,float])
      # Initialise array of non-significant p-values
      pvals = np.ones((7,7))
      # Populate array
      for coords,currval in tukeyres:
         pvals[coords[0],coords[1]] = currval
      # Boolify with function
      return fn(pvals,p)

def saveplot(arr):
      """Saves a heatmap of the desired array"""
      plt.figure()
      plt.imshow(arr,interpolation='None',cmap=cm.hot_r,extent=[1,8,8,1],vmax=47)
      plt.gca().xaxis.set_minor_locator(AutoMinorLocator(2)) # Minor ticks centred between majors
      plt.gca().yaxis.set_minor_locator(AutoMinorLocator(2)) # Minor ticks centred between majors
      for tick in plt.gca().xaxis.get_minor_ticks():
         tick.tick1line.set_markersize(0)
         tick.tick2line.set_markersize(0)
         tick.label1.set_horizontalalignment('center')
      for tick in plt.gca().yaxis.get_minor_ticks():
         tick.tick1line.set_markersize(0)
         tick.tick2line.set_markersize(0)
         tick.label1.set_horizontalalignment('right')
      plt.gca().set_xticklabels(['New params','New polgps','Multifit','Aug-cc','Wat14','OH-scale','GAFF'],rotation=90,minor=True)
      plt.gca().set_yticklabels(['New params','New polgps','Multifit','Aug-cc','Wat14','OH-scale','GAFF'],minor=True)
      plt.gca().set_xticklabels([]) # Switch off majors
      plt.gca().set_yticklabels([]) # Switch off majors
      # Overlay text of cell values
      min_ax, max_ax = 1,8
      idx_array = np.arange(min_ax + 0.5, max_ax + 0.5, 1.0)
      x,y = np.meshgrid(idx_array,idx_array) # Two grids of x & y coords
      for i, (x_val, y_val) in enumerate(zip(x.flatten(), y.flatten())):
         if x_val < y_val:
            continue
         else:
            txt = int(heatmap.flatten()[i]) # Only for top half of fig
            plt.gca().text(x_val, y_val, txt, va='center', ha='center',color='white') # Plot text labels
      plt.suptitle('Tukey HSD significance map, parameter sets 1-7',y=0.98,fontsize=14)
      plt.colorbar()
      plt.savefig('Tukey_HSD_heatmap.eps',bbox_inches='tight',dpi=300)


### main below here ###
if __name__ == "__main__":
   pval_all = np.zeros((7,7))
   boolify = np.vectorize(map_ci)
   for lig in range(1,48):
      pval_curr = tukey_array(lig,boolify,p=0.05)
      pval_all = np.dstack((pval_all,pval_curr))
   heatmap = np.sum(pval_all,axis=2)
   saveplot(heatmap)
   np.savetxt("Data_for_heatmap.txt",heatmap,fmt='%d')

