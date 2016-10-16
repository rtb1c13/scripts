#!/usr/bin/env python
# Author: Richard Bradshaw, R.T.Bradshaw@soton.ac.uk

# Script to create 3D matplotlib plots of e.g. field vectors.
# Vectors may need rotation into reference frame, i.e multiplication
# by a rotation matrix

# Requirements: numpy, matplotlib v2.0.0b1 or later

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from mpl_toolkits.mplot3d import Axes3D
import argparse, itertools

parser = argparse.ArgumentParser()
parser.add_argument("-f","--files",help="Files with fields",nargs='+',type=str,required=True)
parser.add_argument("-c","--cols",help="Column indices to use from fields files. Defaults to 0,1,2",nargs='+',type=int,default=(0,1,2))
parser.add_argument("-r","--rotmats",help="Files with rotation matrices",nargs='+',type=str)



def combine_fields(files,cols):
   """Reads in a list of files and returns concatenated
      Numpy array."""
   cat = np.genfromtxt(files[0],usecols=cols)
   for fn in files[1:]:
      x = np.genfromtxt(fn,usecols=cols)
      cat = np.concatenate((cat,x))
   return cat

def combine_rotmats(files):
   """Requires any number of rotation matrices (as numpy arrays
      and returns concatenated Numpy array of all matrices."""
   cat = np.genfromtxt(files[0],comments="Rotation")
   for fn in files[1:]:
      x = np.genfromtxt(fn,comments="Rotation")
      cat = np.concatenate((cat,x))
   rotmats = np.reshape(cat,(len(cat)/3,3,3))
   return rotmats
   
def rotate_vectors(vecs,mats):
   """Pairwise rotates an array of vectors using an array of
      rotation matrices (performs pairwise dot product of
      vecs . mats). Usage: rotate_vectors(vecs,mats)"""
   newvecs = np.zeros(np.shape(vecs))
   for i,vec in enumerate(vecs):
      newvecs[i] = np.dot(vec,mats[i])
   return newvecs

def plot_quiver(fields):
   """TBC"""
   size = np.shape(fields)
   origins = np.zeros(size)
   X,Y,Z = zip(*origins)
   U,V,W = zip(*fields)
   fig = plt.figure()
   ax = fig.add_subplot(111,projection='3d')
   ax.quiver(X,Y,Z,U,V,W,arrow_length_ratio=0.1,pivot='tail',linewidth=0.2,color='blue')
   ax.quiver(0,0,0,50,0,0,arrow_length_ratio=0.1,pivot='tail',linewidth=3,color='yellow')
   ax.quiver(0,0,0,0,50,0,arrow_length_ratio=0.1,pivot='tail',linewidth=3,color='yellow')
   ax.quiver(0,0,0,0,0,50,arrow_length_ratio=0.1,pivot='tail',linewidth=3,color='yellow')
   ax.set_xlabel("X component / MV cm$^{-1}$")
   ax.set_ylabel("Y component / MV cm$^{-1}$")
   ax.set_zlabel("Z component / MV cm$^{-1}$")
   ax.set_xlim((-50,50))
   ax.set_ylim((-50,50))
   ax.set_zlim((-50,50))
   plt.show()

def plot_scatter(fields,crosses):
   """TBC"""
   fig = plt.figure()
   ax = fig.add_subplot(111,projection='3d')
   normfields = fields/np.reshape(np.linalg.norm(fields,axis=1),(len(fields),1))
   X,Y,Z = zip(*normfields)
   U,V,W = zip(*crosses)
   fig = plt.figure(figsize=(16,12))
   ax = fig.add_subplot(111,projection='3d')
   x_ax = ax.quiver(-1,0,0,2,0,0,arrow_length_ratio=0.1,pivot='tail',linewidth=3,color='yellow')
   y_ax = ax.quiver(0,-1,0,0,2,0,arrow_length_ratio=0.1,pivot='tail',linewidth=3,color='yellow')
   z_ax = ax.quiver(0,0,-1,0,0,2,arrow_length_ratio=0.1,pivot='tail',linewidth=3,color='yellow')
   dots = ax.scatter(U,V,W,color='blue')
   tri = ax.scatter(X,Y,Z,c=np.linalg.norm(fields,axis=1),cmap='inferno',marker='^')
   ax.set_xlabel("Normalized X component")
   ax.set_ylabel("Normalized Y component")
   ax.set_zlabel("Normalized Z component")
   ax.set_xlim((-1,1))
   ax.set_ylim((-1,1))
   ax.set_zlim((-1,1))
   cbar = fig.colorbar(tri)
   cbar.set_label("Total field strength / MV cm$^{-1}$",fontsize=16)
   plt.legend((x_ax,dots,tri),['xyz axes','Pro ring normal','Field'])
   plt.title('Field distribution about proline N',fontsize=20)
   for ang in xrange(0,360,10):
      ax.view_init(azim=ang)
      plt.savefig("Pro_N_"+str(ang)+".png")

def write_avgs(farray,fn="mean_field.txt",ax=0):
   """Writes out mean of array (default across axis 0) to text file"""
#  Reshaped to 2D to print out columns delimited by ' ' properly
   np.savetxt(fn,np.reshape(np.mean(farray,axis=ax),(1,3)),fmt='%10.5f')

def rot_and_plot():
   """Main function - rotates and plots quiver graphs of vectors"""
   fields = combine_fields(args.files,args.cols)
   crosses = combine_fields(args.files,(3,4,5))
   rotmats = combine_rotmats(args.rotmats)
   rotatedfields = rotate_vectors(fields,rotmats)
   rotatedcrosses = rotate_vectors(crosses,rotmats)
   plot_quiver(rotatedfields)
   plot_scatter(rotatedfields,rotatedcrosses)
   write_avgs(rotatedfields,fn="Average_field.txt")
   write_avgs(rotatedcrosses,fn="Average_crossprod.txt")


# Main below here
if __name__ == '__main__':
   global args
   args = parser.parse_args()
   rot_and_plot()


