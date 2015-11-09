#!/usr/bin/env python2.7
# Script to convert Tinker potential or Gaussian cubefile into .dx format
# for visualisation of ele potential around system.
# Prepared by S Genheden, adapted by RTB

import os

import numpy as np
from scipy.stats import norm

def _init_grid(xyz,spacing,padding) :
  """
  Initialize a grid based on a list of x,y,z coordinates

  Parameters
  ----------
  xyz  : Numpy array 
    Cartesian coordinates that should be covered by the grid
  spacing : float 
    the grid spacing
  padding : float 
    the space to add to minimum extent of the coordinates

  Returns
  -------
  Numpy array
    the grid 
  list of Numpy arrays
    the edges of the grid
  """
  
  origin  = np.floor(xyz.min(axis=0))-padding
  tr      = np.ceil(xyz.max(axis=0))+padding
  length  = tr-origin
  shape  =  np.array([int(l/spacing + 0.5) + 1 for l in length],dtype=int)
  grid    = np.zeros(shape)
  edges  = [np.linspace(origin[i],tr[i],shape[i]) for i in range(3)]
  return grid,edges

def _voxel(coord,edges) :
  """
  Wrapper for the numpy digitize function to return the grid coordinates
  """
  return np.array([np.digitize(coord,edges[i])[i] for i in range(3)],dtype=int)
 
def writeDX(grid,origin,spacing,filename) :
  """
  Write the grid to file in DX-format

  Parameters
  ----------
  grid : Numpy array 
    the 3D grid
  origin : NumpyArray 
    the bottom-left coordinate of the grid
  spacing  : float
    the grid spacing
  filename : string
     the name of the DX file
  """
  f = open(filename, 'w')
  f.write("object 1 class gridpositions counts %5d%5d%5d\n"%(grid.shape[0],grid.shape[1],grid.shape[2]))
  f.write("origin %9.4f%9.4f%9.4f\n"%(origin[0],origin[1],origin[2]))
  f.write("delta %10.7f 0.0 0.0\n"%spacing)
  f.write("delta 0.0 %10.7f 0.0\n"%spacing)
  f.write("delta 0.0 0.0 %10.7f\n"%spacing)
  f.write("object 2 class gridconnections counts %5d%5d%5d\n"%(grid.shape[0],grid.shape[1],grid.shape[2]))
  f.write("object 3 class array type double rank 0 items  %10d data follows\n"%(grid.shape[0]*grid.shape[1]*grid.shape[2]))
  cnt = 0
  for x in range(grid.shape[0]) :
    for y in range(grid.shape[1]) :
      for z in range(grid.shape[2]) :
        f.write("%19.10E"%grid[x,y,z])
        cnt = cnt + 1
        if cnt >= 3 :
          cnt = 0
          f.write("\n")
  if cnt > 0 : f.write("\n")
  f.write('attribute "dep" string "positions"\n')
  f.write('object "regular positions regular connections" class field\n')
  f.write('component "positions" value 1\n')
  f.write('component "connections" value 2\n')
  f.write('component "data" value 3\n')
  f.close()

def make_density(xyz,pot,padding=2.0,spacing=0.25) :
  """
  Make a rectangular density of potentials from a set of coordinates
  
  Parameters
  ----------
  xyz : numpy array
    the grid points
  pot : numpy array
    the potential in each grid point
  padding : float, optional
    extra space to add to the minimum extent 
  spacing : float, optional
    the spacing of the grid
 
  Returns
  -------
  NumpyArray 
    the 3D density
  dictionary
    grid properties, keys = spacing, min, and max
  """
  
  # Create a grid
  grid,edges = _init_grid(xyz,spacing,padding)
  
  # Put the potential on the grid.
  # Gaussian filtered using defaults in scipy function
  for point,potential in zip(xyz,pot) :
    v = _voxel(point,edges)
    grid[v[0],v[1],v[2]] = potential
  from scipy.ndimage.filters import gaussian_filter
  grid = gaussian_filter(grid,1.0)

  prop = {}
  prop["spacing"] = spacing
  prop["min"] = [e[0] for e in edges]
  prop["max"] = [e[-1] for e in edges]
  return grid,prop

def read_cubefile(filename) :
  """
  Read a Cube potential file
  
  Parameters
  ----------
  filename : string
    the name of the file to read

  Returns
  -------
  Nx3 NumpyArray
    the grid coordinates
  N NumpyArray
    the potential at each grid coordinate
  """

  xyz = []
  pot = []
  atomxyz = []
  with open(filename,'r') as f :
    line = f.readline() # Skip header #1
    line = f.readline() # Skip header #2
    natom = int(f.readline().strip().split()[0]) # Read number of atoms
    # Skip 3 lines of number of # voxels and axis vector
    line = f.readline()
    line = f.readline()
    line = f.readline()
    # Read atom lines, but ignore the data
    for i in range(natom) : 
      atomxyz.append(f.readline().strip().split()[1:])
    # Now read the interesting data
    line = f.readline()
    while line :
      cols = line.strip().split()
      xyz.append(cols[0:3])
      pot.append(cols[3])
      line = f.readline()
  
  #h,t = os.path.splitext(filename)
  #with open(h+".xyz",'w') as f :
  #  f.write("%d\n1\n"%natom)
  #  for a in atomxyz : f.write("%d %8.3f %8.3f %8.3f\n"%(int(float(a[0])),float(a[1])*0.529,float(a[2])*0.529,float(a[3])*0.529))

  return np.array(xyz,float),np.array(pot,float)

def read_potfile(filename) :
  """
  Read a Tinker potential file
  
  Parameters
  ----------
  filename : string
    the name of the file to read

  Returns
  -------
  Nx3 NumpyArray
    the grid coordinates
  N NumpyArray
    the potential at each grid coordinate
  """
  xyz = []
  pot = []
  with open(filename,'r') as f :
    line = f.readline() # Skip header
    line = f.readline()
    while line :
      cols = line.strip().split()
#     +/-10A box dimensions about origin
      if abs(float(cols[1])) > 10 or abs(float(cols[2])) > 10 or abs(float(cols[3])) > 10 :
        line = f.readline()
      else :
         xyz.append(cols[1:4])
         pot.append(cols[4])
         line = f.readline()
  
  return np.array(xyz,float),np.array(pot,float)
    

#
# If this is run from the command-line
#      
if __name__ == "__main__":

  import argparse

  # Setup a parser of the command-line arguments
  parser = argparse.ArgumentParser(description="Program put cube/pot file on a rectangular grid")
  parser.add_argument('-f','--files',nargs="+",help="the input file(s)")
  parser.add_argument('-p','--padding',type=float,help="the amount to increase the minimum box in each direction, default=2 A",default=2.0)
  parser.add_argument('-s','--spacing',type=float,help="the grid resolution, default=0.25 A",default=0.25)
  args = parser.parse_args()

  if not args.files :
    print "No input files! Nothing to do, so exit."
    quit()

  # Read input files, make the density and write it out in dx-format
  for filename in args.files :
     if filename[-3:].lower() == "pot" :
       xyz,pot = read_potfile(filename)
     elif filename[-4:].lower() == "cube" :
       xyz,pot = read_cubefile(filename)
     print xyz[0,:],pot[0]

     grid,prop = make_density(xyz,pot,padding=args.padding,spacing=args.spacing)

     h,t = os.path.splitext(filename)
     writeDX(grid,prop["min"],args.spacing,h+".dx")


  
