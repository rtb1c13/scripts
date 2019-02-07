#!/usr/bin/env python

# Initial script to analyse IR line shapes using methodology of Corcelli et al.

# IR absorption line shape is "Fourier transform of frequency fluctuations
# averaged over finite time interval, t"

import numpy as np
### For remote running/figure generation without xwindow
import matplotlib as mpl
mpl.use('Agg')
###
import matplotlib.pyplot as plt
import scipy.integrate as igt
import mdtraj as md
import sys,argparse
from subprocess import call

### Argparser ###
def parse():
   parser = argparse.ArgumentParser()
   parser.add_argument("-a","--atoms",help="Atom numbers for analysis (Max 4)",nargs='+',type=int,required=True)
   parser.add_argument("-p","--polarise",help="Polarisabilities for atoms 1 & 2. Defaults to carbonyl C and O.",nargs=2,type=float,default=[1.3340,0.8370])
   parser.add_argument("-t","--traj",help="Coordinate file (Tinker uindxyz, renumbered correctly)",type=str,required=True)
   parser.add_argument("-dt","--diptraj",help="Dipoles file (Tinker uind, numbering unimportant)",type=str)
   parser.add_argument("-g","--gas",help="Flag to calculate gas phase dipoles for environment field calculation. By default this is switched off",dest='gas',action='store_true')
   parser.add_argument("-pre","--prefix",help="Prefix for output filenames. Default 'KSI'", type=str, default='KSI')
   parser.set_defaults(gas=False)

   if len(sys.argv)==1:
      parser.print_help()
      sys.exit(1)

   args = parser.parse_args()
   return args

### Class for atoms ###
class Atom:
   """Class to define properties for atom numbers
      read in from command line arguments"""

   def __init__(self,num):
      """Defines index for atom num"""
      self.idx = num-1

   def polarise(self,num):
      """Defines polarisability for atom num"""
      listidx = args.atoms.index(num)
      if listidx > 1:
         return
      else:
         self.alpha = args.polarise[listidx]

### Class for trajectory ###
class Anal_traj:
   """Class to define all the desired analyses for a
      Tinker trajectory read in using MDTraj."""   

   def __init__(self,arcname):
      """Reads in Tinker trajectory arcname."""
      self.traj = md.load_arc(arcname)

   def getcoords(self,atmlst):
      """Gets coordinates of atom indices
         atm1 & atm2 given in a list. Units in Angstrom."""
      self.coords1 = np.zeros((self.traj.n_frames,3))
      self.coords2 = np.zeros((self.traj.n_frames,3))
      for i in range(0,self.traj.n_frames):
         self.coords1[i] = self.traj[i].xyz[0][atmlst[0].idx] *10
         self.coords2[i] = self.traj[i].xyz[0][atmlst[1].idx] *10

   def getdipls(self,atmlst):
      """Gets dipoles of atom indices
         atm1 & atm2 given in a list. Units in Debye."""
      self.dipls1 = np.zeros((self.traj.n_frames,3))
      self.dipls2 = np.zeros((self.traj.n_frames,3))
      self.field1 = np.zeros((self.traj.n_frames,3))
      self.field2 = np.zeros((self.traj.n_frames,3))
      for i in range(0,self.traj.n_frames):
         self.dipls1[i] = self.traj[i].xyz[0][atmlst[0].idx] *10
         self.dipls2[i] = self.traj[i].xyz[0][atmlst[1].idx] *10
         self.field1[i] = self.dipls1[i]/atmlst[0].alpha
         self.field2[i] = self.dipls2[i]/atmlst[1].alpha

   def vectors(self):
      """Defines the interatomic vector, length and
         unit vector along that path, assuming coordinates
         for atoms are already defined. Returns three arrays:
          1) X/Y/Z vector components
          2) Interatomic vector length
          3) X/Y/Z unit vector components."""
      self.vec = np.zeros((self.traj.n_frames,3))
      self.leng = np.zeros(self.traj.n_frames)
      self.unitvec = np.zeros((self.traj.n_frames,3))
      for i in range(0,self.traj.n_frames):
         self.vec[i] = self.coords2[i] - self.coords1[i]
         self.leng[i] = np.linalg.norm(self.vec[i])
         self.unitvec[i] = self.vec[i] / self.leng[i]

def calcfield(trajobj,dipobj):
   """Calculates average of field projections in MV/cm 
      along unit vector at desired coordinates"""
   fieldproj1 = np.zeros(trajobj.traj.n_frames)
   fieldproj2 = np.zeros(trajobj.traj.n_frames)
   for i in range(0,trajobj.traj.n_frames):
      fieldproj1[i] = np.dot(dipobj.field1[i], trajobj.unitvec[i])
      fieldproj1[i] *= 299.79
      fieldproj2[i] = np.dot(dipobj.field2[i], trajobj.unitvec[i])
      fieldproj2[i] *= 299.79
   avefield = (fieldproj1+fieldproj2)/2
   return avefield

def analyze_dipl_detailed(trajobj,atmlst):
   """Assigns induced dipole components for defined atoms
      using Tinker analyze_dipdebug. This is a modified
      version of analyze that allows the printout of induced
      dipoles only (for speed). Note that the keyfile
      'testgas.key' must be prepared elsewhere, and should
      print out dipoles to 6dp."""

   def _filter_line(fileobj):
      """Generator to filter out lines beginning with
         digit for np.loadtxt"""
      for line in fileobj:
         try:
            if line.split()[0][0].isdigit():
               yield line
         except IndexError:           # Empty line
            continue
 
   atm1=atmlst[0].idx
   atm2=atmlst[1].idx
   alpha1=atmlst[0].alpha
   alpha2=atmlst[1].alpha

   # Call analyze on trajectory
   call('/local/scratch/rtb1c13/Git/tinker/bin/analyze_dipdebug %s -k testgas.key u &> analout.txt ' % args.traj,shell=True)
   call('grep -A%s "Induced Dipole Moments (Debyes) :" analout.txt > dipls.txt && rm analout.txt' % (trajobj.traj.n_atoms+3),shell=True)
   with open('dipls.txt','r') as f:
      dipls = np.loadtxt(_filter_line(f),usecols=(1,2,3))

   # Convert dipoles to fields
   dipls = dipls.reshape((trajobj.traj.n_frames,trajobj.traj.n_atoms,3))
   field = np.zeros((2,trajobj.traj.n_frames,3))
   for i in range(0,trajobj.traj.n_frames):
      field[0][i][0] = dipls[i][atm1][0]
      field[0][i][1] = dipls[i][atm1][1]
      field[0][i][2] = dipls[i][atm1][2]
      field[1][i][0] = dipls[i][atm2][0]
      field[1][i][1] = dipls[i][atm2][1]
      field[1][i][2] = dipls[i][atm2][2]
   field[0] = field[0] / alpha1
   field[1] = field[1] / alpha2
   return field

# Funcs for summing integrals
def integral_maker(freqdevs):
   """Returns adder closed over given array of freqdevs"""
   def integral_adder(integral,step,dt,t_len):
      """Adds integral between step -> step+t_len to current value of integral"""
      _tmpintegral = igt.cumtrapz(freqdevs[step-1:step-1+t_len]*2*np.pi,dx=dt,initial=0.5*freqdevs[step-1]*dt*2*np.pi)
#     No *2*pi
#      _tmpintegral = igt.cumtrapz(freqdevs[step-1:step-1+t_len],dx=dt,initial=0.5*freqdevs[step-1]*dt)
      _tmpintegral[1:] += _tmpintegral[0] # Add initial value ('initial' kwarg doesn't add cumulatively to other cumulative integral vals)
      _tmpintegral = np.exp(_tmpintegral*1j)
      return integral + _tmpintegral
   return integral_adder

# Main...

global args
args = parse()

coordsname = args.traj
dipname = args.diptraj
uindxyz = Anal_traj(coordsname)
uind = Anal_traj(dipname)
atmlst=[]
for i in args.atoms:
   j = Atom(i)
   j.polarise(i)
   atmlst.append(j)
print "Projecting along internuclear vector %d to %d" % (args.atoms[0],args.atoms[1])
uindxyz.getcoords(atmlst)
uindxyz.vectors()
uind.getdipls(atmlst)
if args.gas:
   gasfield = analyze_dipl_detailed(uindxyz,atmlst)
   uind.field1 -= gasfield[0]
   uind.field2 -= gasfield[1]
data = calcfield(uindxyz,uind)
np.savetxt("%s_fields.txt" % args.prefix,data,fmt="%8.3f")


# Read in induced dipoles & ligand coordinates
# For KSI, only 1 polgroup means only complex dipoles are needed (lig dipoles = 0)


# Convert to wavenumbers using quadratic function of Fried et al
# from supp text S2, fig S4F
#freqs = (-0.00125*data**2) + (0.610*data) + 1688.2
# SOLVENTS: Convert to wavenumbers using linear (R2 = 0.98) function of Fried & Wang
# from Fig 1a, JPCB 2013
#freqs = (0.484*data) + 1703.6
# OR: Convert to wavenumbers using linear (R2 = 0.96) function fitted to my field data for Acetophenone
#freqs = (0.4686*data) + 1701.1
# AMOEBA: Convert to wavenumbers using quadratic function fitted to AMOEBA 19-NT solvent fields
freqs = (-0.00131567*data**2) + 0.5496*data + 1699.2
#freqs = data
np.savetxt("%s_wavenums.txt" % args.prefix,freqs,fmt="%10.6f")
freqs = freqs*2.99792458e10 # Convert to Hz


# Times -> Frequencies
dt = 2e-15 # Timestep
t_len = 2**12 # 4096 steps, 8.192 ps
Fs = 1./dt # Sampling rate s-1
n = t_len*(2**4) # Samples in FFT signal inc. padding
k = np.arange(n)
T = n/Fs # Time length of fft signal in s
frq = k/T # Set of frequencies
frq = frq/2.99792458e10 # Frq -> Wavenumbers

# Calculate time averages & deviations
avgfrq_wn = np.mean(freqs)/2.99792458e10
print "avgfreq = %8.3f" % avgfrq_wn
devs = freqs - np.mean(freqs)
integrals = np.zeros(t_len,dtype='complex64')

### Should be able to do this loop with a closure, eg.
#      def make_dfreq_adder(dfreq):
#         num_to_add = igt.cumtrapz(dfreq)
#         print num_to_add
#         def integral_adder(integral):
#            return integral + num_to_add
#         return integral_adder
### This creates an integral adder with a different value of dfreq each time
### Then e.g.:
#      int_sum=0
#      for step in steps:
#      adder = make_dfreq_adder(diffs[step:steps+ncorr])
#      int_sum = adder(int_sum)
#      print "int_sum = ",int_sum


add_step = integral_maker(devs)
for step in range(1,len(devs)+1):
   try:
      integrals = add_step(integrals,step,dt,t_len)
      if step%10000 == 0:
         print "Integration for step %s of %s" % (step,len(devs))
   except ValueError:
      print "Stopping averaging at step %d, not enough timesteps remain for complete integration (t_len = %d, total steps = %d)" % (step,t_len,len(devs))
      integrals = integrals/(step-1)
      break

np.savetxt("integrals_%s.txt" % args.prefix,zip(integrals.real,integrals.imag))

# Do fft with 16*padding
fftout = np.fft.fft(integrals,n=(2**4)*t_len)
n = len(fftout)
# Reorder intensities & frequencies
negfft = fftout[n//2:]
posfft = fftout[:n//2]
negfrq = -1*frq[1:n//2+1] # Same as arange(-len(fftout)/2,len(fftout)/2)
posfrq = frq[:n//2]
reordered = np.concatenate((negfft,posfft))
reorderedfrq = np.concatenate((negfrq[::-1],posfrq)) + avgfrq_wn# Reversed negfrq
plt.plot(reorderedfrq,reordered.real)
#plt.plot(reorderedfrq,abs(reordered))
plt.xlim((avgfrq_wn-100,avgfrq_wn+100))
plt.xlabel('Wavenumber (cm-1)')
plt.title("19-NT in KSI fields, %s" % args.prefix)
plt.ylabel('Intensity (arbitrary units)')
plt.savefig('%s_lineshape.png' % args.prefix, dpi=300)
np.savetxt('%s_lineshape.dat' % args.prefix,zip(reorderedfrq,reordered.real))
#np.savetxt('10ns_lineshape_abs.dat',zip(reorderedfrq,abs(reordered)))


      
