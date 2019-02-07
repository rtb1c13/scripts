#!/usr/bin/env python

# Script to analyse IR line shapes using FFA methodology of Corcelli et al.
# present in 'spectrum.f'. Reads in a file of pre-enerated fluctuating frequencies

# IR absorption line shape is "Fourier transform of frequency fluctuations
# averaged over finite time interval, t"

import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as igt
import sys


# Import raw frequencies
freqs = np.genfromtxt("/scratch2/rtb1c13/KSI/3VSY/analysis/kubo_RTB_lineshape_test/spectrum_prog/omega_biexp.dat",usecols=(0),max_rows=10**5)
# loadtxt for lower memory usage for large files
#freqs = np.loadtxt("omega_biexp.dat",usecols=(0,))


# Time interval between freqs (currently 'biexp_dynamics' prints at
# 10 step intervals, so dt = 10*0.01
dt = 0.1

# Calculate time averages
# t_len defines number of correlated samples over which integration is performed
avgfrq_wn = np.mean(freqs)
print "avgfreq = %8.3f" % avgfrq_wn
devs = freqs - np.mean(freqs)
t_len = 2**12
integrals = np.zeros(t_len,dtype='complex64')

### Integration/summing loop done with a closure, eg.
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

def integral_maker(freqdevs):
   """Returns adder closed over given array of freqdevs"""
   def integral_adder(integral,step,dt,t_len):
      """Adds integral between step -> step+t_len to current value of integral"""
#      For test data, 2pi multiplier not needed
#      _tmpintegral = igt.cumtrapz(freqdevs[step-1:step-1+t_len]*2*np.pi,dx=dt,initial=0.5*freqdevs[step-1]*dt*2*np.pi)
      _tmpintegral = igt.cumtrapz(freqdevs[step-1:step-1+t_len],dx=dt,initial=0.5*freqdevs[step-1]*dt)
      _tmpintegral[1:] += _tmpintegral[0] # Add initial value ('initial' kwarg doesn't add cumulatively to other cumulative integral vals)
      _tmpintegral = np.exp(_tmpintegral*1j)
      return integral + _tmpintegral
   return integral_adder

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

np.savetxt("integrals_omega_python_2e12_10e5.txt",zip(integrals.real,integrals.imag))

# Do fft with padding equivalent to spectrum.f
fftout = np.fft.fft(integrals,n=(2**4)*t_len)
n = len(fftout)
# Reorder intensities & generate a set of associated frequencies
negfft = fftout[n//2:]
posfft = fftout[:n//2]
frq = np.arange(-1*(len(fftout)/2),len(fftout)/2)
frq = frq/(len(fftout)*dt) # Conversion factor in spectrum.f - solely for test data
reordered = np.concatenate((negfft,posfft))

# Plot & save data
plt.plot(frq,reordered.real)
plt.xlim((-0.1,0.1))
plt.xlabel('Frequency (arbitrary units)')
plt.title("IR lineshape test - biexponential FFA, 2**12 sample length, 10**6 snapshots")
plt.ylabel('Intensity (arbitrary units)')
plt.savefig('python_expected_2e12_10e5.png')
np.savetxt('python_expected_2e12_10e5.dat',zip(frq,reordered.real))

