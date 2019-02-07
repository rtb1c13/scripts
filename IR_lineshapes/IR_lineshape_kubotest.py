#!/usr/bin/env python

# Initial script to analyse IR line shapes using methodology of Corcelli et al.

# IR absorption line shape is "Fourier transform of frequency fluctuations
# averaged over finite time interval, t"

import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as igt
import sys


# Import raw frequencies
freqs = np.genfromtxt("omega_mup_parentdir_trimmed.dat",usecols=(0))

#freqs = data
freqs = freqs*2.99792458e10 # Convert to Hz

# Array of times
Fs = 1e15 # Sampling rate s-1
Ts = 1.0/Fs # Time interval
#times = np.arange(0,10e-6,Ts) # time vector
#times = times/100 # snapshots to ns
#times = times*1e-09 # ns to s

# Times -> Frequencies
n = len(freqs)
k = np.arange(n)
T = n/Fs
frq = k/T # Set of frequencies
frq = frq/2.99792458e10 # Frq -> Wavenumbers
dt = 1e-15

# Calculate time averages
#timeaves = np.zeros((len(freqs),len(freqs)))
#integrals = np.zeros((len(freqs),len(freqs)))
#expt = np.zeros((len(freqs),len(freqs)),dtype="complex64")
avgfrq_wn = np.mean(freqs)/2.99792458e10
print "avgfreq = %8.3f" % avgfrq_wn
devs = freqs - np.mean(freqs)
t_len = 2**15
integrals = np.zeros(t_len,dtype='complex64')

### Should be able to do this loop cleverly with a closure, eg.
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
#      _tmpintegral = np.empty(t_len,dtype='complex64')
      _tmpintegral = igt.cumtrapz(freqdevs[step-1:step-1+t_len]*2*np.pi,dx=dt,initial=0.5*freqdevs[step-1]*dt*2*np.pi)
      _tmpintegral[1:] += _tmpintegral[0] # Add initial value ('initial' kwarg doesn't add cumulatively to other cumulative integral vals)
      _tmpintegral = np.exp(_tmpintegral*1j)
      return integral + _tmpintegral
   return integral_adder

add_step = integral_maker(devs)
for step in range(1,len(devs)+1):
   try:
      integrals = add_step(integrals,step,dt,t_len)
      if step%10000 == 0:
         print "Integration step %s of %s" % (step,len(devs))
   except ValueError:
      print "Stopping integration at step %d, not enough timesteps remain for complete integration (t_len = %d, total steps = %d)" % (step,t_len,len(devs))
      integrals = integrals/(step-1)
      break

print "I broke out"
np.savetxt("cumtrapz_omega_kubotest_trimmed.txt",integrals)
# Create FT variable
#for i,ival in enumerate(freqs):
#   if integrals is None:
#      integrals = igt.cumtrapz(devs,x=times)
#      expt = np.exp(1j*integrals)
#   else:
#      tmpintegrals = igt.cumtrapz(devs[i:],x=times[:-i])
#      tmpexpt = np.exp(1j*tmpintegrals)
#      tmpintegrals.resize((7499,))
#      tmpexpt.resize((7499,))
#      integrals = np.vstack((integrals,tmpintegrals))
#      expt = np.vstack((expt,tmpexpt))

# Create fftvar that's average down columns (excluding any 0 values)
#fftvar = np.zeros(len(expt)-1,dtype="complex64") # Final row is all null, hence len(expt)-1
#for i in range(len(expt)-1):
#   fftvar[i] = np.mean(expt[:,i][0:len(expt)-1-i])

# Do fft
fftout = np.fft.fft(integrals)
n = len(fftout)
negfft = fftout[n//2+1:]
posfft = fftout[:n//2]
negfrq = -1*frq[1:n//2]
posfrq = frq[:n//2]
reordered = np.concatenate((negfft,posfft))
reorderedfrq = np.concatenate((negfrq[::-1],posfrq)) + avgfrq_wn# Reversed negfrq
plt.plot(reorderedfrq,abs(reordered))
plt.xlim((avgfrq_wn-10,avgfrq_wn+10))
plt.xlabel('Wavenumber (cm-1)')
plt.title("IR lineshape test - Kubo parent dir, 2**20 corr")
plt.ylabel('Intensity (arbitrary units)')
plt.show()

      
#   for y in range(x,len(freqs)):
#      integrals[x][y] = np.trapz(devs[y:x+y+1],times[y:x+y+1])
#      timeaves[x][y] = 1/times[x] * np.trapz(devs[y:x+y+1],times[y:x+y+1])
#      expt[x][y] = np.exp(1j*integrals[x][y])

# Autocorrelation in devs
auto = np.correlate(devs,devs,mode="full")
plt.plot(auto)
plt.show()


# Do forward fft
#fftout = np.fft.fft(expt)

# Plot all 4
#fig,ax = plt.subplots(4,1)
#ax[0].plot(times,timeaves)
#ax[1].plot(times,integrals)
#ax[2].plot(times,expt)
#ax[3].plot(times,auto[auto.size/2:])

#plt.show()
#n = len(y) # length of the signal
#k = np.arange(n)
#T = n/Fs
#frq = k/T # two sides frequency range
#frq = frq[range(n/2)] # one side frequency range
#
#Y = np.fft.fft(y)/n # fft computing and normalization
#Y = Y[range(n/2)]
#
#fig, ax = plt.subplots(2, 1)
#ax[0].plot(t,y)
#ax[0].set_xlabel('Time')
#ax[0].set_ylabel('Amplitude')
#ax[1].plot(frq,abs(Y),'r') # plotting the spectrum
#ax[1].set_xlabel('Freq (Hz)')
#ax[1].set_ylabel('|Y(freq)|')
