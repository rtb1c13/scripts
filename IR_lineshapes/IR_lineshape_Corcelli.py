#!/usr/bin/env python

# Initial script to analyse IR line shapes using methodology of Corcelli et al.

# IR absorption line shape is "Fourier transform of frequency fluctuations
# averaged over finite time interval, t"

import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as igt

#global sinint,cosint,totsin,totcos,dint,nc
ncorr = 500
sinint,cosint=(np.empty(ncorr),np.empty(ncorr))
totsin,totcos=(np.zeros(ncorr),np.zeros(ncorr))
dint=np.empty(ncorr)
nc=0

# Function for doing integral:
def ls0(istep,dfreq,nfreq,tstep):
   global sinint,cosint,totsin,totcos,dint,nc
   ncorr = 500
   clight=2.99792458E10
   dt = tstep*1E-15
   # Above is setup from Corcelli
   if istep == 1:
      for i in range(ncorr):
         sinint[i] = 0.0
         cosint[i] = 0.0
   if istep < ncorr:
      return

   nc += 1
   index0 = (istep%ncorr)+1
   dint=np.zeros(ncorr)
   for j in range(1,ncorr+1):
      index1 = (((j-1)+istep)%ncorr)+1
      index2 = (((j-2)+istep)%ncorr)+1
#      with open("indices.txt",'a') as f:
#         f.write("%s %s %s %s\n" % (index1,index2,dfreq[index1-1],dfreq[index2-1]))
      if j == 1:
         dint[j-1] = 0.5 * dfreq[index1-1] * clight * 2.0*np.pi * dt
      else:
         dint[j-1] = (dint[j-2] + (0.5 * dfreq[index1-1] + 
                      0.5 * dfreq[index2-1])*clight*2.0*np.pi* dt)

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


#      dint = igt.cumtrapz(dfreq,dx=clight*2.0*np.arccos(-1.)*dt,initial=0.5*dfreq[index1-1]*clight*2.0*np.arccos(-1.)*dt)
      dot = 1.
      cosint[j-1] = cosint[j-1] + dot*np.cos(dint[j-1])
      sinint[j-1] = sinint[j-1] + dot*np.sin(dint[j-1])
#      with open("sinint_python",'a') as f:
#         f.write("Step = %s\n" % nc)
#         f.write("%s" % sinint)

   if istep != nfreq:
   #   print nc
      return
   cosint = cosint/nc
   sinint = sinint/nc
   np.savetxt("KSI_out.real",cosint)
   np.savetxt("dint_KSI",dint)
   np.savetxt("KSI_out.imag",sinint)
   return 

### Main
# Import raw frequencies
#data = np.genfromtxt("matlab_Ornstein-Uhlenbeck.txt",usecols=(0),max_rows=10000)
# KSI data below here
data = np.genfromtxt("3VSY_Amoeba_Runs1-3.txt",usecols=(0),max_rows=7500)
freqs = (-0.002915*data**2) + (0.9344*data) + 1735.3
data = freqs
#freqs = freqs*2.998e08 # Convert to Hz

# Define ncorr(?)
ncorr = 500
diffs = np.empty(ncorr)

# Calc average
avgfreq = np.mean(data)
print "avgfreq = %8.3f" % avgfreq

# Calc diffs
for val in range(1,len(data)+1):
   ncorr_idx = ((val-1)%ncorr) +1
   diffs[ncorr_idx-1] = data[val-1] - avgfreq
   print val
   ls0(val,diffs,7500,1.0)
### Fin here?
### so far this recreates the Corcelli code exactly (real & imag outputs)

ncut = 128
dt = 1E-15
clight=2.99792458E10
# Introduce some padding
reals = np.zeros(ncut*72)
imags = np.zeros(ncut*72)

reals[0:ncut] = cosint[0:ncut]
imags[0:ncut] = sinint[0:ncut]

fftvar=np.zeros(len(reals)*2-1,dtype='complex64')
for idx,val in enumerate(reals):
   fftvar[idx] = val - imags[idx]*1j

fftout = np.fft.fft(fftvar)
dom=2*np.pi/(len(fftvar)*dt)
om=np.zeros(len(fftvar)/2)
for i in range(1,(len(fftvar)/2)+1):
   ffttmp = fftout[i-1]
   fftout[i-1] = fftout[(len(fftvar)/2+i)-1]
   fftout[(len(fftvar)/2+i)-1] = ffttmp
   om[i-1] = dom*(i-1-(len(fftvar)/2)/2/np.pi/clight)




# Convert to wavenumbers using quadratic function of Fried et al
#freqs = (-0.002915*data**2) + (0.9344*data) + 1735.3
#freqs = freqs*2.998e08 # Convert to Hz

# Array of times
#Fs = 1e11 # Sampling rate s-1
#Ts = 1.0/Fs # Time interval
#times = np.arange(0,75e-09,Ts) # time vector
#times = times/100 # snapshots to ns
#times = times*1e-09 # ns to s

# Times -> Frequencies
#n = len(freqs)
#k = np.arange(n)
#T = n/Fs
#frq = k/T # Set of frequencies
#frq = frq/2.998e08 # Frq -> Wavenumbers

# Calculate time averages
#timeaves = np.zeros((len(freqs),len(freqs)))
#integrals = np.zeros((len(freqs),len(freqs)))
#expt = np.zeros((len(freqs),len(freqs)),dtype="complex64")
#devs = freqs - np.mean(freqs)
#integrals = None

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
#
# Create fftvar that's average down columns (excluding any 0 values)
#fftvar = np.zeros(len(expt)-1,dtype="complex64") # Final row is all null, hence len(expt)-1
#for i in range(len(expt)-1):
#   fftvar[i] = np.mean(expt[:,i][0:len(expt)-1-i])

# Do fft
#fftout = np.fft.fft(fftvar)
#plt.plot(frq[1:],abs(fftout))
#plt.show()

      
#   for y in range(x,len(freqs)):
#      integrals[x][y] = np.trapz(devs[y:x+y+1],times[y:x+y+1])
#      timeaves[x][y] = 1/times[x] * np.trapz(devs[y:x+y+1],times[y:x+y+1])
#      expt[x][y] = np.exp(1j*integrals[x][y])

# Autocorrelation in devs
#auto = np.correlate(freqs,freqs,mode="full")

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
