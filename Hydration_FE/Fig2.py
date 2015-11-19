#!/usr/bin/env python

# Example Matplotlib commands for scatter plots & overlays

import numpy as np
import matplotlib.pyplot as plt

# Some inputs here
fn = ["New_params.txt","New_polgps.txt"]
ligs = [

# Separate raw data
data1 = np.genfromtxt(fn[0])
data2 = np.genfromtxt(fn[0])
data1 = 
x = data[:,0]
y = np.zeros(len(data))
sdev = np.zeros(len(data))
for i,j in enumerate(data[:,1:4]):
   y[i] = np.mean(j)
   sdev[i] = np.std(j)
print y


# Plotting - scatter with error bars
fig,sfigs=plt.subplots(nrows=2,ncols=1,sharex=True,sharey=True,figsize=(8,10))
fig.gca().invert_yaxis()
fig.gca().invert_xaxis()
sfig=sfigs[0]
sfig.errorbar(x,y,yerr=sdev,fmt='x',ms=6,mew=1.0,ecolor='black',color='blue',capthick=2)
sfig.ylim=(-25,5)
#sfig.gca().invert_yaxis()
#sfig.gca().invert_xaxis()
# Regression
m,c = np.polyfit(x,y,1)
xreg = np.arange(-25,6)
sfig.plot(xreg,m*xreg+c,'-',color='black',lw=1,label="Linear fit")
sfig.plot(xreg,xreg,'--',color='black',lw=1,label="y = x")
# Labels
sfig.legend(loc=2)
sfig.set_title("HFEs of parameter set 1: Poltype",fontweight='bold')
sfig.set_xlabel(r"Experimental $\Delta$G / kcal mol$^{-1}$")
sfig.set_ylabel(r"Predicted $\Delta$G / kcal mol$^{-1}$")

#Sfig 2
data = np.genfromtxt(fn[1])
x = data[:,0]
y = np.zeros(len(data))
sdev = np.zeros(len(data))
for i,j in enumerate(data[:,1:4]):
   y[i] = np.mean(j)
   sdev[i] = np.std(j)
print y
# Scatter with error bars
sfig=sfigs[1]
sfig.errorbar(x,y,yerr=sdev,fmt='x',ms=6,mew=1.0,ecolor='black',color='blue',capthick=2)
sfig.ylim=(-25,5)
fig.gca().invert_yaxis()
fig.gca().invert_xaxis()
#sfig.gca().invert_yaxis()
#sfig.gca().invert_xaxis()
# Regression
m,c = np.polyfit(x,y,1)
xreg = np.arange(-25,6)
sfig.plot(xreg,m*xreg+c,'-',color='black',lw=1,label="Linear fit")
sfig.plot(xreg,xreg,'--',color='black',lw=1,label="y = x")
# Labels
sfig.legend(loc=2)
sfig.set_title("HFEs of parameter set 2: New polgps",fontweight='bold')
sfig.set_xlabel(r"Experimental $\Delta$G / kcal mol$^{-1}$")
sfig.set_ylabel(r"Predicted $\Delta$G / kcal mol$^{-1}$")
# Show/save


plt.savefig("Fig1.eps",dpi=300,format="eps")
