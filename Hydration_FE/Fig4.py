#!/usr/bin/env python

# Example Matplotlib commands for scatter plots & overlays

import numpy as np
import matplotlib.pyplot as plt

# Some inputs here
fn = ["New_params.txt","New_polgps.txt","Multifit.txt","Aug_cc.txt","Wat14.txt","OH_scale.txt","Gaff.txt"]

# Separate raw data
data = np.genfromtxt(fn[0])
x = data[:,0]
y = np.zeros(len(data))
sdev = np.zeros(len(data))
for i,j in enumerate(data[:,1:4]):
   y[i] = np.mean(j)
   sdev[i] = np.std(j)


# Plotting - scatter with error bars
fig,sfigs=plt.subplots(nrows=4,ncols=2,sharex=True,sharey=True,figsize=(12,16))
#sfigs[0,0].invert_yaxis()
#sfigs[0,0].invert_xaxis()
sfig=sfigs[(0,0)]
sfig.errorbar(x,y,yerr=sdev,fmt='x',ms=6,mew=1.0,ecolor='black',color='blue',capthick=2)
sfig.ylim=(-35,5)
# Regression
m,c = np.polyfit(x,y,1)
xreg = np.arange(-35,6)
sfig.plot(xreg,m*xreg+c,'-',color='black',lw=1,label="Linear fit")
sfig.plot(xreg,xreg,'--',color='black',lw=1,label="y = x")
# Labels
sfig.legend(loc=2,fontsize=12)
sfig.set_title("Parameter set 1: Poltype",fontweight='bold')
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
# Scatter with error bars
sfig=sfigs[(0,1)]
sfig.errorbar(x,y,yerr=sdev,fmt='x',ms=6,mew=1.0,ecolor='black',color='blue',capthick=2)
sfig.ylim=(-35,5)
# Regression
m,c = np.polyfit(x,y,1)
xreg = np.arange(-35,6)
sfig.plot(xreg,m*xreg+c,'-',color='black',lw=1,label="Linear fit")
sfig.plot(xreg,xreg,'--',color='black',lw=1,label="y = x")
# Labels
sfig.legend(loc=2,fontsize=12)
sfig.set_title("Parameter set 2: New polgps",fontweight='bold')
sfig.set_xlabel(r"Experimental $\Delta$G / kcal mol$^{-1}$")
sfig.set_ylabel(r"Predicted $\Delta$G / kcal mol$^{-1}$")

#Sfig 3
data = np.genfromtxt(fn[2])
x = data[:,0]
y = np.zeros(len(data))
sdev = np.zeros(len(data))
for i,j in enumerate(data[:,1:4]):
   y[i] = np.mean(j)
   sdev[i] = np.std(j)
# Scatter with error bars
sfig=sfigs[(1,0)]
sfig.errorbar(x,y,yerr=sdev,fmt='x',ms=6,mew=1.0,ecolor='black',color='blue',capthick=2)
sfig.ylim=(-35,5)
# Regression
m,c = np.polyfit(x,y,1)
xreg = np.arange(-35,6)
sfig.plot(xreg,m*xreg+c,'-',color='black',lw=1,label="Linear fit")
sfig.plot(xreg,xreg,'--',color='black',lw=1,label="y = x")
# Labels
sfig.legend(loc=2,fontsize=12)
sfig.set_title("Parameter set 3: Multi-fit",fontweight='bold')
sfig.set_xlabel(r"Experimental $\Delta$G / kcal mol$^{-1}$")
sfig.set_ylabel(r"Predicted $\Delta$G / kcal mol$^{-1}$")

#Sfig 4
data = np.genfromtxt(fn[3])
x = data[:,0]
y = np.zeros(len(data))
sdev = np.zeros(len(data))
for i,j in enumerate(data[:,1:4]):
   y[i] = np.mean(j)
   sdev[i] = np.std(j)
# Scatter with error bars
sfig=sfigs[(1,1)]
sfig.errorbar(x,y,yerr=sdev,fmt='x',ms=6,mew=1.0,ecolor='black',color='blue',capthick=2)
sfig.ylim=(-35,5)
# Regression
m,c = np.polyfit(x,y,1)
xreg = np.arange(-35,6)
sfig.plot(xreg,m*xreg+c,'-',color='black',lw=1,label="Linear fit")
sfig.plot(xreg,xreg,'--',color='black',lw=1,label="y = x")
# Labels
sfig.legend(loc=2,fontsize=12)
sfig.set_title("Parameter set 4: Aug-cc-PVTZ",fontweight='bold')
sfig.set_xlabel(r"Experimental $\Delta$G / kcal mol$^{-1}$")
sfig.set_ylabel(r"Predicted $\Delta$G / kcal mol$^{-1}$")
#Sfig 5
data = np.genfromtxt(fn[4])
x = data[:,0]
y = np.zeros(len(data))
sdev = np.zeros(len(data))
for i,j in enumerate(data[:,1:4]):
   y[i] = np.mean(j)
   sdev[i] = np.std(j)
# Scatter with error bars
sfig=sfigs[(2,0)]
sfig.errorbar(x,y,yerr=sdev,fmt='x',ms=6,mew=1.0,ecolor='black',color='blue',capthick=2)
sfig.ylim=(-35,5)
# Regression
m,c = np.polyfit(x,y,1)
xreg = np.arange(-35,6)
sfig.plot(xreg,m*xreg+c,'-',color='black',lw=1,label="Linear fit")
sfig.plot(xreg,xreg,'--',color='black',lw=1,label="y = x")
# Labels
sfig.legend(loc=2,fontsize=12)
sfig.set_title("Parameter set 5: Wat14",fontweight='bold')
sfig.set_xlabel(r"Experimental $\Delta$G / kcal mol$^{-1}$")
sfig.set_ylabel(r"Predicted $\Delta$G / kcal mol$^{-1}$")

#Sfig 6
data = np.genfromtxt(fn[5])
x = data[:,0]
y = np.zeros(len(data))
sdev = np.zeros(len(data))
for i,j in enumerate(data[:,1:4]):
   y[i] = np.mean(j)
   sdev[i] = np.std(j)
# Scatter with error bars
sfig=sfigs[(2,1)]
sfig.errorbar(x,y,yerr=sdev,fmt='x',ms=6,mew=1.0,ecolor='black',color='blue',capthick=2)
sfig.ylim=(-35,5)
# Regression
m,c = np.polyfit(x,y,1)
xreg = np.arange(-35,6)
sfig.plot(xreg,m*xreg+c,'-',color='black',lw=1,label="Linear fit")
sfig.plot(xreg,xreg,'--',color='black',lw=1,label="y = x")
# Labels
sfig.legend(loc=2,fontsize=12)
sfig.set_title("Parameter set 6: OH scaled",fontweight='bold')
sfig.set_xlabel(r"Experimental $\Delta$G / kcal mol$^{-1}$")
sfig.set_ylabel(r"Predicted $\Delta$G / kcal mol$^{-1}$")
plt.setp(sfig.get_xticklabels(),visible=True)
#Sfig 7
data = np.genfromtxt(fn[6])
x = data[:,0]
y = np.zeros(len(data))
sdev = np.zeros(len(data))
for i,j in enumerate(data[:,1:4]):
   y[i] = np.mean(j)
   sdev[i] = np.std(j)
# Scatter with error bars
sfig=sfigs[(3,0)]
sfig.errorbar(x,y,yerr=sdev,fmt='x',ms=6,mew=1.0,ecolor='black',color='blue',capthick=2)
sfig.ylim=(-35,5)
# Regression
m,c = np.polyfit(x,y,1)
xreg = np.arange(-35,6)
sfig.plot(xreg,m*xreg+c,'-',color='black',lw=1,label="Linear fit")
sfig.plot(xreg,xreg,'--',color='black',lw=1,label="y = x")
# Labels
sfig.legend(loc=2,fontsize=12)
sfig.set_title("Parameter set 7: GAFF",fontweight='bold')
sfig.set_xlabel(r"Experimental $\Delta$G / kcal mol$^{-1}$")
sfig.set_ylabel(r"Predicted $\Delta$G / kcal mol$^{-1}$")

# Show/save
sfigs[-1,-1].axis('off')
plt.xlim(-35,5)
plt.ylim(-35,5)
plt.gca().invert_xaxis()
plt.gca().invert_yaxis()

plt.savefig("Fig4.eps",dpi=300,format="eps")
