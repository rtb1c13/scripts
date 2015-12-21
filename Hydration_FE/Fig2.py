#!/usr/bin/env python

# Example Matplotlib commands for scatter plots & overlays

import numpy as np
import matplotlib.pyplot as plt

# Some inputs here
fn = ["New_params.txt","New_polgps.txt"]
ligindex = [4,5,6,7,8,27,28,29,30,31,32,33,34]

# Separate raw data
data1 = np.genfromtxt(fn[0])
data2 = np.genfromtxt(fn[1])
data1 = data1[ligindex]
data2 = data2[ligindex]
x = data1[:,0]
y = np.zeros(len(data1))
sdev = np.zeros(len(data1))
for i,j in enumerate(data1[:,1:4]):
   y[i] = np.mean(j)
   sdev[i] = np.std(j)
x2 = data2[:,0]
y2 = np.zeros(len(data2))
sdev2 = np.zeros(len(data2))
for i,j in enumerate(data2[:,1:4]):
   y2[i] = np.mean(j)
   sdev2[i] = np.std(j)


# Plotting - scatter with error bars
plt.figure(figsize=(8,5))
plt.xlim(-15,5)
plt.ylim(-15,5)
plt.errorbar(x,y,yerr=sdev,fmt='x',ms=6,mew=1.0,ecolor='black',color='blue',capthick=2)
plt.errorbar(x2,y2,yerr=sdev2,fmt='x',ms=6,mew=1.0,ecolor='black',color='red',capthick=2)
#sfig.gca().invert_yaxis()
#sfig.gca().invert_xaxis()
# Regression
m,c = np.polyfit(x,y,1)
m2,c2 = np.polyfit(x2,y2,1)
xreg = np.arange(-15,6)
plt.plot(xreg,m*xreg+c,'-',color='blue',lw=1,label="Linear fit to Poltype")
plt.plot(xreg,m2*xreg+c2,'-',color='red',lw=1,label="Linear fit to New polgps")
plt.plot(xreg,xreg,'--',color='black',lw=1,label="y = x")
# Labels
plt.legend(loc=2,fontsize=12)
plt.suptitle("Comparison of HFEs for substituted benzenes subgroup",fontweight='bold')
plt.xlabel(r"Experimental $\Delta$G / kcal mol$^{-1}$")
plt.ylabel(r"Predicted $\Delta$G / kcal mol$^{-1}$")
# Show/save
plt.gca().invert_yaxis()
plt.gca().invert_xaxis()


plt.savefig("Fig2.eps",dpi=300,format="eps")
