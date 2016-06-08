#!/usr/bin/env python

# Example Matplotlib commands for scatter plots & overlays

import numpy as np
import matplotlib.pyplot as plt

# Some inputs here
fn = ["New_params.txt","New_polgps.txt"]

# Separate raw data
data = np.genfromtxt(fn[0])
x = data[:,0]
y = np.zeros(len(data))
sdev = np.zeros(len(data))
for i,j in enumerate(data[:,1:4]):
   y[i] = np.mean(j)
   sdev[i] = np.std(j)


# Plotting - scatter with error bars
plt.errorbar(x,y,yerr=sdev,fmt='x',ms=6,mew=1.0,ecolor='black',color='blue',capthick=2)
plt.ylim=(-25,5)
plt.xlim=(-25,5)
# Regression
m,c = np.polyfit(x,y,1)
xreg = np.arange(-25,6)
plt.plot(xreg,m*xreg+c,'-',color='black',lw=1,label="Linear fit")
plt.plot(xreg,xreg,'--',color='black',lw=1,label="y = x")
# Labels
plt.legend(loc=2,fontsize=12)
plt.xlabel(r"Experimental $\Delta$G / kcal mol$^{-1}$")
plt.ylabel(r"Predicted $\Delta$G / kcal mol$^{-1}$")
plt.gca().invert_yaxis()
plt.gca().invert_xaxis()
# Annotations of 25 & 42 
plt.gca().annotate('Solute 25', style='italic', xy=(-4.15,-6.40),xycoords='data',
                   xytext=(-100,20),textcoords='offset points',
                   arrowprops=dict(arrowstyle="->",
                                   connectionstyle="angle3,angleA=-90,angleB=0"),
                   )
plt.gca().annotate('Solute 42', style='italic', xy=(-14.21,-11.20),xycoords='data',
                   xytext=(20,-30),textcoords='offset points',
                   arrowprops=dict(arrowstyle="->",
                                   connectionstyle="angle3,angleA=0,angleB=-90"),
                   )

# Show/save

plt.suptitle("Hydration FE predictions with parameter set 1: Poltype",fontweight='bold')
plt.savefig("Fig1.eps",dpi=300,format="eps")
