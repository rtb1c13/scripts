#!/usr/bin/env python

# Example Matplotlib commands for scatter plots of all datasets

import numpy as np
import matplotlib.pyplot as plt

# Some inputs here
fn = ["New_params.txt","New_polgps.txt","Multifit.txt","Aug_cc.txt","Wat14.txt","OH_scale.txt","Gaff.txt"]

# Separate raw data
# Be sensible and put this all in a proper loop at some point *sigh*
data1 = np.genfromtxt(fn[0])
data2 = np.genfromtxt(fn[1])
data3 = np.genfromtxt(fn[2])
data4 = np.genfromtxt(fn[3])
data5 = np.genfromtxt(fn[4])
data6 = np.genfromtxt(fn[5])
data7 = np.genfromtxt(fn[6])
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
x3 = data3[:,0]
y3 = np.zeros(len(data3))
sdev3 = np.zeros(len(data3))
for i,j in enumerate(data3[:,1:4]):
   y3[i] = np.mean(j)
   sdev3[i] = np.std(j)
x4 = data4[:,0]
y4 = np.zeros(len(data4))
sdev4 = np.zeros(len(data4))
for i,j in enumerate(data4[:,1:4]):
   y4[i] = np.mean(j)
   sdev4[i] = np.std(j)
x5 = data5[:,0]
y5 = np.zeros(len(data5))
sdev5 = np.zeros(len(data5))
for i,j in enumerate(data5[:,1:4]):
   y5[i] = np.mean(j)
   sdev5[i] = np.std(j)
x6 = data6[:,0]
y6 = np.zeros(len(data6))
sdev6 = np.zeros(len(data6))
for i,j in enumerate(data6[:,1:4]):
   y6[i] = np.mean(j)
   sdev6[i] = np.std(j)
x7 = data7[:,0]
y7 = np.zeros(len(data7))
sdev7 = np.zeros(len(data7))
for i,j in enumerate(data7[:,1:4]):
   y7[i] = np.mean(j)
   sdev7[i] = np.std(j)


# Plotting - scatter with error bars
plt.figure(figsize=(8,5))
plt.gca().invert_yaxis()
plt.gca().invert_xaxis()
plt.xlim(-25,5)
plt.ylim(-25,5)
#plt.errorbar(x,y,yerr=sdev,fmt='x',ms=6,mew=1.0,ecolor='black',capthick=2)
#plt.errorbar(x2,y2,yerr=sdev2,fmt='x',ms=6,mew=1.0,ecolor='black',capthick=2)
#plt.errorbar(x3,y3,yerr=sdev3,fmt='x',ms=6,mew=1.0,ecolor='black',capthick=2)
#plt.errorbar(x4,y4,yerr=sdev4,fmt='x',ms=6,mew=1.0,ecolor='black',capthick=2)
#plt.errorbar(x5,y5,yerr=sdev5,fmt='x',ms=6,mew=1.0,ecolor='black',capthick=2)
#plt.errorbar(x6,y6,yerr=sdev6,fmt='x',ms=6,mew=1.0,ecolor='black',capthick=2)
#plt.errorbar(x7,y7,yerr=sdev7,fmt='x',ms=6,mew=1.0,ecolor='black',capthick=2)
#sfig.gca().invert_yaxis()
#sfig.gca().invert_xaxis()
# Regression
m,c = np.polyfit(x,y,1)
m2,c2 = np.polyfit(x2,y2,1)
m3,c3 = np.polyfit(x3,y3,1)
m4,c4 = np.polyfit(x4,y4,1)
m5,c5 = np.polyfit(x5,y5,1)
m6,c6 = np.polyfit(x6,y6,1)
m7,c7 = np.polyfit(x7,y7,1)
xreg = np.arange(-25,6)
plt.plot(xreg,m*xreg+c,'-',lw=1,label="Poltype")
plt.plot(xreg,m2*xreg+c2,'-',lw=1,label="New polgps")
plt.plot(xreg,m3*xreg+c3,'-',lw=1,label="Multi-fit")
plt.plot(xreg,m4*xreg+c4,'-',lw=1,label="Aug-cc-pVTZ")
plt.plot(xreg,m5*xreg+c5,'-',lw=1,label="Wat14")
plt.plot(xreg,m6*xreg+c6,'-',lw=1,label="OH scaled")
plt.plot(xreg,m7*xreg+c7,'-',lw=1,label="GAFF")
plt.plot(xreg,xreg,'--',color='black',lw=1,label="y = x")
# Labels
plt.legend(loc=2,fontsize=12)
plt.suptitle("Comparison of HFE linear regression for all parameter sets",fontweight='bold')
plt.xlabel(r"Experimental $\Delta$G / kcal mol$^{-1}$")
plt.ylabel(r"Predicted $\Delta$G / kcal mol$^{-1}$")
# Show/save
plt.gca().invert_yaxis()
plt.gca().invert_xaxis()


plt.savefig("Fig3.eps",dpi=300,format="eps")
