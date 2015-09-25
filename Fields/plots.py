#!/usr/bin/env python

# Matplotlib histogram of fields

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab

x = np.loadtxt("fields_atN_newZ_run1.txt")
n, bins, patches = plt.hist(x, 50, normed=1, facecolor='green')
plt.show()
print np.mean(x)

#y = mlab.normpdf( bins, mu, sigma)
#l = plt.plot(bins, y, 'r--', linewidth=1)

#From tuts
#plt.xlabel('Smarts')
#plt.ylabel('Probability')
#plt.title(r'$\mathrm{Histogram\ of\ IQ:}\ \mu=100,\ \sigma=15$')
#plt.axis([40, 160, 0, 0.03])
#plt.grid(True)

#plt.show()
