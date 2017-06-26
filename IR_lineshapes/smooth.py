#!/usr/bin/env python

# Author: Richard Bradshaw, R.T.Bradshaw@soton.ac.uk
# Module to smooth curves using univariate spline method

# Requirements: scipy

import scipy.interpolate as inter
import matplotlib.pyplot as plt

class SmoothError(Exception):
    """Exception in curve smoothing"""

class Smoothed():
    """Class to contain methods for curve smoothing from
       scipy.interpolate package."""
    def __init__(self, data):
        """Usage: Smoothed(data). Initialises data as xs and ys attributes.
           Data should be a 2d numpy array of x and y values."""
        if len(data) != 2:
            raise SmoothError("""Your data is formatted incorrectly -  
                                 it should be a 2D array of all x-,
                                 then all y-values""")
        self.xs = data[0]
        self.ys = data[1]

    def do_smooth(self, **kwargs):
        """1D spline fit of data. Usage: self.do_smooth([**kwargs])
           kwargs are passed to scipy.interpolate.UnivariateSpline.

           By default this is a cubic spline fit, with smoothing
           factor s = len(x). In general, larger values of s
           increase smoothing, s = 0 corresponds to interpolating spline"""

        self.spl = inter.UnivariateSpline(self.xs,self.ys,**kwargs)
        self.smoothed_ys = self.spl(self.xs)

    def plot(self, **kwargs):
        """Returns matplotlib figure with original data & smoothed function.
           Usage: self.plots([**kwargs])
           kwargs are passed to matplotlib.Axes.plot for BOTH subplot axes"""
        try:
            fig, ax = plt.subplots()
            ax.plot(self.xs, self.ys, linestyle='dashed', linewidth=1, color='black', **kwargs)
            ax.plot(self.xs, self.smoothed_ys, color='blue', linewidth=2, **kwargs)
            ax.set_xlim((1500,1700))
            ax.set_xlabel('Wavenumber (cm-1)')
            ax.set_ylabel('Intensity (arbitrary units)')
            ax.legend(['Original data','Smoothed data'])
            fig.suptitle("Comparison of original & smoothed data", fontsize=16)
            return fig
        except AttributeError:
            raise FitError("No fit yet performed to plot!")

